# ============================================================
# Paper: Taxa de coleta de resíduos sólidos
# Estimação - Full Sample vs. Only Treated (Janela Balanceada -3 a +3)
# Autoria: Cosmo Hugo da Silva (JEEM Pipeline)
# ============================================================


# ============================================================
# 1. CARREGAR PACOTES
# ============================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  did,
  showtext,
  gt,
  WeightIt,
  kableExtra
)


# ============================================================
# 2. RENDERIZAÇÃO DE FONTES
# ============================================================

font_family <- "STIX Two Text"
font_add_google(font_family)
showtext_auto()


# ============================================================
# 3. PROPENSITY SCORE / ENTROPY BALANCING
# ============================================================

dados_ps <- dados_saneamento |>
  filter(ano == 2009) |>
  mutate(
    tratado = if_else(primeiro_tratamento > 0, 1, 0)
  )

ps_weight <- weightit(
  tratado ~
    tx_pop_acesso_agua +
    taxa_cob_imun +
    pib_pc +
    tx_inter_feco_oral +
    tx_inter_inseto_vetor,
  data = dados_ps,
  method = "ebal",
  estimand = "ATT"
)

summary(ps_weight)

# Adicionar pesos à base
dados_ps <- dados_ps |>
  mutate(ps_weight = ps_weight$weights) |>
  select(cod_mun, ps_weight)


# ============================================================
# 4. ADICIONAR PESOS À BASE PRINCIPAL
# ============================================================

dados_saneamento <- dados_saneamento |>
  left_join(dados_ps, by = "cod_mun")


# ============================================================
# 5. CRIAR PASTA DE RESULTADOS ESCALONADOS
# ============================================================

pasta_destino <- "resultados_principais/resultados_escalonados/"

if (!dir.exists(pasta_destino)) {
  dir.create(pasta_destino, recursive = TRUE)
}


# ============================================================
# 6. COVARIÁVEIS
# ============================================================

covariates <- c(
  "tx_pop_acesso_agua",
  "taxa_cob_imun",
  "pib_pc"
)

xformula_str <- paste("~", paste(covariates, collapse = " + "))
xformula <- as.formula(xformula_str)


# ============================================================
# 7. FUNÇÃO PARA SALVAR TABELA EM TXT
# ============================================================

salvar_tabela_txt <- function(resultado, outcome, path = pasta_destino) {
  
  nome_arquivo <- paste0(path, outcome, ".txt")
  
  tabela <- resultado |>
    select(amostra, ATT, SE, N, baseline) |>
    mutate(
      ATT = round(ATT, 4),
      SE = round(SE, 4),
      baseline = round(baseline, 4)
    )
  
  texto <- paste0(
    "============================================================\n",
    "Outcome: ", outcome, "\n",
    "Janela de Estudo de Eventos: -3 a +3 (Sample Balanceada)\n",
    "============================================================\n\n",
    "Especificação          ATT         SE          N        Baseline\n",
    "------------------------------------------------------------\n"
  )
  
  for (i in seq_len(nrow(tabela))) {
    linha <- sprintf(
      "%-22s %-11.4f %-11.4f %-8d %-11.4f\n",
      tabela$amostra[i],
      tabela$ATT[i],
      tabela$SE[i],
      tabela$N[i],
      tabela$baseline[i]
    )
    texto <- paste0(texto, linha)
  }
  
  texto <- paste0(
    texto,
    "\n",
    "============================================================\n",
    "Covariáveis:\n",
    "tx_pop_acesso_agua + taxa_cob_imun + pib_pc\n",
    "Método: Callaway & Sant'Anna (DR)\n",
    "Grupo de controle: Not-yet-treated\n",
    "Pesos: ps_weight (Entropy Balancing)\n",
    "============================================================\n"
  )
  
  writeLines(texto, nome_arquivo, useBytes = TRUE)
  return(nome_arquivo)
}


# ============================================================
# 8. FUNÇÃO DE ESTIMAÇÃO DID (JANELA BALANCEADA -3 A +3)
# ============================================================

estimacao_did <- function(
    outcome,
    data = dados_saneamento,
    xformla = xformula,
    tname = "ano",
    idname = "cod_mun",
    gname = "primeiro_tratamento",
    weightsname = "ps_weight",
    clustervars = "cod_mun",
    est_method = "dr",
    control_group = "notyettreated",
    base_period = "varying",
    min_e = -3,
    max_e = 3,
    path = pasta_destino
) {
  
  y_sym <- rlang::sym(outcome)
  
  # ----------------------------------------------------------
  # FUNÇÃO INTERNA COMPATÍVEL COM GRUPOS DINÂMICOS
  # ----------------------------------------------------------
  rodar_did <- function(base, nome_amostra, var_gname) {
    
    cat("\n------------------------------------------\n")
    cat("Outcome:", outcome, "\n")
    cat("Amostra:", nome_amostra, "\n")
    cat("Janela: -3 a +3\n")
    cat("------------------------------------------\n")
    
    modelo_cs <- did::att_gt(
      yname = outcome,
      tname = tname,
      idname = idname,
      gname = var_gname,
      xformla = xformla,
      panel = FALSE,
      allow_unbalanced_panel = FALSE,
      control_group = control_group,
      weightsname = weightsname,
      clustervars = clustervars,
      est_method = est_method,
      base_period = base_period,
      data = base
    )
    
    output_cs <- did::aggte(
      modelo_cs,
      type = "dynamic",
      na.rm = TRUE,
      min_e = min_e,
      max_e = max_e
    )
    
    crit_val <- output_cs$crit.val.egt
    
    # Cálculo do Baseline
    baseline <- base |>
      mutate(
        tratado = if_else(.data[[var_gname]] > 0, 1, 0),
        ano_relativo = .data[[var_gname]] - .data[[tname]]
      ) |>
      filter(ano_relativo == -1, tratado == 1) |>
      pull(!!y_sym) |>
      mean(na.rm = TRUE)
    
    # Event Study Data Frame
    df <- data.frame(
      t_label = output_cs$egt,
      coeficientes = output_cs$att.egt,
      se = output_cs$se.egt
    ) |>
      tidyr::drop_na(se) |>
      filter(t_label >= -3 & t_label <= 3) |>
      mutate(
        ymin = coeficientes - crit_val * se,
        ymax = coeficientes + crit_val * se,
        t_label_str = as.character(t_label),
        amostra = nome_amostra
      ) |>
      mutate(
        t_label = factor(
          t_label_str,
          levels = c("-3", "-2", "-1", "0", "1", "2", "3")
        )
      )
    
    resultados <- tibble::tibble(
      amostra = nome_amostra,
      outcome = outcome,
      ATT = output_cs$overall.att,
      SE = output_cs$overall.se,
      N = as.integer(output_cs$DIDparams$n),
      baseline = baseline
    )
    
    return(list(df = df, res = resultados))
  }
  
  # ----------------------------------------------------------
  # 1. FULL SAMPLE
  # ----------------------------------------------------------
  full <- rodar_did(data, "Full sample", var_gname = gname)
  
  # ----------------------------------------------------------
  # 2. ONLY TREATED (AJUSTADO: COORTES BALANCEADAS 2012-2015)
  # ----------------------------------------------------------
  dados_so_tratados <- data |>
    # Converte a coorte de 2018 para Grupo 0 (controle)
    mutate(
      primeiro_tratamento_only = if_else(
        .data[[gname]] == 2018,
        0L,
        as.integer(.data[[gname]])
      )
    ) |>
    # Mantém apenas o Grupo 0 + coortes balanceadas conforme pedido do revisor (2012 a 2015)
    filter(
      primeiro_tratamento_only == 0L | 
        (primeiro_tratamento_only >= 2012 & primeiro_tratamento_only <= 2015)
    )
  
  treated <- rodar_did(dados_so_tratados, "Only treated", var_gname = "primeiro_tratamento_only")
  
  # ----------------------------------------------------------
  # 3. SALVAR TABELA TXT
  # ----------------------------------------------------------
  resultados_outcome <- bind_rows(full$res, treated$res)
  salvar_tabela_txt(resultados_outcome, outcome, path)
  
  # ----------------------------------------------------------
  # 4. PLOT EVENT STUDY COMBINADO
  # ----------------------------------------------------------
  df_plot <- bind_rows(full$df, treated$df) |>
    mutate(
      x_num = as.numeric(t_label),
      x_plot = ifelse(amostra == "Full sample", x_num - 0.1, x_num + 0.1)
    )
  
  plot <- ggplot(
    df_plot,
    aes(x = x_plot, y = coeficientes, color = amostra, shape = amostra)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
    geom_point(size = 3.5) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
    scale_x_continuous(
      breaks = 1:7,
      labels = levels(df_plot$t_label)
    ) +
    scale_color_manual(
      values = c("Full sample" = "#006D77", "Only treated" = "#E29578")
    ) +
    scale_shape_manual(
      values = c("Full sample" = 16, "Only treated" = 15)
    ) +
    xlab("Relative time to treatment") +
    ylab("Coefficients") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 25),
      legend.text = element_text(size = 27)
    )
  
  ggsave(
    filename = paste0(outcome, "_full_vs_treated_balanced.pdf"),
    plot = plot,
    device = "pdf",
    path = path,
    width = 14,
    height = 8.5
  )
  
  print(plot)
  return(resultados_outcome)
}


# ============================================================
# 9. OUTCOMES - RESÍDUOS SÓLIDOS
# ============================================================

outcomes <- c(
  "d_plano_gestao_residuos",
  "existe_lixao",
  "d_coleta_seletiva",
  "tx_pop_resid_solidos",
  "tx_pop_coleta_diaria",
  "tx_pop_coleta_2_3_semana"
)

resultados1 <- purrr::map_dfr(
  outcomes,
  ~ estimacao_did(.x)
)


# ============================================================
# 10. OUTCOMES - SAÚDE
# ============================================================

outcomes_saude <- c(
  "tx_inter_feco_oral",
  "tx_inter_contato_agua",
  "tx_inter_higiente"
)

resultados2 <- purrr::map_dfr(
  outcomes_saude,
  ~ estimacao_did(.x)
)


# ============================================================
# 11. RESULTADOS FINAIS
# ============================================================

tabela_final <- bind_rows(resultados1, resultados2)
tabela_final
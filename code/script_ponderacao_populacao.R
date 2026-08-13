# Paper: Taxa de coleta de resíduos sólidos
# Estimação - Only Treated ponderado pela população
# ============================================================


# ============================================================
# 1. CARREGAR PACOTES
# ============================================================

pacman::p_load(
  tidyverse,
  did,
  showtext,
  gt,
  kableExtra
)


# ============================================================
# 2. RENDERIZAÇÃO DE FONTES
# ============================================================

font_family <- "STIX Two Text"

font_add_google(font_family)

showtext_auto()


# ============================================================
# 3. CRIAR PASTA DE RESULTADOS
# ============================================================

if (!dir.exists("resultados_ponderados_pela_populacao")) {
  dir.create(
    "resultados_ponderados_pela_populacao",
    recursive = TRUE
  )
}


# ============================================================
# 4. COVARIÁVEIS
# ============================================================

covariates <- c(
  "tx_pop_acesso_agua",
  "taxa_cob_imun",
  "pib_pc"
)

xformula_str <- paste(
  "~",
  paste(
    covariates,
    collapse = " + "
  )
)

xformula <- as.formula(
  xformula_str
)


# ============================================================
# 5. FUNÇÃO PARA SALVAR TABELA EM TXT
# ============================================================

salvar_tabela_txt <- function(
    resultado,
    outcome,
    path = "resultados_ponderados_pela_populacao/"
) {
  
  nome_arquivo <- paste0(
    path,
    outcome,
    ".txt"
  )
  
  
  # ----------------------------------------------------------
  # Organizar resultados
  # ----------------------------------------------------------
  
  tabela <- resultado |>
    select(
      amostra,
      ATT,
      SE,
      N,
      baseline
    ) |>
    mutate(
      ATT = round(ATT, 4),
      SE = round(SE, 4),
      baseline = round(baseline, 4)
    )
  
  
  # ----------------------------------------------------------
  # Texto inicial
  # ----------------------------------------------------------
  
  texto <- paste0(
    "============================================================\n",
    "Outcome: ", outcome, "\n",
    "============================================================\n\n",
    
    "Especificação          ATT          SE          N       Baseline\n",
    "------------------------------------------------------------\n"
  )
  
  
  # ----------------------------------------------------------
  # Adicionar cada resultado
  # ----------------------------------------------------------
  
  for (i in seq_len(nrow(tabela))) {
    
    linha <- sprintf(
      "%-22s %-11.4f %-11.4f %-8d %-11.4f\n",
      tabela$amostra[i],
      tabela$ATT[i],
      tabela$SE[i],
      tabela$N[i],
      tabela$baseline[i]
    )
    
    texto <- paste0(
      texto,
      linha
    )
  }
  
  
  # ----------------------------------------------------------
  # Informações da estimação
  # ----------------------------------------------------------
  
  texto <- paste0(
    texto,
    "\n",
    "============================================================\n",
    "Covariáveis:\n",
    "tx_pop_acesso_agua + taxa_cob_imun + pib_pc\n",
    "Método: Callaway & Sant'Anna (DR)\n",
    "Amostra: Only treated\n",
    "Ponderação: pop_total\n",
    "Grupo de controle: Not-yet-treated\n",
    "============================================================\n"
  )
  
  
  # ----------------------------------------------------------
  # Salvar
  # ----------------------------------------------------------
  
  writeLines(
    texto,
    nome_arquivo,
    useBytes = TRUE
  )
  
  
  return(
    nome_arquivo
  )
}


# ============================================================
# 6. FUNÇÃO DE ESTIMAÇÃO DID
# ============================================================

estimacao_did <- function(
    outcome,
    data = dados_saneamento,
    xformla = xformula,
    tname = "ano",
    idname = "cod_mun",
    gname = "primeiro_tratamento",
    weightsname = "pop_total",
    clustervars = "cod_mun",
    est_method = "dr",
    control_group = "notyettreated",
    base_period = "varying",
    min_e = -7,
    max_e = 7,
    path = "resultados_ponderados_pela_populacao/"
) {
  
  
  y_sym <- rlang::sym(
    outcome
  )
  
  
  # ==========================================================
  # ONLY TREATED
  # ==========================================================
  
  dados_so_tratados <- data |>
    filter(
      primeiro_tratamento != 0
    ) |>
    mutate(
      primeiro_tratamento = if_else(
        primeiro_tratamento == 2018,
        0L,
        primeiro_tratamento
      )
    )
  
  
  # ==========================================================
  # INFORMAÇÕES NO CONSOLE
  # ==========================================================
  
  cat("\n")
  cat("------------------------------------------\n")
  cat("Outcome:", outcome, "\n")
  cat("Amostra: Only treated\n")
  cat(
    "Covariáveis:",
    paste(covariates, collapse = ", "),
    "\n"
  )
  cat("Ponderação: pop_total\n")
  cat("------------------------------------------\n")
  
  
  # ==========================================================
  # CALLAWAY & SANT'ANNA
  # ==========================================================
  
  modelo_cs <- did::att_gt(
    yname = outcome,
    tname = tname,
    idname = idname,
    gname = gname,
    
    # COVARIÁVEIS
    xformla = xformla,
    
    panel = FALSE,
    allow_unbalanced_panel = FALSE,
    control_group = control_group,
    
    # PONDERAÇÃO PELA POPULAÇÃO
    weightsname = weightsname,
    
    clustervars = clustervars,
    est_method = est_method,
    base_period = base_period,
    data = dados_so_tratados
  )
  
  
  # ==========================================================
  # AGREGAÇÃO DINÂMICA
  # ==========================================================
  
  output_cs <- did::aggte(
    modelo_cs,
    type = "dynamic",
    na.rm = TRUE,
    min_e = min_e,
    max_e = max_e
  )
  
  
  crit_val <- output_cs$crit.val.egt
  
  
  # ==========================================================
  # BASELINE
  # ==========================================================
  
  baseline <- dados_so_tratados |>
    mutate(
      tratado = if_else(
        primeiro_tratamento > 0,
        1,
        0
      ),
      
      ano_relativo =
        primeiro_tratamento - ano
    ) |>
    filter(
      ano_relativo == -1,
      tratado == 1
    ) |>
    pull(
      !!y_sym
    ) |>
    mean(
      na.rm = TRUE
    )
  
  
  # ==========================================================
  # EVENT STUDY
  # ==========================================================
  
  df <- data.frame(
    t_label = output_cs$egt,
    coeficientes = output_cs$att.egt,
    se = output_cs$se.egt
  ) |>
    
    tidyr::drop_na(
      se
    ) |>
    
    mutate(
      t_group = case_when(
        t_label <= -5 ~ -5,
        t_label >= 5 ~ 5,
        TRUE ~ t_label
      )
    ) |>
    
    group_by(
      t_group
    ) |>
    
    summarise(
      coeficientes = mean(
        coeficientes,
        na.rm = TRUE
      ),
      
      se = sqrt(
        mean(
          se^2,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
    ) |>
    
    mutate(
      ymin =
        coeficientes -
        crit_val * se,
      
      ymax =
        coeficientes +
        crit_val * se,
      
      t_label = case_when(
        t_group == -5 ~ "-5",
        t_group == 5 ~ "+5",
        TRUE ~ as.character(t_group)
      ),
      
      amostra = "Only treated"
    ) |>
    
    mutate(
      t_label = factor(
        t_label,
        levels = c(
          "-5",
          "-4",
          "-3",
          "-2",
          "-1",
          "0",
          "1",
          "2",
          "3",
          "4",
          "+5"
        )
      )
    )
  
  
  # ==========================================================
  # RESULTADO AGREGADO
  # ==========================================================
  
  resultados <- tibble::tibble(
    amostra = "Only treated",
    outcome = outcome,
    ATT = output_cs$overall.att,
    SE = output_cs$overall.se,
    N = as.integer(
      output_cs$DIDparams$n
    ),
    baseline = baseline
  )
  
  
  # ==========================================================
  # SALVAR TABELA
  # ==========================================================
  
  salvar_tabela_txt(
    resultados,
    outcome,
    path
  )
  
  
  # ==========================================================
  # PLOT
  # ==========================================================
  
  plot <- ggplot(
    df,
    aes(
      x = as.numeric(t_label),
      y = coeficientes
    )
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 1.2
    ) +
    
    geom_point(
      size = 3
    ) +
    
    geom_errorbar(
      aes(
        ymin = ymin,
        ymax = ymax
      ),
      width = 0
    ) +
    
    scale_x_continuous(
      breaks = 1:11,
      labels = levels(
        df$t_label
      )
    ) +
    
    xlab(
      "Relative time to treatment"
    ) +
    
    ylab(
      "Coefficients"
    ) +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      axis.text = element_text(
        size = 22
      ),
      axis.title = element_text(
        size = 25
      )
    )
  
  
  # ==========================================================
  # SALVAR PLOT
  # ==========================================================
  
  ggsave(
    filename = paste0(
      outcome,
      "_pop_weighted.pdf"
    ),
    plot = plot,
    device = "pdf",
    path = path,
    width = 14,
    height = 8.5
  )
  
  
  print(plot)
  
  
  # ==========================================================
  # RETORNAR RESULTADOS
  # ==========================================================
  
  return(
    resultados
  )
}


# ============================================================
# 7. OUTCOMES - RESÍDUOS SÓLIDOS
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
# 8. OUTCOMES - SAÚDE
# ============================================================

outcomes_saude <- c(
  "tx_inter_feco_oral",
  "tx_inter_contato_agua",
  "tx_inter_higiene"
)


resultados2 <- purrr::map_dfr(
  outcomes_saude,
  ~ estimacao_did(.x)
)


# ============================================================
# 9. RESULTADOS FINAIS
# ============================================================

tabela_final <- bind_rows(
  resultados1,
  resultados2
)


# Visualizar
tabela_final
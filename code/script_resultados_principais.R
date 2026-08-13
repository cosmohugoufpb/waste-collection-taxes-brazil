# Paper: Taxa de coleta de resíduos sólidos
# Estimação - Full Sample vs. Only Treated
# ============================================================


# ============================================================
# 1. CARREGAR PACOTES
# ============================================================

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
  filter(
    ano == 2009
  ) |>
  mutate(
    tratado = if_else(
      primeiro_tratamento > 0,
      1,
      0
    )
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


# Adicionar pesos
dados_ps <- dados_ps |>
  mutate(
    ps_weight = ps_weight$weights
  ) |>
  select(
    cod_mun,
    ps_weight
  )


# ============================================================
# 4. ADICIONAR PESOS À BASE PRINCIPAL
# ============================================================

dados_saneamento <- dados_saneamento |>
  left_join(
    dados_ps,
    by = "cod_mun"
  )


# ============================================================
# 5. CRIAR PASTA DE RESULTADOS
# ============================================================

if (!dir.exists("resultados_principais")) {
  dir.create(
    "resultados_principais",
    recursive = TRUE
  )
}


# ============================================================
# 6. COVARIÁVEIS
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
# 7. FUNÇÃO PARA SALVAR TABELA EM TXT
# ============================================================

salvar_tabela_txt <- function(
    resultado,
    outcome,
    path = "resultados_principais/"
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
  # Adicionar cada linha
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
# 8. FUNÇÃO DE ESTIMAÇÃO DID
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
    min_e = -7,
    max_e = 7,
    path = "resultados_principais/"
) {
  
  
  y_sym <- rlang::sym(
    outcome
  )
  
  
  # ==========================================================
  # FUNÇÃO INTERNA PARA RODAR O DID
  # ==========================================================
  
  rodar_did <- function(
    base,
    nome_amostra
  ) {
    
    
    cat("\n")
    cat("------------------------------------------\n")
    cat("Outcome:", outcome, "\n")
    cat("Amostra:", nome_amostra, "\n")
    cat("------------------------------------------\n")
    
    
    # --------------------------------------------------------
    # Callaway & Sant'Anna
    # --------------------------------------------------------
    
    modelo_cs <- did::att_gt(
      yname = outcome,
      tname = tname,
      idname = idname,
      gname = gname,
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
    
    
    # --------------------------------------------------------
    # Agregação dinâmica
    # --------------------------------------------------------
    
    output_cs <- did::aggte(
      modelo_cs,
      type = "dynamic",
      na.rm = TRUE,
      min_e = min_e,
      max_e = max_e
    )
    
    
    crit_val <- output_cs$crit.val.egt
    
    
    # --------------------------------------------------------
    # Baseline
    # --------------------------------------------------------
    
    baseline <- base |>
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
    
    
    # --------------------------------------------------------
    # Event study
    # --------------------------------------------------------
    
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
        
        amostra = nome_amostra
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
    
    
    # --------------------------------------------------------
    # Resultado agregado
    # --------------------------------------------------------
    
    resultados <- tibble::tibble(
      amostra = nome_amostra,
      outcome = outcome,
      ATT = output_cs$overall.att,
      SE = output_cs$overall.se,
      N = as.integer(
        output_cs$DIDparams$n
      ),
      baseline = baseline
    )
    
    
    return(
      list(
        df = df,
        res = resultados
      )
    )
  }
  
  
  # ==========================================================
  # FULL SAMPLE
  # ==========================================================
  
  full <- rodar_did(
    data,
    "Full sample"
  )
  
  
  # ==========================================================
  # ONLY TREATED
  # ==========================================================
  #
  # Mesmo filtro utilizado anteriormente:
  #   primeiro_tratamento != 0
  #
  # E 2018 é transformado em 0.
  #
  
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
  
  
  treated <- rodar_did(
    dados_so_tratados,
    "Only treated"
  )
  
  
  # ==========================================================
  # SALVAR TABELA TXT
  # ==========================================================
  
  resultados_outcome <- bind_rows(
    full$res,
    treated$res
  )
  
  
  salvar_tabela_txt(
    resultados_outcome,
    outcome,
    path
  )
  
  
  # ==========================================================
  # JUNTAR EVENT STUDIES
  # ==========================================================
  
  df_plot <- bind_rows(
    full$df,
    treated$df
  ) |>
    
    mutate(
      x_num = as.numeric(
        t_label
      ),
      
      x_plot = ifelse(
        amostra == "Full sample",
        x_num - 0.15,
        x_num + 0.15
      )
    )
  
  
  # ==========================================================
  # PLOT
  # ==========================================================
  
  plot <- ggplot(
    df_plot,
    aes(
      x = x_plot,
      y = coeficientes,
      color = amostra,
      shape = amostra
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
        df_plot$t_label
      )
    ) +
    
    scale_color_manual(
      values = c(
        "Full sample" = "#006D77",
        "Only treated" = "#E29578"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Full sample" = 16,
        "Only treated" = 15
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
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text = element_text(
        size = 22
      ),
      axis.title = element_text(
        size = 25
      ),
      legend.text = element_text(
        size = 27
      )
    )
  
  
  # ==========================================================
  # SALVAR PLOT
  # ==========================================================
  
  ggsave(
    filename = paste0(
      outcome,
      "_full_vs_treated.pdf"
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
    resultados_outcome
  )
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
  "tx_inter_higiene"
)


resultados2 <- purrr::map_dfr(
  outcomes_saude,
  ~ estimacao_did(.x)
)


# ============================================================
# 11. RESULTADOS FINAIS
# ============================================================

tabela_final <- bind_rows(
  resultados1,
  resultados2
)


# Visualizar
tabela_final
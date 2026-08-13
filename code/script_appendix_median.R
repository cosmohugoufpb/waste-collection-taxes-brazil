# Paper: Taxa de coleta de resíduos sólidos
# Estimação - Heterogeneidade por população média pré-tratamento
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
  tratado ~ tx_pop_acesso_agua +
    taxa_cob_imun +
    pib_pc +
    tx_inter_feco_oral +
    tx_inter_inseto_vetor,
  data = dados_ps,
  method = "ebal",
  estimand = "ATT"
)

summary(ps_weight)


dados_ps <- dados_ps |>
  mutate(
    ps_weight = ps_weight$weights
  ) |>
  select(
    cod_mun,
    ps_weight
  )


# Adicionar pesos à base principal
dados_saneamento <- dados_saneamento |>
  left_join(
    dados_ps,
    by = "cod_mun"
  )


# ============================================================
# 4. POPULAÇÃO MÉDIA PRÉ-TRATAMENTO
# ============================================================
#
# Utiliza a variável POP_TOTAL_PRE.
#
# Para cada município:
#   1. calcula a média de pop_total_pre;
#   2. calcula a mediana das médias municipais;
#   3. classifica os municípios acima/abaixo da mediana.
#
# ============================================================

media_pop_municipio <- dados_saneamento |>
  group_by(cod_mun) |>
  summarise(
    pop_media = mean(
      pop_total_pre,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    pop_media = if_else(
      is.nan(pop_media),
      NA_real_,
      pop_media
    )
  )


# ============================================================
# 5. MEDIANA DAS MÉDIAS MUNICIPAIS
# ============================================================

mediana_pop <- median(
  media_pop_municipio$pop_media,
  na.rm = TRUE
)

print(mediana_pop)


# ============================================================
# 6. CLASSIFICAÇÃO ACIMA / ABAIXO DA MEDIANA
# ============================================================

media_pop_municipio <- media_pop_municipio |>
  mutate(
    pop_above = if_else(
      pop_media > mediana_pop,
      1,
      0,
      missing = NA_real_
    ),
    
    pop_below = if_else(
      pop_media <= mediana_pop,
      1,
      0,
      missing = NA_real_
    )
  )


# Conferência
table(
  media_pop_municipio$pop_above,
  useNA = "ifany"
)

table(
  media_pop_municipio$pop_below,
  useNA = "ifany"
)


# ============================================================
# 7. ADICIONAR CLASSIFICAÇÃO À BASE PRINCIPAL
# ============================================================

dados_saneamento <- dados_saneamento |>
  left_join(
    media_pop_municipio,
    by = "cod_mun"
  )


# ============================================================
# 8. CRIAR PASTA APPENDIX
# ============================================================

if (!dir.exists("appendix")) {
  dir.create(
    "appendix",
    recursive = TRUE
  )
}


# ============================================================
# 9. COVARIÁVEIS
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

xformula <- as.formula(xformula_str)


# ============================================================
# 10. FUNÇÃO PARA SALVAR RESULTADOS EM TXT
# ============================================================

salvar_tabela_txt <- function(
    resultado,
    outcome,
    grupo,
    path = "appendix/"
) {
  
  nome_arquivo <- paste0(
    path,
    outcome,
    "_",
    grupo,
    ".txt"
  )
  
  
  # Criar tabela simples
  tabela <- resultado |>
    select(
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
  
  
  # Criar texto
  texto <- paste0(
    "====================================================\n",
    "Outcome: ", outcome, "\n",
    "Grupo: ", grupo, "\n",
    "====================================================\n\n",
    
    "ATT        SE        N        Baseline\n",
    "----------------------------------------------------\n",
    
    sprintf(
      "%-10.4f %-10.4f %-8d %-10.4f\n",
      tabela$ATT,
      tabela$SE,
      tabela$N,
      tabela$baseline
    )
  )
  
  
  # Salvar como TXT
  writeLines(
    texto,
    nome_arquivo,
    useBytes = TRUE
  )
  
  
  return(nome_arquivo)
}


# ============================================================
# 11. FUNÇÃO DE ESTIMAÇÃO DID
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
    path = "appendix/"
) {
  
  
  y_sym <- rlang::sym(outcome)
  
  
  # ==========================================================
  # 11.1 FILTRO GLOBAL
  # ==========================================================
  #
  # Somente municípios tratados.
  #
  # 2018 é transformado em 0 e utilizado como controle.
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
  
  
  # ==========================================================
  # 11.2 FUNÇÃO INTERNA PARA RODAR O DID
  # ==========================================================
  
  rodar_did <- function(
    base,
    nome_amostra
  ) {
    
    
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
    
    
    output_cs <- did::aggte(
      modelo_cs,
      type = "dynamic",
      na.rm = TRUE,
      min_e = min_e,
      max_e = max_e
    )
    
    
    crit_val <- output_cs$crit.val.egt
    
    
    # ========================================================
    # BASELINE
    # ========================================================
    
    baseline <- base |>
      mutate(
        tratado = if_else(
          primeiro_tratamento > 0,
          1,
          0
        ),
        
        ano_relativo = primeiro_tratamento - ano
      ) |>
      filter(
        ano_relativo == -1,
        tratado == 1
      ) |>
      pull(!!y_sym) |>
      mean(
        na.rm = TRUE
      )
    
    
    # ========================================================
    # EVENT STUDY
    # ========================================================
    
    df <- data.frame(
      t_label = output_cs$egt,
      coeficientes = output_cs$att.egt,
      se = output_cs$se.egt
    ) |>
      tidyr::drop_na(se) |>
      mutate(
        t_group = case_when(
          t_label <= -5 ~ -5,
          t_label >= 5 ~ 5,
          TRUE ~ t_label
        )
      ) |>
      group_by(t_group) |>
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
        ymin = coeficientes - crit_val * se,
        ymax = coeficientes + crit_val * se,
        
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
    
    
    # ========================================================
    # RESULTADOS AGREGADOS
    # ========================================================
    
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
  # 11.3 ABOVE MEDIAN
  # ==========================================================
  
  dados_above <- dados_so_tratados |>
    filter(
      pop_above == 1
    )
  
  
  above <- rodar_did(
    dados_above,
    "Above median"
  )
  
  
  # ==========================================================
  # 11.4 BELOW MEDIAN
  # ==========================================================
  
  dados_below <- dados_so_tratados |>
    filter(
      pop_below == 1
    )
  
  
  below <- rodar_did(
    dados_below,
    "Below median"
  )
  
  
  # ==========================================================
  # 11.5 SALVAR TXT INDIVIDUAIS
  # ==========================================================
  
  salvar_tabela_txt(
    above$res,
    outcome,
    "pop_above",
    path
  )
  
  
  salvar_tabela_txt(
    below$res,
    outcome,
    "pop_below",
    path
  )
  
  
  # ==========================================================
  # 11.6 JUNTAR PARA PLOT
  # ==========================================================
  
  df_plot <- bind_rows(
    above$df,
    below$df
  ) |>
    mutate(
      x_num = as.numeric(t_label),
      
      x_plot = ifelse(
        amostra == "Above median",
        x_num - 0.15,
        x_num + 0.15
      )
    )
  
  
  # ==========================================================
  # 11.7 PLOT FINAL
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
        "Above median" = "#006D77",
        "Below median" = "#E29578"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Above median" = 16,
        "Below median" = 15
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
  # 11.8 SALVAR PLOT
  # ==========================================================
  
  ggsave(
    filename = paste0(
      outcome,
      "_pop_comparison.pdf"
    ),
    plot = plot,
    device = "pdf",
    path = path,
    width = 14,
    height = 8.5
  )
  
  
  print(plot)
  
  
  # ==========================================================
  # 11.9 RETORNAR RESULTADOS
  # ==========================================================
  
  return(
    bind_rows(
      above$res |>
        mutate(
          grupo_pop = "pop_above"
        ),
      
      below$res |>
        mutate(
          grupo_pop = "pop_below"
        )
    )
  )
}


# ============================================================
# 12. OUTCOMES - RESÍDUOS SÓLIDOS
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
# 13. OUTCOMES - SAÚDE
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
# 14. RESULTADOS FINAIS
# ============================================================

tabela_final <- bind_rows(
  resultados1,
  resultados2
) |>
  pivot_wider(
    names_from = grupo_pop,
    values_from = c(
      ATT,
      SE,
      N,
      baseline
    ),
    names_sep = "_"
  )


# Visualizar
tabela_final


# ============================================================
# 15. SALVAR CLASSIFICAÇÃO DOS MUNICÍPIOS
# ============================================================

write.csv(
  media_pop_municipio,
  file = "appendix/classificacao_populacao_municipios.csv",
  row.names = FALSE
)


# ============================================================
# 16. SALVAR MEDIANA
# ============================================================

writeLines(
  paste0(
    "Mediana da população média pré-tratamento = ",
    mediana_pop
  ),
  "appendix/mediana_populacao.txt"
)
```


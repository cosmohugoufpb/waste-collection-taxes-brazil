# Paper: Taxa de coleta de residuos solidos
# Estimacao


# Carregar pacotes ----
pacman::p_load(
  tidyverse, did, showtext, gt, WeightIt, kableExtra
)

# Renderizacao de fontes
font_family <- "STIX Two Text"
sysfonts::font_add_google(font_family)
showtext::showtext_auto()


# Ler dados
dados_saneamento <- readr::read_csv(
  "data/dados_saneamento.csv.gz"
)

# Calcular o propensity score -----
dados_ps <- dados_saneamento |>
  filter(
    ano == 2009
  ) |> 
  mutate(
    tratado = if_else(primeiro_tratamento > 0, 1, 0)
  )

ps_weight <- weightit(
  tratado ~ tx_pop_acesso_agua + taxa_cob_imun + pib_pc  + tx_inter_feco_oral + tx_inter_inseto_vetor, 
  #tx_inter_contato_agua + tx_inter_higiente + tx_inter_teniase,
  data = dados_ps,
  method = "ebal",
  estimand = "ATT"
) 
summary(ps_weight)


dados_ps <- dados_ps |> 
  mutate(ps_weight = ps_weight$weights) |> 
  select(cod_mun, ps_weight)

# Salvar dados
dados_saneamento <- dados_saneamento |> 
  left_join(
    dados_ps
  )

# Estimacao: Tratados e nao tratados -----

# Covariaveis
covariates <- c("tx_pop_acesso_agua", "taxa_cob_imun", "pib_pc")
xformula_str <- paste("~", paste(covariates, collapse = " + "))
xformula <- as.formula(xformula_str)


# Funcao para estimacao

estimacao_did <- function(outcome,
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
                          min_e = -7, max_e = 7,
                          path = "figuras/") {
  
  y_sym <- rlang::sym(outcome)
  
  # =====================================================
  # FUNÇÃO INTERNA (RODA DID)
  # =====================================================
  
  rodar_did <- function(base, nome_amostra){
    
    modelo_cs <- did::att_gt(
      yname = outcome,
      tname = tname,
      idname = idname,
      gname = gname,
      xformla = xformula, 
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
    
    # =========================
    # BASELINE
    # =========================
    
    baseline <- base |> 
      mutate(
        tratado = if_else(primeiro_tratamento > 0, 1, 0),
        ano_relativo = primeiro_tratamento - ano
      ) |> 
      filter(ano_relativo == -1, tratado == 1) |> 
      pull(!!y_sym) %>%
      mean(na.rm = TRUE)
    
    # =========================
    # EVENT STUDY (PLOT)
    # =========================
    
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
        coeficientes = mean(coeficientes, na.rm = TRUE),
        se = sqrt(mean(se^2, na.rm = TRUE)),
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
          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4", "+5")
        )
      )
    
    # =========================
    # RESULTADOS AGREGADOS
    # =========================
    
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
  
  # =====================================================
  # FULL SAMPLE
  # =====================================================
  
  full <- rodar_did(data, "Full sample")
  
  # =====================================================
  # ONLY TREATED
  # =====================================================
  
  dados_treated <- data %>%
    filter(primeiro_tratamento != 0) %>%
    mutate(primeiro_tratamento = ifelse(primeiro_tratamento == 2018, 0, primeiro_tratamento))
  
  treated <- rodar_did(dados_treated, "Only treated")
  
  # =====================================================
  # JUNTAR PARA PLOT
  # =====================================================
  
  df_plot <- bind_rows(full$df, treated$df) |>
    
    # 🔥 deslocamento lateral (evita sobreposição)
    mutate(
      x_num = as.numeric(t_label),
      x_plot = ifelse(amostra == "Full sample",
                      x_num - 0.15,
                      x_num + 0.15)
    )
  
  # =====================================================
  # PLOT FINAL
  # =====================================================
  
  plot <- ggplot(df_plot, aes(x = x_plot, y = coeficientes,
                              color = amostra, shape = amostra)) +
    
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
    
    geom_point(size = 3) +
    
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
    
    scale_x_continuous(
      breaks = 1:11,
      labels = levels(df_plot$t_label)
    ) +
    
    scale_color_manual(values = c(
      "Full sample" = "#006D77",
      "Only treated" = "#E29578"
    )) +
    
    scale_shape_manual(values = c(
      "Full sample" = 16,
      "Only treated" = 15
    )) +
    
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
  
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  
  ggsave(
    filename = paste0(outcome, "_duplo.pdf"),
    plot = plot,
    device = "pdf",
    path = path,
    width = 14, height = 8.5
  )
  
  print(plot)
  
  # =====================================================
  # RETORNAR RESULTADOS PARA TABELA
  # =====================================================
  
  return(bind_rows(full$res, treated$res))
}


# =====================================================
# RESULTADOS (TABELAS)
# =====================================================

outcomes <- c(
  "d_plano_gestao_residuos", 
  "existe_lixao", 
  "d_coleta_seletiva", 
  "tx_pop_resid_solidos",
  "tx_pop_coleta_diaria",
  "tx_pop_coleta_2_3_semana"
)

resultados1 <- purrr::map_dfr(outcomes, ~ estimacao_did(.x))


outcomes_saude <- c(
  "tx_inter_feco_oral", 
  "tx_inter_inseto_vetor", 
  "tx_inter_contato_agua", 
  "tx_inter_higiene"
)

resultados2 <- purrr::map_dfr(outcomes_saude, ~ estimacao_did(.x))


# =====================================================
# TABELA FINAL
# =====================================================

library(tidyr)

tabela_final <- bind_rows(resultados1, resultados2) |>
  pivot_wider(
    names_from = amostra,
    values_from = c(ATT, SE, N, baseline)
  )

tabela_final

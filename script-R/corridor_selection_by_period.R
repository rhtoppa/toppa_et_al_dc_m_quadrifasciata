# ===============================================================
# TITLE: Seasonal Corridor Prioritization for Melipona quadrifasciata Movement 
#        Using GLMM-derived Metrics
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script implements a corridor prioritization system based on seasonal periods
# using significant metrics identified through GLMM analysis. The composite score
# combines MAX_NDVI_CORR and RANGE_NDWI_CORR with predefined weights to identify
# top-N corridors per season and provides interseasonal summaries for 
# Melipona quadrifasciata conservation planning.
# ===============================================================

# Packages -----------------------------------------------------
req <- c("readxl","dplyr","tidyr","stringr","readr","tibble","scales","purrr")
new <- req[!(req %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(dplyr); library(tidyr); library(stringr)
library(readr);  library(tibble); library(scales); library(purrr)

# Adjust the parameters for data input and output----------------------
XLSX    <- "C:/r_studio/data_statistic.xlsx"
OUTDIR  <- "C:/r_studio/corridors_by_period"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Score weights (adjusted according to GLMM: fixed or standardized z-beta values)
W_NDVI  <- 0.60
W_NDWI  <- 0.40

# Maximum practical distance for eligibility (K=5.5 km)
DIST_MAX <- 5.5

# Quantile for "high" NDVI and NDWI filters (can be disabled below)
QTILE_FILTER <- 0.75
APLICAR_FILTRO_QTILE <- TRUE   #Set to FALSE to not filter by quantiles.

# Number of priority lanes per period
TOP_N <- 10

# Read data and check columns. --------------------------------
raw <- read_excel(XLSX)

need <- c("period","release_point","distance_km","MAX_NDVI_CORR","RANGE_NDWI_CORR")
miss <- setdiff(need, names(raw))
if (length(miss)) stop("Faltam colunas na planilha: ", paste(miss, collapse=", "))

dat <- raw %>%
  mutate(
    period         = factor(period),
    release_point  = as.factor(release_point),
    distance_km    = suppressWarnings(as.numeric(distance_km)),
    MAX_NDVI_CORR  = suppressWarnings(as.numeric(MAX_NDVI_CORR)),
    RANGE_NDWI_CORR= suppressWarnings(as.numeric(RANGE_NDWI_CORR))
  ) %>%
  drop_na(period, release_point, distance_km, MAX_NDVI_CORR, RANGE_NDWI_CORR)

# Function to prioritize by period. --------------------------
priorizar_por_periodo <- function(df_period,
                                  w_ndvi = W_NDVI,
                                  w_ndwi = W_NDWI,
                                  dist_max = DIST_MAX,
                                  qtile = QTILE_FILTER,
                                  aplicar_filtro = APLICAR_FILTRO_QTILE,
                                  top_n = TOP_N) {
  
  # Calculates thresholds by quantile (within the same period).
  thr_ndvi <- stats::quantile(df_period$MAX_NDVI_CORR, probs = qtile, na.rm = TRUE)
  thr_ndwi <- stats::quantile(df_period$RANGE_NDWI_CORR, probs = qtile, na.rm = TRUE)
  
  # Rescale 0–1 for score comparability.
  df_sc <- df_period %>%
    mutate(
      z_ndvi = scales::rescale(MAX_NDVI_CORR,  to = c(0,1)),
      z_ndwi = scales::rescale(RANGE_NDWI_CORR, to = c(0,1)),
      score  = w_ndvi * z_ndvi + w_ndwi * z_ndwi
    )
  
  # Practical filters
  df_filt <- df_sc %>%
    filter(!is.na(score),
           !is.na(distance_km),
           distance_km <= dist_max)
  
  if (aplicar_filtro) {
    df_filt <- df_filt %>%
      filter(MAX_NDVI_CORR  >= thr_ndvi,
             RANGE_NDWI_CORR >= thr_ndwi)
  }
  
  # Aggregation by release_point within the period
  out <- df_filt %>%
    group_by(release_point) %>%
    summarise(
      n_obs       = n(),
      dist_km_med = median(distance_km, na.rm = TRUE),
      max_NDVI    = max(MAX_NDVI_CORR, na.rm = TRUE),
      rng_NDWI    = median(RANGE_NDWI_CORR, na.rm = TRUE),
      score_mean  = mean(score, na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    arrange(desc(score_mean)) %>%
    mutate(rank_period = row_number()) %>%
    slice_head(n = top_n)   # ✅ versão corrigida
}
  
  # Execute by period ---------------------------------------
  periodos <- levels(dat$period)
  if (is.null(periodos)) periodos <- sort(unique(dat$period))
  
  lista_por_periodo <- lapply(periodos, function(per) {
    sub <- dat %>% filter(period == per)
    top_tbl <- priorizar_por_periodo(sub)
    if (nrow(top_tbl) > 0) {
      top_tbl %>%
        mutate(period = per) %>%
        relocate(period, .before = release_point)
    } else {
      tibble(period = per,
             release_point = factor(NA),
             n_obs = NA_integer_,
             dist_km_med = NA_real_,
             max_NDVI = NA_real_,
             rng_NDWI = NA_real_,
             score_mean = NA_real_,
             rank_period = NA_integer_)
    }
  })
  
  prioritarios_por_periodo <- bind_rows(lista_por_periodo)
  
  # Save and display in console --------------------------------
  out_per_file <- file.path(OUTDIR, "candidate_corridors_TOPN_by_period.csv")
  write_csv(prioritarios_por_periodo, out_per_file)
  
  cat("\n=== TOP-N por período (salvo em):\n", out_per_file, "\n")
  print(prioritarios_por_periodo, n = Inf)
  
  # Interseasonal summary (frequency, average, standard deviation, best placement) -------
  # Consider only valid lines (with non-NA release_point)
  priori_valid <- prioritarios_por_periodo %>%
    filter(!is.na(release_point))
  
  summary_all <- priori_valid %>%
    group_by(release_point) %>%
    summarise(
      periods_in_topN    = n(),                              # In how many periods did it enter the TOP-N
      mean_score_across  = mean(score_mean, na.rm = TRUE),   # average score
      sd_score_across    = sd(score_mean, na.rm = TRUE),     # variability
      best_rank_across   = min(rank_period, na.rm = TRUE),   # best placement (lowest rank)
      median_dist_km     = median(dist_km_med, na.rm = TRUE),
      max_NDVI_overall   = max(max_NDVI, na.rm = TRUE),
      median_rng_NDWI    = median(rng_NDWI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(periods_in_topN), best_rank_across, desc(mean_score_across)) %>%
    mutate(global_rank = row_number())
  
  out_summary_file <- file.path(OUTDIR, "candidate_corridors_INTERSEASON_summary.csv")
  write_csv(summary_all, out_summary_file)
  
  cat("\n=== Sumário intersazonal (salvo em):\n", out_summary_file, "\n")
  print(summary_all, n = Inf)
  
  # List of "core" and "seasonal" periods ----------------------
    # Core = appears in >= 3 periods; Seasonal = appears in 1-2 periods
  core_thresh <- 3
  core_corridors <- summary_all %>% filter(periods_in_topN >= core_thresh)
  seasonal_corridors <- summary_all %>% filter(periods_in_topN < core_thresh)
  
  write_csv(core_corridors,    file.path(OUTDIR, "candidate_corridors_CORE.csv"))
  write_csv(seasonal_corridors,file.path(OUTDIR, "candidate_corridors_SEASONAL.csv"))
  
  cat("\nArquivos gerados em:\n", OUTDIR, "\n")

  

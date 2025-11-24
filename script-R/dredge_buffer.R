# =============================================================================
# TITLE: GLMM Model Selection for Buffer Landscape Metrics Using Dredge Algorithm
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script performs automated model selection using the dredge algorithm
# on GLMM binomial models with buffer landscape metrics. Distance is forced
# in all models, and random effects account for colony, batch, and release point
# structure in Melipona quadrifasciata return probability analysis.
# ==========================================

# Packages ---------------------------------------------------------------
req <- c(
  "readxl","dplyr","tidyr","glmmTMB","MuMIn","DHARMa",
  "broom.mixed","performance","ggplot2","stringr","ggeffects",
  "scales","readr","tibble"
)
new <- req[!(req %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(dplyr); library(tidyr); library(glmmTMB); library(MuMIn)
library(DHARMa); library(broom.mixed); library(performance); library(ggplot2)
library(stringr); library(ggeffects); library(scales); library(readr); library(tibble)

# Adjust the parameters for data input and output---------------------------------------------
xlsx   <- "C:/r_studio/statistic_data.xlsx"
OUTDIR <- "C:/r_studio/dredge_buffer"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# 2) Read data and minimum typing -------------------------------------------
dat <- read_excel(xlsx)

df <- dat %>%
  mutate(
    event_use     = as.integer(event_use),
    distance_km   = suppressWarnings(as.numeric(distance_km)),
    colony        = factor(colony),
    release_point = factor(release_point),
    release_batch = as.Date(release_batch)
  )

# Define the candidates using the BUFFER (adjust freely) -----------------------
CAND_BUFFER <- c(
  "MEAN_NDVI_BUFFER","MAX_NDVI_BUFFER","RANGE_NDVI_BUFFER",
  "MEAN_NDWI_BUFFER","MAX_NDWI_BUFFER","RANGE_NDWI_BUFFER",
  "MEAN_MR_BUFFER","MAX_MR_BUFFER","RANGE_MR_BUFFER"
)

# It only keeps the ones that exist and are numerical.
CAND_BUFFER <- CAND_BUFFER[CAND_BUFFER %in% names(df)]
CAND_BUFFER <- CAND_BUFFER[sapply(df[CAND_BUFFER], is.numeric)]
if (length(CAND_BUFFER) == 0) stop("Nenhuma variável *_BUFFER numérica encontrada nas colunas.")

# Remove candidates with variance ~ 0 (avoid rank-deficient)
nzv <- sapply(df[CAND_BUFFER], function(x) stats::sd(x, na.rm = TRUE))
CAND_BUFFER <- CAND_BUFFER[nzv > 0]
if (length(CAND_BUFFER) == 0) stop("Todos os *_BUFFER têm variância ~ 0 após checagem.")

# Formula FULL (distance + all candidates) ---------------------------
rhs_fixed <- paste(c("distance_km", CAND_BUFFER), collapse = " + ")
rhs_rand  <- "(1|colony) + (1|release_batch) + (1|release_point)"
FORM_FULL <- as.formula(paste0("event_use ~ ", rhs_fixed, " + ", rhs_rand))

cat("\nFórmula FULL (BUFFER) para dredge:\n"); print(FORM_FULL)

# Subset without NA in the used columns
vars_used <- unique(c("event_use","distance_km", CAND_BUFFER, "colony","release_batch","release_point"))
dsub <- df %>% dplyr::select(dplyr::all_of(vars_used)) %>% tidyr::drop_na()
cat("N linhas após remoção de NAs:", nrow(dsub), "\n")

# Adjust FULL and run dredge (forced distance) -------------------------
m_full <- glmmTMB(FORM_FULL, family = binomial(link="logit"), data = dsub)

# see collinearity in FULL
suppressWarnings({
  colin_full <- performance::check_collinearity(m_full)
  capture.output(colin_full, file = file.path(OUTDIR,"collinearity_FULL_BUFFER.txt"))
})

old_na <- getOption("na.action")
options(na.action = "na.fail")  # exigido pelo dredge
dd <- MuMIn::dredge(m_full, fixed = "distance_km", rank = "AICc")
options(na.action = old_na)

# saves complete model table
write.csv(as.data.frame(dd), file.path(OUTDIR,"AICc_table_BUFFER_dredge.csv"), row.names = FALSE)

cat("\nTop-5 modelos (BUFFER) por AICc:\n")
print(as_tibble(head(dd, 5)))

# Select the top-1 model and also mark the models with ΔAICc ≤ 2
best      <- get.models(dd, 1)[[1]]
delta2_ix <- which(dd$delta <= 2)
write.csv(as.data.frame(dd[delta2_ix,]),
          file.path(OUTDIR,"AICc_models_delta2_BUFFER.csv"), row.names = FALSE)

# Main outputs of the best model ------------------------------------
sink(file.path(OUTDIR,"best_model_summary_BUFFER.txt"))
cat("BEST MODEL (BUFFER)\n\n")
print(summary(best))
cat(sprintf("\nAICc = %.2f\n", MuMIn::AICc(best)))
print(performance::r2(best))   # R2 marginal/condicional (Nakagawa)
sink()

cat("\nResumo do melhor modelo gravado em: best_model_summary_BUFFER.txt\n")
cat(sprintf("AICc (best) = %.2f\n", MuMIn::AICc(best)))
print(performance::r2(best))

# Odds Ratios (distance and indexes present) ---------------------------
fx <- broom.mixed::tidy(best, effects="fixed", conf.int=TRUE)

# distance
if ("distance_km" %in% fx$term) {
  rr <- fx[fx$term=="distance_km", ]
  cat(sprintf("\nDistance OR per 1 km = %.3f (%.3f–%.3f)\n",
              exp(rr$estimate), exp(rr$conf.low), exp(rr$conf.high)))
  cat(sprintf("Distance OR per 0.5 km = %.3f (%.3f–%.3f)\n",
              exp(0.5*rr$estimate), exp(0.5*rr$conf.low), exp(0.5*rr$conf.high)))
}

# LULC indices (only those that entered the best model)
idx_terms <- fx$term[grepl("NDVI|NDWI|MR", fx$term, ignore.case = TRUE)]
if (length(idx_terms)) {
  cat("\nIndex effects (OR per +0.1):\n")
  or_table <- lapply(idx_terms, function(t) {
    r <- fx[fx$term==t,]
    data.frame(
      term = t,
      OR   = exp(r$estimate*0.1),
      LCL  = exp(r$conf.low*0.1),
      UCL  = exp(r$conf.high*0.1)
    )
  }) %>% dplyr::bind_rows()
  print(as_tibble(or_table))
  write.csv(or_table, file.path(OUTDIR,"OR_index_plus0.1_best_BUFFER.csv"), row.names = FALSE)
}

# Forest plot of the fixed elements -------------------------------------------------
fx_plot <- fx %>%
  mutate(
    OR  = exp(estimate),
    LCL = exp(conf.low),
    UCL = exp(conf.high),
    label = dplyr::case_when(
      term == "distance_km" ~ "Distance (per 1 km)",
      TRUE ~ stringr::str_replace_all(term, "_", " ")
    )
  ) %>% filter(term != "(Intercept)")

p_forest <- ggplot(fx_plot, aes(x = reorder(label, OR), y = OR)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = .15) +
  coord_flip() +
  labs(x = NULL, y = "Odds ratio (95% CI)",
       title = "Fixed effects on return probability (BUFFER)") +
  theme_minimal(base_size = 12)

ggsave(file.path(OUTDIR,"forest_OR_best_BUFFER.png"), p_forest, width=7, height=5, dpi=300)
print(p_forest)

# Predictions by distance (0–7.5) under P25/P75 scenarios ------------------
# creates P25/P75 scenarios for EACH index present in the best model (one at a time)
pred_list <- list()
vars_in_best <- intersect(idx_terms, names(dsub))

if (length(vars_in_best)) {
  for (v in vars_in_best) {
    p25 <- stats::quantile(dsub[[v]], 0.25, na.rm=TRUE)
    p75 <- stats::quantile(dsub[[v]], 0.75, na.rm=TRUE)
    
    g_lo <- ggpredict(best, terms = c(
      "distance_km [0:7.5 by=0.25]",
      paste0(v," [", round(p25,3), "]")
    )) %>% as.data.frame() %>% mutate(scenario = paste(v, "low (P25)"))
    
    g_hi <- ggpredict(best, terms = c(
      "distance_km [0:7.5 by=0.25]",
      paste0(v," [", round(p75,3), "]")
    )) %>% as.data.frame() %>% mutate(scenario = paste(v, "high (P75)"))
    
    pred_list[[v]] <- bind_rows(g_lo, g_hi)
  }
  
  pred_df <- bind_rows(pred_list, .id = "index")
  keyd <- c(0,2,4,6,7.5)
  tab_preds <- pred_df %>%
    filter(x %in% keyd) %>%
    transmute(index,
              scenario,
              distance_km = x,
              prob = predicted,
              lo = conf.low,
              hi = conf.high) %>%
    arrange(index, scenario, distance_km)
  
  write_csv(tab_preds, file.path(OUTDIR, "predictions_best_BUFFER_P25P75.csv"))
  print(as_tibble(tab_preds), n = Inf)
}

# DHARMa Diagnosis ---------------------------------------------------
png(file.path(OUTDIR,"DHARMa_best_BUFFER.png"), width=900, height=700, res=120)
plot(DHARMa::simulateResiduals(best, n = 1000))
dev.off()

# Collinearity -----------------------------
suppressWarnings({
  colin <- performance::check_collinearity(best)
  capture.output(colin, file = file.path(OUTDIR,"collinearity_best_BUFFER.txt"))
})

cat("\nArquivos (BUFFER) gerados em:\n", OUTDIR, "\n")


# ==========================================
# TITLE: GLMM Analysis of Corridor Landscape Metrics on 
#        Melipona quadrifasciata Return Probability
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script fits GLMM binomial models to assess the effects of corridor
# landscape metrics on Melipona quadrifasciata return probability. Includes 
# comprehensive model outputs, predictions, and diagnostic checks for
# conservation corridor planning.
# ==========================================

# Packages ---------------------------------------------------------------
req <- c(
  "readxl","dplyr","tidyr","glmmTMB","MuMIn","DHARMa","broom.mixed","performance",
  "ggplot2","stringr","ggeffects","purrr","scales","readr","tibble"
)
new <- req[!(req %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(dplyr); library(tidyr); library(glmmTMB); library(MuMIn)
library(DHARMa); library(broom.mixed); library(performance); library(ggplot2)
library(stringr); library(ggeffects); library(purrr); library(scales); library(readr); library(tibble)

# Adjust the parameters for data input and output -------------------------------------------------------
xlsx   <- "C:/r_studio/statistic_data.xlsx"
OUTDIR <- "C:/r_studio/corridor_glmm"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Read data and minimum typing -------------------------------------------
dat <- read_excel(xlsx)

df <- dat %>%
  mutate(
    event_use     = as.integer(event_use),
    distance_km   = suppressWarnings(as.numeric(distance_km)),
    period        = factor(period),
    colony        = factor(colony),
    release_point = factor(release_point),
    release_batch = as.Date(release_batch)
  )

# WRITE YOUR FORMULA (use the names from the spreadsheet) ----------------------------------
FORM <- event_use ~ distance_km + MAX_NDVI_CORR + RANGE_NDWI_CORR +
  (1|colony) + (1|release_batch) + (1|release_point)

cat("\nFórmula usada:\n"); print(FORM)

# Subset without NA in the used columns --------------------------------
vars_used <- unique(c("event_use", all.vars(FORM)))
missing_vars <- setdiff(vars_used, names(df))
if (length(missing_vars)) stop(paste("Variáveis ausentes na planilha:", paste(missing_vars, collapse=", ")))

dsub <- df %>% select(all_of(vars_used)) %>% tidyr::drop_na()
cat("N linhas após remoção de NAs:", nrow(dsub), "\n")

# GLMM Adjustment ------------------------------------------------------
m <- glmmTMB(FORM, family = binomial(link="logit"), data = dsub)

# It saves to a file and displays in the console
sinkfile <- file.path(OUTDIR, "model_summary.txt")
zz <- file(sinkfile, open = "wt"); sink(zz); on.exit({sink(); close(zz)}, add = TRUE)
print(summary(m))
cat(sprintf("\nAICc = %.2f\n", MuMIn::AICc(m)))
print(performance::r2(m))
sink(); close(zz); on.exit(NULL, add = FALSE)

cat("\nResumo também salvo em:", sinkfile, "\n")
print(summary(m))
cat(sprintf("\nAICc = %.2f\n", MuMIn::AICc(m)))
print(performance::r2(m))

# 6) Odds ratios ------------------------------------
fx <- broom.mixed::tidy(m, effects="fixed", conf.int=TRUE)

# Distance (if present)
if ("distance_km" %in% all.vars(FORM)) {
  rr <- fx[fx$term=="distance_km", ]
  if (nrow(rr)==1) {
    cat(sprintf("\nDistance OR per 1 km = %.3f (%.3f–%.3f)\n",
                exp(rr$estimate), exp(rr$conf.low), exp(rr$conf.high)))
    cat(sprintf("Distance OR per 0.5 km = %.3f (%.3f–%.3f)\n",
                exp(0.5*rr$estimate), exp(0.5*rr$conf.low), exp(0.5*rr$conf.high)))
  }
}

# LULC indices (OR per +0.1) 
in_model  <- fx$term[fx$term %in% all.vars(FORM)]
idx_terms <- in_model[grepl("NDVI|NDWI|NDBI|MR", in_model, ignore.case = TRUE)]

if (length(idx_terms)) {
  cat("\nIndex effects (OR per +0.1):\n")
  or_table <- lapply(idx_terms, function(t) {
    r <- fx[fx$term==t,]
    tibble(term = t,
           OR  = exp(r$estimate*0.1),
           LCL = exp(r$conf.low*0.1),
           UCL = exp(r$conf.high*0.1))
  }) %>% bind_rows()
  print(or_table)
  write_csv(or_table, file.path(OUTDIR,"OR_index_plus0.1.csv"))
}

# Forest plot de ORs  --------------------------------
# Choose the scale for the LULC indices: "sd" (recommended) or "0.1"
OR_MODE <- "sd"  # troque para "0.1" se quiser por +0.1 na escala 0–1

# SDs of the numeric variables used (for OR_MODE == "sd")
num_sd <- sapply(dsub, function(v) if (is.numeric(v)) stats::sd(v, na.rm = TRUE) else NA_real_)

fx_plot <- fx %>%
  filter(term != "(Intercept)") %>%
  mutate(
    step = case_when(
      term == "distance_km" ~ 1,                          # distance by +1 km
      OR_MODE == "0.1"      ~ 0.1,                        # indexes by +0.1
      OR_MODE == "sd"       ~ as.numeric(num_sd[term]),   # indexes by +1 SD
      TRUE ~ NA_real_
    ),
    OR  = exp(estimate * step),
    LCL = exp(conf.low * step),
    UCL = exp(conf.high * step),
    label = case_when(
      term == "distance_km" ~ "Distance (per +1 km)",
      TRUE ~ str_replace_all(term, "_", " ")
    )
  )

subt <- if (OR_MODE == "sd") {
  "Distance: per +1 km | Indices: per +1 SD (z-score)"
} else {
  "Distance: per +1 km | Indices: per +0.1 (on 0–1 scale)"
}

p_forest <- ggplot(fx_plot, aes(y = reorder(label, OR), x = OR)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_point() +
  geom_errorbar(aes(xmin = LCL, xmax = UCL), height = .12) +
  scale_x_log10() +
  labs(x = "Odds ratio (95% CI) [log scale]",
       y = NULL,
       title = "Fixed effects on return probability",
       subtitle = subt) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUTDIR,"forest_OR_fixed_effects.png"), p_forest, width=7, height=5, dpi=300)
print(p_forest)

# FORECASTS AT KEY DISTANCES (LULC SCENARIOS) ------------------------
candidate_vars <- intersect(all.vars(FORM), names(df))
candidate_vars <- candidate_vars[grepl("NDVI|NDWI|NDBI|MR", candidate_vars, ignore.case = TRUE)]
candidate_vars <- candidate_vars[sapply(df[candidate_vars], is.numeric)]

if (length(candidate_vars)) {
  qs <- tibble(
    var = candidate_vars,
    p25 = sapply(candidate_vars, function(v) stats::quantile(df[[v]], 0.25, na.rm=TRUE)),
    p75 = sapply(candidate_vars, function(v) stats::quantile(df[[v]], 0.75, na.rm=TRUE))
  )
  
  mk_scen <- function(var, p25, p75) {
    list(
      ggpredict(m, terms = c(
        "distance_km [0:7.5 by=0.25]",
        paste0(var, " [", round(p25, 3), "]")
      )) |> as.data.frame() |> mutate(scenario = paste(var, "low (P25)")),
      ggpredict(m, terms = c(
        "distance_km [0:7.5 by=0.25]",
        paste0(var, " [", round(p75, 3), "]")
      )) |> as.data.frame() |> mutate(scenario = paste(var, "high (P75)"))
    )
  }
  
  pred_list <- purrr::pmap(qs, mk_scen) |> unlist(recursive = FALSE)
  pred_df   <- dplyr::bind_rows(pred_list)
  
  keyd <- c(0,2,4,6,7.5)
  tab_preds_all <- pred_df %>%
    dplyr::filter(x %in% keyd) %>%
    dplyr::transmute(
      index = sub(" .*","", scenario),
      scenario,
      distance_km = x,
      prob = predicted,
      lo = conf.low,
      hi = conf.high
    ) %>%
    dplyr::arrange(index, scenario, distance_km)
  
  write_csv(tab_preds_all, file.path(OUTDIR, "predictions_by_distance_scenarios.csv"))
  print(tab_preds_all)
  
} else {
  gg_dist <- ggeffects::ggpredict(m, terms = c("distance_km [0:7.5 by=0.25]")) |> as.data.frame()
  keyd <- c(0,2,4,6,7.5)
  tab_preds_all <- gg_dist %>%
    dplyr::filter(x %in% keyd) %>%
    dplyr::transmute(distance_km = x, prob = predicted, lo = conf.low, hi = conf.high)
  
  write_csv(tab_preds_all, file.path(OUTDIR, "predictions_by_distance_only.csv"))
  print(tab_preds_all)
}

# 10) DHARMa Diagnosis ----------------------------------------
png(file.path(OUTDIR,"DHARMa_residuals.png"), width=900, height=700, res=120)
plot(DHARMa::simulateResiduals(m, n = 1000))
dev.off()

cat("\nArquivos salvos em:", OUTDIR, "\n")


# ===============================================================
# TITLE: Cox Proportional Hazards Regression with Frailty for Melipona quadrifasciata 
#        Return Time Analysis
# AUTHOR: Toppa RH, et al.
# AFFILIATION: Federal University of São Carlos, Campus Sorocaba
# CONTACT: toppa@ufscar.br
# DESCRIPTION:
# This script performs survival analysis using Cox regression with frailty terms
# to account for clustered data (colony × release_batch × release_point).
# Includes robust standard errors, model diagnostics, hazard ratios, and
# Kaplan-Meier curves with 35-hour censoring for Melipona quadrifasciata homing analysis.
# ===============================================================

# Packages -----------------------------------------------------
req <- c(
  "readxl","dplyr","tidyr","stringr",
  "survival","coxme","broom","broom.mixed",
  "survminer","ggplot2","performance","readr","tibble"
)
new <- req[!(req %in% rownames(installed.packages()))]
if (length(new)) install.packages(new, dependencies = TRUE)

library(readxl); library(dplyr); library(tidyr); library(stringr)
library(survival); library(coxme); library(broom); library(broom.mixed)
library(survminer); library(ggplot2); library(performance); library(readr); library(tibble)

# Adjust the parameters for data input and output--------------------------
xlsx   <- "C:/r_studio/data_statistic.xlsx"
OUTDIR <- "C:/r_studio/cox_regression"

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Read data and minimum typing — **INSERT PERIOD HERE** --------------
raw <- read_excel(xlsx)

# Time/event names (adjust if necessary)
guess_col <- function(cands, data) cands[cands %in% names(data)][1]
time_col  <- guess_col(c("tempo_retorno","time_return_min","ret_time_min","time_min","tempo_min"), raw)
event_col <- guess_col(c("retorno","event_use","returned","event"), raw)
if (is.na(time_col) || is.na(event_col)) stop("Ajuste os nomes de tempo/retorno em 'time_col' e 'event_col'.")

# Required variables (include 'period' here to avoid misalignment)
need_cols <- c(time_col, event_col, "distance_km","colony","release_point","release_batch",
               "MAX_NDVI_CORR","RANGE_NDWI_CORR","MAX_NDVI_BUFFER","period")
need_cols <- unique(need_cols[need_cols %in% names(raw)])

dat <- raw %>%
  dplyr::select(all_of(need_cols)) %>%
  mutate(
    !!time_col    := suppressWarnings(as.numeric(.data[[time_col]])),
    !!event_col   := as.integer(.data[[event_col]]),
    distance_km   = suppressWarnings(as.numeric(distance_km)),
    colony        = factor(colony),
    release_point = factor(release_point),
    # batch como fator estável (data pode estar em Date/char)
    release_batch = factor(as.character(release_batch)),
    period        = if ("period" %in% names(.)) factor(period) else NULL
  )

# Censorship (35 h = 2100 min)
CENSOR_MAX <- 2100
dat <- dat %>%
  mutate(
    time_censor  = pmin(.data[[time_col]], CENSOR_MAX),
    event_censor = ifelse(.data[[time_col]] > CENSOR_MAX | is.na(.data[[time_col]]), 0L, .data[[event_col]])
  ) %>%
  drop_na(time_censor, event_censor, distance_km)

# Single cluster for coxph (robust)
dat <- dat %>% mutate(cluster_id = interaction(colony, release_batch, release_point, drop = TRUE))

# Distance strata (for KM and possible adjustments) --------
dat <- dat %>%
  mutate(
    dist_strata = cut(distance_km,
                      breaks = c(-Inf, 2.0, 4.5, Inf),
                      labels = c("short (<=2 km)", "medium (2.5–4.5 km)", "long (>=5 km)"),
                      right = TRUE)
  )

# Formulas ----------------------------------------------------
SFORM <- as.formula(Surv(time_censor, event_censor) ~ 1)

# coxme (frailty): modelos principal/distância, corredor e buffer
FORM_DIST_COXME <- update(SFORM, . ~ distance_km + (1|colony) + (1|release_batch) + (1|release_point))

PRED_CORR <- c("distance_km","MAX_NDVI_CORR","RANGE_NDWI_CORR")
PRED_CORR <- PRED_CORR[PRED_CORR %in% names(dat)]
FORM_CORR_COXME <- as.formula(
  paste("Surv(time_censor, event_censor) ~",
        paste(PRED_CORR, collapse=" + "),
        "+ (1|colony) + (1|release_batch) + (1|release_point)")
)

PRED_BUF <- c("distance_km","MAX_NDVI_BUFFER")
PRED_BUF <- PRED_BUF[PRED_BUF %in% names(dat)]
FORM_BUF_COXME <- as.formula(
  paste("Surv(time_censor, event_censor) ~",
        paste(PRED_BUF, collapse=" + "),
        "+ (1|colony) + (1|release_batch) + (1|release_point)")
)

# Coxme adjustment (frailty) -------------------------------------
fit_dist <- coxme(FORM_DIST_COXME, data = dat)
fit_corr <- coxme(FORM_CORR_COXME, data = dat)
fit_buff <- coxme(FORM_BUF_COXME,  data = dat)

cat("\n=== coxme: DISTANCE-ONLY ===\n"); print(fit_dist)
cat("\n=== coxme: CORRIDOR ===\n");     print(fit_corr)
cat("\n=== coxme: BUFFER ===\n");       print(fit_buff)

# Coxph parallel models (1 combined cluster) ---------------
make_coxph <- function(terms, data){
  as.formula(
    paste("Surv(time_censor, event_censor) ~",
          paste(terms, collapse = " + "),
          "+ cluster(cluster_id)")
  )
}
if (length(PRED_CORR)==0) PRED_CORR <- "distance_km"
if (length(PRED_BUF)==0)  PRED_BUF  <- "distance_km"

ph_dist <- coxph(make_coxph("distance_km", dat), data = dat, ties = "efron", x=TRUE, y=TRUE)
ph_corr <- coxph(make_coxph(PRED_CORR,      dat), data = dat, ties = "efron", x=TRUE, y=TRUE)
ph_buff <- coxph(make_coxph(PRED_BUF,       dat), data = dat, ties = "efron", x=TRUE, y=TRUE)

cat("\n=== coxph (robusto): DISTANCE-ONLY ===\n"); print(summary(ph_dist))
cat("\n=== coxph (robusto): CORRIDOR ===\n");     print(summary(ph_corr))
cat("\n=== coxph (robusto): BUFFER ===\n");       print(summary(ph_buff))

# C-index (concordance) --------------------------------------
cidx <- tibble(
  model = c("DISTANCE","CORRIDOR","BUFFER"),
  concordance = c(summary(ph_dist)$concordance[1],
                  summary(ph_corr)$concordance[1],
                  summary(ph_buff)$concordance[1]),
  se = c(summary(ph_dist)$concordance[2],
         summary(ph_corr)$concordance[2],
         summary(ph_buff)$concordance[2])
)
print(cidx)
write_csv(cidx, file.path(OUTDIR,"concordance_cindex.csv"))

# Proportional Hazards (cox.zph) + grafics -------------------
zph_to_tbl <- function(z){ as.data.frame(z$table) %>% tibble::rownames_to_column("term") }

zph_dist <- cox.zph(ph_dist);  cat("\nPH (DIST):\n"); print(zph_dist)
zph_corr <- cox.zph(ph_corr);  cat("\nPH (CORR):\n");  print(zph_corr)
zph_buff <- cox.zph(ph_buff);  cat("\nPH (BUFF):\n");  print(zph_buff)

# Save PH tables
write_csv(zph_to_tbl(zph_dist), file.path(OUTDIR,"PH_test_DIST.csv"))
write_csv(zph_to_tbl(zph_corr), file.path(OUTDIR,"PH_test_CORR.csv"))
write_csv(zph_to_tbl(zph_buff), file.path(OUTDIR,"PH_test_BUFF.csv"))

# graphs in R (ggcoxzph) and PNG base
print(ggcoxzph(zph_dist))
print(ggcoxzph(zph_corr))
print(ggcoxzph(zph_buff))

png(file.path(OUTDIR,"PH_schoenfeld_DIST.png"), 900, 700, res=120);  plot(zph_dist); dev.off()
png(file.path(OUTDIR,"PH_schoenfeld_CORR.png"), 900, 700, res=120);  plot(zph_corr); dev.off()
png(file.path(OUTDIR,"PH_schoenfeld_BUFF.png"), 900, 700, res=120);  plot(zph_buff); dev.off()

# HR (Hazard Ratios) with IC95% (coxph) -------------------
tidy_hr <- function(fit){
  broom::tidy(fit, conf.int=TRUE, exponentiate = TRUE) |>
    dplyr::filter(!grepl("^cluster\\(", term))
}
hr_dist <- tidy_hr(ph_dist); cat("\nHR DISTANCE:\n"); print(hr_dist)
hr_corr <- tidy_hr(ph_corr); cat("\nHR CORRIDOR:\n"); print(hr_corr)
hr_buff <- tidy_hr(ph_buff); cat("\nHR BUFFER:\n");   print(hr_buff)

write_csv(hr_dist, file.path(OUTDIR,"HR_distance.csv"))
write_csv(hr_corr, file.path(OUTDIR,"HR_corridor.csv"))
write_csv(hr_buff, file.path(OUTDIR,"HR_buffer.csv"))

# KM curves with numbers-at-risk ---------------------------------------
km_dist <- survfit(Surv(time_censor, event_censor) ~ dist_strata, data = dat)
p_km <- ggsurvplot(
  km_dist, data = dat,
  risk.table = TRUE, conf.int = TRUE, pval = TRUE,
  xlab = "Time since release (min)", ylab = "Return probability (survival)",
  legend.title = "Distance strata",
  ggtheme = theme_minimal(base_size = 12)
)
print(p_km)  # show in R
ggsave(file.path(OUTDIR,"KM_by_distance.png"), p_km$plot,  width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTDIR,"KM_by_distance_risktable.png"), p_km$table, width = 7, height = 2.8, dpi = 300)

# (Optional) per period, if applicable  -------------------------------------
if ("period" %in% names(dat)) {
  km_per <- survfit(Surv(time_censor, event_censor) ~ period, data = dat)
  p_km2 <- ggsurvplot(
    km_per, data = dat,
    risk.table = TRUE, conf.int = TRUE, pval = TRUE,
    xlab = "Time since release (min)", ylab = "Return probability (survival)",
    legend.title = "Period",
    ggtheme = theme_minimal(base_size = 12)
  )
  print(p_km2)
  ggsave(file.path(OUTDIR,"KM_by_period.png"), p_km2$plot,  width = 7, height = 5, dpi = 300)
  ggsave(file.path(OUTDIR,"KM_by_period_risktable.png"), p_km2$table, width = 7, height = 2.8, dpi = 300)
}

# Summary table ----------------------------------------
summ_tbl <- bind_rows(
  hr_dist %>% mutate(model = "DISTANCE"),
  hr_corr %>% mutate(model = "CORRIDOR"),
  hr_buff %>% mutate(model = "BUFFER")
) %>%
  select(model, term, estimate, conf.low, conf.high, p.value) %>%
  rename(HR = estimate, LCL = conf.low, UCL = conf.high, p = p.value) %>%
  left_join(cidx, by = join_by(model == model))

print(summ_tbl)
write_csv(summ_tbl, file.path(OUTDIR,"cox_summary_for_MS.csv"))

# Functional form and influence -------------------
# Martingale residuals (DISTANCE model)
mres <- residuals(ph_dist, type = "martingale")
p_mres <- ggplot(data.frame(distance_km = dat$distance_km, mres),
                 aes(distance_km, mres)) +
  geom_point(alpha=.35) +
  geom_smooth(se=FALSE, method="loess") +
  theme_minimal(base_size = 12) +
  labs(x="Distance (km)", y="Martingale residuals",
       title="Functional form check (DISTANCE model)")
print(p_mres)
ggsave(file.path(OUTDIR,"martingale_distance.png"), p_mres, width=7, height=5, dpi=300)

# Influence (dfbeta) – save summary
dfb <- residuals(ph_dist, type = "dfbeta")
dfb_sum <- as.data.frame(dfb) %>% summarise(across(everything(), ~max(abs(.), na.rm=TRUE)))
write_csv(dfb_sum, file.path(OUTDIR,"influence_dfbeta_DIST.csv"))

# Summary of events/censorship
evt_tbl <- dat %>%
  summarise(n = n(),
            events = sum(event_censor==1, na.rm=TRUE),
            censored = sum(event_censor==0, na.rm=TRUE),
            median_followup_min = median(time_censor, na.rm=TRUE))
print(evt_tbl)
write_csv(evt_tbl, file.path(OUTDIR,"events_censoring_summary.csv"))

cat("\nArquivos gerados em:\n", OUTDIR, "\n")



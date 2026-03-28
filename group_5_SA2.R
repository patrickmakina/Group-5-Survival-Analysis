# Survival Analysis Group Project 2: Advanced Survival Analysis
# Acute Myelogenous Leukemia (AML) Remission Data
# Group 5: Kutlo, Patrick, Sylvestre, Oliver, Faida
# Project 2

# Clean environment
rm(list = ls())

# Load required packages
library(survival)
library(flexsurv)
library(knitr)
library(kableExtra)
library(dplyr)
library(muhaz)

# STEP 0: Data Preparation from Project 1

data(cancer)
str(myeloid)

# Recode variables
myeloid = myeloid %>%
  mutate(
    trt = factor(trt, levels = c("A", "B"),
                 labels = c("Standard", "Experimental")),
    sex = factor(sex, levels = c("f", "m"),
                 labels = c("Female", "Male")),
    flt3 = factor(flt3, levels = c("A", "B", "C"),
                  labels = c("Low", "Moderate", "High")))

str(myeloid)
summary(myeloid)

# Create survival object using full dataset with covariates
SurvObj = Surv(time = myeloid$futime, event = myeloid$death)

# STEP 1: Non-Parametric Analysis (Kaplan-Meier Estimator)

# Overall KM estimate
km_overall = survfit(SurvObj ~ 1, data = myeloid)
summary(km_overall)

# Median survival time and 95% CI
km_overall

# Survival probabilities at 25th, 50th, 75th percentiles of observed time
quantile_times = quantile(myeloid$futime, probs = c(0.25, 0.50, 0.75))
km_summary = summary(km_overall, times = quantile_times)
km_prob_table = data.frame(
  Percentile    = c("25th", "50th", "75th"),
  Time_days     = round(quantile_times, 2),
  Survival_Prob = round(km_summary$surv,  4),
  Lower_95CI    = round(km_summary$lower, 4),
  Upper_95CI    = round(km_summary$upper, 4))

kable(km_prob_table, row.names = FALSE,
      col.names = c("Percentile", "Time (days)",
                    "Survival Probability S(t)",
                    "Lower 95% CI", "Upper 95% CI"),
      caption = "<span style='color:black; font-weight:bold;'>Kaplan-Meier Survival Probabilities at Key Time Percentiles</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  row_spec(0, color = "black")

# Overall KM plot
par(mfrow = c(1, 1))
plot(km_overall,
     xlab = "Time (days)", ylab = "Survival Probability S(t)",
     main = "Kaplan-Meier Overall Survival Estimate",
     col = c("black", "blue", "blue"), lwd = 2, conf.int = TRUE,
     lty = c(1, 2, 2))
abline(h = 0.5, col = "magenta", lwd = 2.5, lty = 3)
legend("topright",
       legend = c("KM Estimate", "95% CI", "50% survival"),
       col = c("black", "blue", "magenta"),
       lty = c(1, 2, 3), lwd = 2, cex = 0.85)

# Stratified KM by treatment
km_trt = survfit(SurvObj ~ trt, data = myeloid)
km_trt
par(mfrow = c(1, 1))
plot(km_trt,
     xlab = "Time (days)", ylab = "Survival Probability S(t)",
     main = "Kaplan-Meier Curves by Treatment",
     col = c("brown", "purple"), lwd = 2, conf.int = FALSE)
abline(h = 0.5, col = "magenta", lwd =2.5, lty = 3)
legend("topright",
       legend = c("Standard", "Experimental", "50% survival"),
       col = c("brown", "purple", "magenta"),
       lty = c(1, 1, 3), lwd = 2, cex = 0.85)

# Stratified KM by sex
km_sex = survfit(SurvObj ~ sex, data = myeloid)
km_sex
par(mfrow = c(1, 1))
plot(km_sex,
     xlab = "Time (days)", ylab = "Survival Probability S(t)",
     main = "Kaplan-Meier Curves by Sex",
     col = c("turquoise", "darkorange"), lwd = 2, conf.int = FALSE)
abline(h = 0.5, col = "magenta", lwd = 2.5, lty = 3)
legend("topright",
       legend = c("Female", "Male", "50% survival"),
       col = c("turquoise", "darkorange", "magenta"),
       lty = c(1, 1, 3), lwd = 2, cex = 0.85)

# Stratified KM by FLT3 mutation level
km_flt3 = survfit(SurvObj ~ flt3, data = myeloid)
km_flt3
par(mfrow = c(1, 1))
plot(km_flt3,
     xlab = "Time (days)", ylab = "Survival Probability S(t)",
     main = "Kaplan-Meier Curves by FLT3 Mutation Level",
     col = c("forestgreen", "dodgerblue", "firebrick"), lwd = 2, conf.int = FALSE)
abline(h = 0.5, col = "magenta", lwd = 2.5, lty = 3)
legend("topright",
       legend = c("Low", "Moderate", "High", "50% survival"),
       col = c("forestgreen", "dodgerblue", "firebrick", "magenta"),
       lty = c(1, 1, 1, 3), lwd = 2, cex = 0.85)

# STEP 2: Non-Parametric Analysis (Nelson-Aalen Estimator)

# Nelson-Aalen cumulative hazard estimate
na_fit = survfit(SurvObj ~ 1, data = myeloid, type = "fleming-harrington")

# Cumulative hazard H(t) = -log(S(t))
na_time = na_fit$time
na_cumhaz = na_fit$cumhaz

# Plot cumulative hazard
par(mfrow = c(1, 1))
plot(na_time, na_cumhaz,
     type = "l", col = "darkred", lwd = 2,
     xlab = "Time (days)", ylab = "Cumulative Hazard H(t)",
     main = "Nelson-Aalen Cumulative Hazard Estimate")

# Nelson-Aalen cumulative hazard and derived hazard at key time points
key_times = c(30, 90, 180, 365, 730, 1095)
na_at_keys = summary(na_fit, times = key_times)

# Non-parametric hazard derived from Nelson-Aalen: h(t) = dH(t)/dt
# Approximate using differences
na_h_t = diff(c(0, na_at_keys$cumhaz)) / diff(c(0, key_times))

# Gompertz hazard at same time points

# Fit best model from Project 1 (Gompertz)
best_model = flexsurvreg(SurvObj ~ 1, dist = "gompertz")
time_seq   = seq(1, max(myeloid$futime), by = 1)
gompertz_h_keys = summary(best_model, t = key_times,
                          type = "hazard", ci = FALSE)[[1]]$est
na_table = data.frame(
  Time_days = key_times,
  Cumulative_Hazard = round(na_at_keys$cumhaz, 4),
  NA_Hazard = round(na_h_t, 6),
  Gompertz_Hazard = round(gompertz_h_keys, 6))

kable(na_table, row.names = FALSE,
      col.names = c("Time (days)", "Cumulative Hazard H(t)",
                    "Non-Parametric Hazard h(t)", "Gompertz Hazard h(t)"),
      caption = "<span style='color:black; font-weight:bold;'>Nelson-Aalen Cumulative Hazard and Hazard Estimates vs Gompertz Parametric Hazard</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  row_spec(0, color = "black")

# Non-parametric hazard estimate using kernel smoothing (muhaz)
haz_smooth = muhaz(myeloid$futime, myeloid$death)

par(mfrow = c(1, 1))
plot(haz_smooth,
     col = "darkblue", lwd = 2,
     xlab = "Time (days)", ylab = "Hazard Rate h(t)",
     main = "Non-Parametric Hazard Estimate (Kernel Smoothed)")

# Compare with Gompertz parametric hazard (best model from Project 1)
gompertz_h = summary(best_model, t = time_seq,
                     type = "hazard", ci = FALSE)[[1]]$est

lines(time_seq, gompertz_h, col = "maroon", lwd = 2.5, lty = 2)
legend("topright",
       legend = c("Non-parametric (Kernel)", "Gompertz (Parametric)"),
       col = c("darkblue", "maroon"), lty = c(1, 2), lwd = 2, cex = 0.85)
# STEP 3: Non-Parametric Hypothesis Testing

# Log-Rank and Wilcoxon tests by treatment
logrank_trt = survdiff(SurvObj ~ trt, data = myeloid, rho = 0)
logrank_trt
wilcoxon_trt = survdiff(SurvObj ~ trt, data = myeloid, rho = 1)
wilcoxon_trt

# Log-Rank and Wilcoxon tests by sex
logrank_sex = survdiff(SurvObj ~ sex, data = myeloid, rho = 0)
logrank_sex
wilcoxon_sex = survdiff(SurvObj ~ sex, data = myeloid, rho = 1)
wilcoxon_sex

# Log-Rank and Wilcoxon tests by FLT3 mutation level
logrank_flt3 = survdiff(SurvObj ~ flt3, data = myeloid, rho = 0)
logrank_flt3
wilcoxon_flt3 = survdiff(SurvObj ~ flt3, data = myeloid, rho = 1)
wilcoxon_flt3

# Summary table of test results

fmt_pval = function(p) ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))
test_table = data.frame(
  Covariate = c("Treatment", "Sex", "FLT3 Mutation Level"),
  LogRank_Chi2 = c(round(logrank_trt$chisq, 3),
                   round(logrank_sex$chisq, 3),
                   round(logrank_flt3$chisq, 3)),
  LogRank_p = c(fmt_pval(logrank_trt$pvalue),
                fmt_pval(logrank_sex$pvalue),
                fmt_pval(logrank_flt3$pvalue)),
  Wilcoxon_Chi2 = c(round(wilcoxon_trt$chisq, 3),
                    round(wilcoxon_sex$chisq, 3),
                    round(wilcoxon_flt3$chisq, 3)),
  Wilcoxon_p = c(fmt_pval(wilcoxon_trt$pvalue),
                 fmt_pval(wilcoxon_sex$pvalue),
                 fmt_pval(wilcoxon_flt3$pvalue)))

kable(test_table, row.names = FALSE,
      col.names = c("Covariate", "Chi-Square", "p-value", "Chi-Square", "p-value"),
      caption = "<span style='color:black; font-weight:bold;'>Log-Rank and Wilcoxon Test Results</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  add_header_above(c(" " = 1, "Log-Rank" = 2, "Wilcoxon" = 2)) %>%
  row_spec(0, color = "black")

# STEP 4: Semi-Parametric Analysis (Cox PH Model)

# Univariable Cox models
cox_trt = coxph(SurvObj ~ trt,  data = myeloid)
summary(cox_trt)
cox_sex = coxph(SurvObj ~ sex,  data = myeloid)
summary(cox_sex)
cox_flt3 = coxph(SurvObj ~ flt3, data = myeloid)
summary(cox_flt3)

# Univariable summary table
fmt_pval = function(p) ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))

covariate_labels = c(
  "trtExperimental" = "Treatment: Experimental vs Standard",
  "sexMale" = "Sex: Male vs Female",
  "flt3Moderate" = "FLT3 Mutation: Moderate vs Low",
  "flt3High" = "FLT3 Mutation: High vs Low")

extract_cox = function(model, label) {
  s = summary(model)
  coef_table = s$coefficients
  conf_table = s$conf.int
  raw_names = rownames(coef_table)
  data.frame(
    Covariate = covariate_labels[raw_names],
    HR = round(conf_table[, "exp(coef)"], 3),
    Lower_95 = round(conf_table[, "lower .95"], 3),
    Upper_95 = round(conf_table[, "upper .95"], 3),
    p_value = fmt_pval(coef_table[, "Pr(>|z|)"]))
}

univar_table = rbind(
  extract_cox(cox_trt, "Treatment"),
  extract_cox(cox_sex, "Sex"),
  extract_cox(cox_flt3, "FLT3"))

kable(univar_table, row.names = FALSE,
      col.names = c("Covariate", "Hazard Ratio", "Lower 95% CI", "Upper 95% CI", "p-value"),
      caption = "<span style='color:black; font-weight:bold;'>Univariable Cox PH Model Results</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  row_spec(0, color = "black")

# Multivariable Cox model
cox_multi = coxph(SurvObj ~ trt + sex + flt3, data = myeloid)
summary(cox_multi)
multi_table = extract_cox(cox_multi, "Multivariable")

kable(multi_table, row.names = FALSE,
      col.names = c("Covariate", "Hazard Ratio", "Lower 95% CI", "Upper 95% CI", "p-value"),
      caption = "<span style='color:black; font-weight:bold;'>Multivariable Cox PH Model Results</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  row_spec(0, color = "black")

# Proportional Hazards Assumption (Schoenfeld Residuals)
ph_test = cox.zph(cox_multi)
ph_test

par(mfrow = c(2, 2))
plot(ph_test,
     var = 1, main = "PH Assumption: Treatment")
plot(ph_test,
     var = 2, main = "PH Assumption: Sex")
plot(ph_test,
     var = 3, main = "PH Assumption: FLT3 Mutation Level")

# Influential Observations (dfbeta)
dfbeta_vals = residuals(cox_multi, type = "dfbeta")

dfbeta_labels = c(
  "Treatment: Experimental vs Standard",
  "Sex: Male vs Female",
  "FLT3 Mutation: Moderate vs Low",
  "FLT3 Mutation: High vs Low")

par(mfrow = c(2,2))
for (j in 1:ncol(dfbeta_vals)) {
  plot(myeloid$futime, dfbeta_vals[, j],
       xlab = "Follow-up Time (days)",
       ylab = "dfbeta",
       main = paste("Influential Observations:", dfbeta_labels[j]),
       pch = 20, col = "darkblue",
       cex.main = 0.85)
  abline(h = 0, col = "red", lwd = 2.5, lty = 2)
}


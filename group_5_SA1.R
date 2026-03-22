# Survival Analysis Group Project: Parametric Survival Analysis of 
# Acute Myelogenous Leukemia (AML) Remission Data
# Group 5: Kutlo, Patrick, Sylvestre, Oliver, Faida

# Clean environment
rm(list=ls())

# Load required packages
library(survival)
library(flexsurv)
library(knitr)
library(kableExtra)
library(dplyr)

# Load dataset
data(cancer)

# 1: Data Exploration and 2: Data Preparation
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

# Check distribution of numeric variables
par(mfrow = c(2, 2))

# Follow-up time
h1 = hist(myeloid$futime, plot = FALSE)
d1 = density(myeloid$futime, na.rm = TRUE)
hist(myeloid$futime, probability = TRUE,
     main = "Follow-up time", xlab = "Days", col = "lightblue",
     ylim = c(0, max(h1$density, d1$y)))
lines(d1, lwd = 2)

# Transplant time
x2 = na.omit(myeloid$txtime)
h2 = hist(x2, plot = FALSE)
d2 = density(x2)
hist(x2, probability = TRUE,
     main = "Transplant time", xlab = "Days", col = "orange",
     ylim = c(0, max(h2$density, d2$y)))
lines(d2, lwd = 2)

# Complete remission time
x3 = na.omit(myeloid$crtime)
h3 = hist(x3, plot = FALSE)
d3 = density(x3)
hist(x3, probability = TRUE,
     main = "Complete remission time", xlab = "Days", col = "purple",
     ylim = c(0, max(h3$density, d3$y)))
lines(d3, lwd = 2)

# Relapse time
x4 = na.omit(myeloid$rltime)
h4 = hist(x4, plot = FALSE)
d4 = density(x4)
hist(x4, probability = TRUE,
     main = "Relapse time", xlab = "Days", col = "coral",
     ylim = c(0, max(h4$density, d4$y)))
lines(d4, lwd = 2)

# Baseline characteristics table
baseline_table = data.frame(
  Characteristic = c(
    "Total Participants",
    "Events (deaths)",
    "Censored",
    "Follow-up time (days)",
    "Sex",
    "",
    "Treatment",
    "",
    "FLT3",
    "",
    "",
    "Transplant time (days)",
    "Complete remission time (days)",
    "Relapse time (days)"),
  `Level / Statistic` = c(
    "N (%)",
    "N (%)",
    "N (%)",
    "Median (IQR)",
    "Female",
    "Male",
    "Standard",
    "Experimental",
    "Low",
    "Moderate",
    "High",
    "Median (IQR)",
    "Median (IQR)",
    "Median (IQR)"),
  `Summary` = c(
    paste0(nrow(myeloid), " (100.0%)"),
    paste0(sum(myeloid$death == 1, na.rm = TRUE), " (",
           round(mean(myeloid$death == 1, na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$death == 0, na.rm = TRUE), " (",
           round(mean(myeloid$death == 0, na.rm = TRUE) * 100, 2), "%)"),
    paste0(round(median(myeloid$futime, na.rm = TRUE), 2), " (",
           round(quantile(myeloid$futime, 0.25, na.rm = TRUE), 2), " - ",
           round(quantile(myeloid$futime, 0.75, na.rm = TRUE), 2), ")"),
    paste0(sum(myeloid$sex == "Female", na.rm = TRUE), " (",
           round(mean(myeloid$sex == "Female", na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$sex == "Male", na.rm = TRUE), " (",
           round(mean(myeloid$sex == "Male", na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$trt == "Standard", na.rm = TRUE), " (",
           round(mean(myeloid$trt == "Standard", na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$trt == "Experimental", na.rm = TRUE), " (",
           round(mean(myeloid$trt == "Experimental", na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$flt3 == "Low", na.rm = TRUE), " (",
           round(mean(myeloid$flt3 == "Low", na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$flt3 == "Moderate", na.rm = TRUE), " (",
           round(mean(myeloid$flt3 == "Moderate", na.rm = TRUE) * 100, 2), "%)"),
    paste0(sum(myeloid$flt3 == "High", na.rm = TRUE), " (",
           round(mean(myeloid$flt3 == "High", na.rm = TRUE) * 100, 2), "%)"),
    paste0(round(median(myeloid$txtime, na.rm = TRUE), 2), " (",
           round(quantile(myeloid$txtime, 0.25, na.rm = TRUE), 2), " - ",
           round(quantile(myeloid$txtime, 0.75, na.rm = TRUE), 2), ")"),
    paste0(round(median(myeloid$crtime, na.rm = TRUE), 2), " (",
           round(quantile(myeloid$crtime, 0.25, na.rm = TRUE), 2), " - ",
           round(quantile(myeloid$crtime, 0.75, na.rm = TRUE), 2), ")"),
    paste0(round(median(myeloid$rltime, na.rm = TRUE), 2), " (",
           round(quantile(myeloid$rltime, 0.25, na.rm = TRUE), 2), " - ",
           round(quantile(myeloid$rltime, 0.75, na.rm = TRUE), 2), ")")),
  `Missing, N (%)` = c(
    "0 (0.00%)",
    paste0(sum(is.na(myeloid$death)), " (", round(sum(is.na(myeloid$death)) / nrow(myeloid) * 100, 2), "%)"),
    paste0(sum(is.na(myeloid$death)), " (", round(sum(is.na(myeloid$death)) / nrow(myeloid) * 100, 2), "%)"),
    paste0(sum(is.na(myeloid$futime)), " (", round(sum(is.na(myeloid$futime)) / nrow(myeloid) * 100, 2), "%)"),
    paste0(sum(is.na(myeloid$sex)), " (", round(sum(is.na(myeloid$sex)) / nrow(myeloid) * 100, 2), "%)"),
    "",
    paste0(sum(is.na(myeloid$trt)), " (", round(sum(is.na(myeloid$trt)) / nrow(myeloid) * 100, 2), "%)"),
    "",
    paste0(sum(is.na(myeloid$flt3)), " (", round(sum(is.na(myeloid$flt3)) / nrow(myeloid) * 100, 2), "%)"),
    "",
    "",
    paste0(sum(is.na(myeloid$txtime)), " (", round(sum(is.na(myeloid$txtime)) / nrow(myeloid) * 100, 2), "%)"),
    paste0(sum(is.na(myeloid$crtime)), " (", round(sum(is.na(myeloid$crtime)) / nrow(myeloid) * 100, 2), "%)"),
    paste0(sum(is.na(myeloid$rltime)), " (", round(sum(is.na(myeloid$rltime)) / nrow(myeloid) * 100, 2), "%)")),
  `Min - Max` = c(
    "",
    "",
    "",
    paste0(round(min(myeloid$futime, na.rm = TRUE), 2), " - ",
           round(max(myeloid$futime, na.rm = TRUE), 2)),
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    paste0(round(min(myeloid$txtime, na.rm = TRUE), 2), " - ",
           round(max(myeloid$txtime, na.rm = TRUE), 2)),
    paste0(round(min(myeloid$crtime, na.rm = TRUE), 2), " - ",
           round(max(myeloid$crtime, na.rm = TRUE), 2)),
    paste0(round(min(myeloid$rltime, na.rm = TRUE), 2), " - ",
           round(max(myeloid$rltime, na.rm = TRUE), 2))),
  check.names = FALSE)

kable(baseline_table, row.names = FALSE,
      caption = "<span style='color:black; font-weight:bold;'>Baseline Characteristics and Survival Outcome Summary</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  row_spec(0, color = "black")

# Remove missing values on time and event
myeloid = myeloid[complete.cases(myeloid[, c("futime", "death")]), ]

# Select only relevant columns
myeloid_clean = myeloid[, c("futime", "death")]

# Create survival object
SurvObj_myeloid = Surv(time = myeloid_clean$futime, event = myeloid_clean$death)
# 3: Fit Multiple Parametric Models
fits = list(
  Exponential = flexsurvreg(SurvObj_myeloid ~ 1, dist = "exp"),
  Weibull     = flexsurvreg(SurvObj_myeloid ~ 1, dist = "weibull"),
  LogNormal   = flexsurvreg(SurvObj_myeloid ~ 1, dist = "lnorm"),
  LogLogistic = flexsurvreg(SurvObj_myeloid ~ 1, dist = "llogis"),
  Gamma       = flexsurvreg(SurvObj_myeloid ~ 1, dist = "gamma"),
  GenGamma    = flexsurvreg(SurvObj_myeloid ~ 1, dist = "gengamma"),
  Gompertz    = flexsurvreg(SurvObj_myeloid ~ 1, dist = "gompertz"))

# Table of MLE estimates, log-likelihood, AIC
get_param = function(fit, i) {
  if (nrow(fit$res) >= i) {
    paste0(rownames(fit$res)[i], " = ", format(round(fit$res[i, "est"], 4), scientific = FALSE))
  } else {
    NA
  }
}

model_table = data.frame(
  Model = names(fits),
  `Parameter 1 Estimate` = sapply(fits, function(f) get_param(f, 1)),
  `Parameter 2 Estimate` = sapply(fits, function(f) get_param(f, 2)),
  `Parameter 3 Estimate` = sapply(fits, function(f) get_param(f, 3)),
  LogLik = sapply(fits, function(f) round(f$loglik, 3)),
  AIC = sapply(fits, function(f) round(AIC(f), 3)),
  row.names = NULL,
  check.names = FALSE)

kable(model_table, row.names = FALSE,
      caption = "<span style='color:black; font-weight:bold;'>Parametric Model Comparison (MLE, LogLik, AIC)</span>",
      escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  row_spec(0, color = "black")

# 4: Model Selection
aic_values = sapply(fits, AIC)
best_model_name = names(which.min(aic_values))
best_model = fits[[best_model_name]]

cat("Best model according to AIC:", best_model_name, "\n")

# Visual comparison: all models vs Kaplan-Meier
par(mfrow = c(1, 1))
plot(survfit(SurvObj_myeloid ~ 1),
     xlab = "Time (days)", ylab = "S(t)",
     main = "All Parametric Models vs Kaplan-Meier",
     col = "black", lwd = 2.5, conf.int = FALSE)

colors = c("pink", "blue", "green", "orange", "purple", "brown", "red")
for (i in seq_along(fits)) {
  lines(fits[[i]], col = colors[i], lwd = 2, ci = FALSE)
}
legend("topright",
       legend = c("Kaplan-Meier", names(fits)),
       col = c("black", colors),
       lwd = 2, cex = 0.75)
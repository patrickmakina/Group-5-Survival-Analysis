# Group 5 Survival Analysis Project
## Project 1: Parametric Survival Analysis of Acute Myelogenous Leukaemia (AML) Remission Data
This project applies parametric survival analysis to the `myeloid` dataset from the `survival` R package, containing data on 646 AML patients enrolled in a randomised clinical trial.

Seven parametric distributions - Exponential, Weibull, Log-Normal, Log-Logistic, Gamma, Generalised Gamma, and Gompertz - were fitted to the survival data using maximum likelihood estimation. Models were compared using AIC, with the Gompertz distribution identified as the best-fitting model.

Life functions of the best model - survival function S(t), hazard function h(t), density function f(t) and cumulative distribution function F(t) - were computed and visualised. Key survival probabilities, mean survival time and variance were also reported and interpreted.

#### Packages used in R: `survival`, `flexsurv`, `knitr`, `kableExtra`, `dplyr`


## Project 2: Advanced Survival Analysis of Acute Myelogenous Leukaemia (AML) Remission Data

This project extends the parametric analysis from Project 1 by applying non-parametric and semi-parametric methods to the myeloid dataset from the survival R package, containing data on 646 AML patients enrolled in a randomised clinical trial.

Non-parametric survival and hazard functions were estimated using the Kaplan-Meier and Nelson-Aalen estimators. Stratified survival curves were generated to compare survival experiences across treatment groups, sex, and FLT3 mutation levels. Log-Rank and Wilcoxon tests were performed to formally assess differences in survival distributions between groups.

A Cox proportional hazards model was fitted to evaluate the simultaneous effects of treatment, sex, and FLT3 mutation level on survival. Univariable and multivariable hazard ratios with 95% confidence intervals were reported. Model diagnostics were conducted to assess the proportional hazards assumption using Schoenfeld residuals and identify influential observations using dfbeta residuals.

Refinement strategies including stratification and time-dependent coefficients were demonstrated to address potential assumption violations. The hazard trend was characterised using the best-fitting parametric model from Project 1, and clinical implications of the final model were interpreted.

### Packages used in R: survival, flexsurv, knitr, kableExtra, dplyr, muhaz

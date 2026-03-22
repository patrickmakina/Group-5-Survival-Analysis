# Group 5 Survival Analysis Project
## Project Title: Parametric Survival Analysis of Acute Myelogenous Leukaemia (AML) Remission Data

This project applies parametric survival analysis to the `myeloid` dataset from the `survival` R package, containing data on 646 AML patients enrolled in a randomised clinical trial.

Seven parametric distributions - Exponential, Weibull, Log-Normal, Log-Logistic, Gamma, Generalised Gamma, and Gompertz - were fitted to the survival data using maximum likelihood estimation. Models were compared using AIC, with the Gompertz distribution identified as the best-fitting model.

Life functions of the best model - survival function S(t), hazard function h(t), density function f(t) and cumulative distribution function F(t) - were computed and visualised. Key survival probabilities, mean survival time and variance were also reported and interpreted.

#### Packages used in R: `survival`, `flexsurv`, `knitr`, `kableExtra`, `dplyr`

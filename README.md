# skmle <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The **skmle** package implements robust and efficient inference for a general class of **transformed hazards models** with sparse, intermittently observed longitudinal covariates. 

Standard survival analysis methods typically require the availability of the entire trajectory of time-dependent covariates. In practice, however, covariates are often observed sparsely and irregularly over time (e.g., occasional clinical visits or vital sign measurements). The `skmle` package resolves this issue using a novel combination of **kernel weighting** and **sieve maximum log-likelihood estimation (SMKLE)**. 

The core estimation computations are implemented via a highly efficient C++ backend utilizing `Rcpp`, `RcppArmadillo`, and the `nlopt` C API, enabling fast and scalable model fitting.

## Features

- **`skmle()`**: Fits flexible transformed hazards models (utilizing the Box-Cox transformation family) using the Sieve Maximum Kernel-weighted Log-likelihood Estimator. This accommodates both the proportional hazards model (`s = 0`) and additive hazards model (`s = 1`), among others.
- **`kee_cox()`**: Fits proportional hazards models using the Kernel Estimating Equations (KEE) approach.
- **`kee_additive()`**: Fits additive hazards models using the KEE approach.

## Installation

You can install the development version of `skmle` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("dayusun/skmle")
```

*(Note: Since the package depends on C++ code, you will need appropriate C++ compilers installed on your system—e.g., Rtools for Windows, or Xcode Command Line Tools for macOS).*

## Usage Example

Below is a brief outline of how to use the `skmle` function to fit a model:

```r
library(skmle)
library(survival)

# Assume `dat` is a long-format dataset containing:
# - id: Subject identifier
# - obs_times: The times at which the covariates were measured
# - time: The observed survival or censoring time
# - status: The event indicator (1 = event, 0 = censored)
# - x1, x2: Time-dependent covariates

# 1. Fit the proportional hazards model (Box-Cox parameter s = 0)
fit_ph <- skmle(Surv(time, status) ~ x1 + x2, 
                data = dat, 
                id = id, 
                obs_times = obs_times,
                s = 0,        # s=0 corresponds to proportional hazards
                h = 0.5,      # Kernel bandwidth
                nknots = 3)   # Number of knots for the B-spline sieve

# Check the summary of the model
summary(fit_ph)

# Plot the estimated nonparametric baseline function
plot(fit_ph)
```

## References

The theoretical foundation and methodology for this package are thoroughly detailed in:

> **Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."** *Journal of the American Statistical Association* (2025): 1-12. [DOI: 10.1080/01621459.2025.2476781](https://doi.org/10.1080/01621459.2025.2476781)

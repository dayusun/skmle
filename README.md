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
- **`sim_skmle_data()`**: A comprehensive set of tools to simulate survival dataset objects equipped with sparse, intermittently observed longitudinal covariate processes (powered by Non-Homogeneous Poisson Processes) compatible for testing the estimators.

## Installation

You can install the development version of `skmle` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("dayusun/skmle")
```

_(Note: Since the package depends on C++ code, you will need appropriate C++ compilers installed on your system—e.g., Rtools for Windows, or Xcode Command Line Tools for macOS)._

## Usage Example

Below is a brief outline of how to use the `sim_skmle_data` function to generate data and the `skmle` function to fit a model:

```r
library(skmle)
library(survival)

# 1. Simulate data for 200 subjects

*Note:* subject identifiers (`id`) may be non‑numeric; the functions will internally convert them to integer codes.
set.seed(123)
dat <- sim_skmle_data(
    n = 200,
    mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
    mu_bar = 8,
    alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
    beta = c(1, -0.5), # True coefficients
    s = 0,             # proportional hazards
    cen = 0.7          # censoring parameter
)

# 2. Fit the proportional hazards model (Box-Cox parameter s = 0)
# 'covariates' column from sim_skmle_data is a matrix containing both covariates
fit_ph <- skmle(Surv(X, delta) ~ covariates,
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

## Benchmarking

We benchmarked the additive and transformed hazards implementations in `skmle` against the equivalent estimating and joint-optimization techniques in the `SurvSparse` package using a simulated dataset of $N=200$ subjects containing sparsely measured longitudinal variables (`h = n^(-0.5)`).

`skmle` drastically outperforms baseline R solvers due to its highly optimized Armadillo C++ backend.

| Model Evaluation                | `SurvSparse` | `skmle` | Speedup         |
| :------------------------------ | :----------- | :------ | :-------------- |
| **Additive Hazards** (`s=1`)    | ~229 ms      | ~80 ms  | **~3x Faster**  |
| **Transformed Hazards** (`s=0`) | ~2.54 s      | ~244 ms | **~10x Faster** |

For the complete breakdown of evaluation iterations and syntax, refer to the package vignette `vignettes("benchmark_survsparse")`.

## References

The theoretical foundation and methodology for this package are thoroughly detailed in:

> **Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."** _Journal of the American Statistical Association_ (2025): 1-12. [DOI: 10.1080/01621459.2025.2476781](https://doi.org/10.1080/01621459.2025.2476781)

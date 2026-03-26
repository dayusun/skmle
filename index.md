# skmle

`skmle` fits transformed hazards survival models when longitudinal
covariates are observed sparsely and intermittently over time. The
package combines kernel weighting and sieve maximum likelihood
estimation to make semiparametric regression feasible in settings where
the full covariate trajectory is not observed.

The core numerical work is implemented in C++ through `Rcpp` and
`RcppArmadillo`, so the user-facing interface follows ordinary R
modeling conventions while the expensive likelihood and
estimating-equation calculations run in compiled code.

The model class is indexed by the Box-Cox transformation parameter `s`.
In practice:

- `s = 0` corresponds to the proportional hazards model.
- `s = 1` corresponds to the additive hazards model.
- other values of `s` interpolate within the transformed hazards family.

This is one of the main motivations for `skmle`: it provides a single
interface for a broader transformed hazards family rather than treating
proportional and additive hazards as completely separate modeling
frameworks.

## Main Functions

- [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md): fit
  the general transformed hazards model by sieve maximum kernel-weighted
  likelihood.
- [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md):
  choose the kernel bandwidth by subject-level cross-validation and
  refit the final model.
- [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md):
  fit the proportional hazards estimating-equation model.
- [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md):
  fit the additive hazards estimating-equation model.
- [`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md):
  simulate sparse longitudinal survival data in long format.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("dayusun/skmle")
```

Because the package compiles C++ code, you need a working toolchain such
as `Rtools` on Windows or the Xcode command line tools on macOS.

## Typical Workflow

``` r
library(skmle)
library(survival)

set.seed(123)

dat <- sim_skmle_data(
  n = 200,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,   # Box-Cox transformation parameter; s = 0 is proportional hazards
  cen = 0.7
)

fit <- skmle(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,   # proportional hazards model
  h = 0.5,
  nknots = 3
)

fit
summary(fit)
plot(fit)
```

The returned objects use a standard R model style:

- [`print()`](https://rdrr.io/r/base/print.html) shows the call and
  fitted coefficients.
- [`summary()`](https://rdrr.io/r/base/summary.html) reports estimates,
  standard errors, z statistics, and p-values.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) is available
  for `skmle` fits to visualize the estimated baseline component.

If you want a specialized estimator instead of the full SMKLE fit:

``` r
fit_cox <- kee_cox(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

summary(fit_cox)
```

Bandwidth selection is available through
[`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md):

``` r
cv_fit <- skmle_cv(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,   # proportional hazards model
  K = 3,
  h_grid = c(0.3, 0.4, 0.5),
  quiet = TRUE
)

cv_fit$h_cv
cv_fit$cv_results
summary(cv_fit$fit)
```

## Benchmarking Against SurvSparse

The package includes a benchmark vignette comparing `skmle` with
`SurvSparse` on matched sparse longitudinal survival-data settings.

In the current benchmark sweep:

| Scenario                       | `SurvSparse` median |          `skmle` median | Relative result     |
|:-------------------------------|--------------------:|------------------------:|:--------------------|
| Additive hazards, `n = 200`    |             ~288 ms |  `kee_additive`: ~87 ms | `skmle` faster      |
| Additive hazards, `n = 500`    |             ~835 ms | `kee_additive`: ~116 ms | `skmle` faster      |
| Transformed hazards, `n = 100` |            ~1096 ms | `skmle(s = 0)`: ~137 ms | `skmle` much faster |
| Transformed hazards, `n = 200` |            ~2238 ms | `skmle(s = 0)`: ~216 ms | `skmle` much faster |

The additive comparison also includes the general spline-based
`skmle(s = 1)` fit, which remains competitive but is slower than
[`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
because it solves the broader joint optimization problem rather than a
specialized estimating equation.

For the full code and benchmark setup, see the vignette:

``` r
vignette("benchmark_survsparse", package = "skmle")
```

For a package tutorial, see:

``` r
vignette("tutorial", package = "skmle")
```

## Box-Cox Transformation in This Package

The transformed hazards model used by `skmle` is indexed by the Box-Cox
parameter `s`. This parameter determines how the hazard-scale model is
linked to the linear predictor and the nonparametric baseline component.

From a user perspective, the key interpretation is:

- choose `s = 0` if you want the proportional hazards model;
- choose `s = 1` if you want the additive hazards model;
- choose another fixed value of `s` if your application calls for an
  intermediate transformed hazards specification.

This gives a coherent way to compare and fit related hazard models
within one estimation framework, using the same long-format sparse
longitudinal data structure.

## References

Sun, D., Sun, Z., Zhao, X., & Cao, H. (2025). Kernel Meets Sieve:
Transformed Hazards Models with Sparse Longitudinal Covariates. *Journal
of the American Statistical Association, 120*(552), 2580-2591.
<https://doi-org.proxy.ulib.uits.iu.edu/10.1080/01621459.2025.2476781>

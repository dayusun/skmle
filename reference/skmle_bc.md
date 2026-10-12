# Bias-Corrected SMKLE by Extrapolation in the Bandwidth

Fits [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) over a
grid of bandwidths, chooses a bandwidth by a plug-in estimate of the
mean squared error of the coefficients, and reports the coefficients
extrapolated to bandwidth zero. The standard error accounts for the
correlation between fits at different bandwidths.

## Usage

``` r
# S3 method for class 'bc'
print(x, ...)

# S3 method for class 'bc'
summary(object, ...)

# S3 method for class 'summary.bc'
print(x, ...)

skmle_bc(
  formula,
  data,
  id,
  obs_times,
  s = 0,
  h_grid = NULL,
  n_h = 12,
  degree = 1,
  n_boot = 2000,
  seed = NULL,
  quiet = FALSE,
  ...
)
```

## Arguments

- x, object:

  A bias-corrected fit from `skmle_bc()`,
  [`kee_cox_bc()`](https://www.sundayu.me/skmle/reference/kee_bc.md) or
  [`kee_async_bc()`](https://www.sundayu.me/skmle/reference/kee_bc.md).

- ...:

  Further arguments passed to
  [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md), such as
  `nknots`, `maxeval` or `optimizer`.

- formula:

  A model formula. The left-hand side must be a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  response and the right-hand side must contain at least one covariate.

- data:

  Data frame containing the variables used in `formula`, `id`, and
  `obs_times`.

- id:

  Subject identifier. Non-numeric identifiers are allowed and are
  internally converted to integer subject codes.

- obs_times:

  Longitudinal observation times aligned row-wise with `data`. Times may
  be on any scale; the sieve basis and the cumulative-hazard quadrature
  are built on the observed follow-up, so there is no need to rescale to
  the unit interval first. `h` must be on the same scale.

- s:

  Box-Cox transformation parameter, defaulting to `0`. `s = 0` is the
  proportional hazards model, `s = 1` is the additive hazards model, and
  values in between interpolate. If you do not have a reason to choose
  otherwise, the default is the familiar Cox model.

- h_grid:

  Optional increasing vector of candidate bandwidths. If `NULL`, `n_h`
  values log-spaced from \\\max\\\min(X_i - T\_{ij} : X_i \> T\_{ij}),\\
  \tau n^{-0.6}\\\\ to \\\tau = \max_i X_i\\.

- n_h:

  Number of candidate bandwidths when `h_grid` is `NULL`.

- degree:

  Degree in \\h\\ of the bias model and of the extrapolation: `1` (the
  default) removes the first-order bias of a one-sided kernel, and `2`
  also removes a term in \\h^2\\.

- n_boot:

  Number of simulated grids used to include in `var` the variability
  from the data-driven choice of bandwidth. With `0`, `var` is the
  variance that treats the selected bandwidth and window as fixed.

- seed:

  Optional integer seed for those simulated grids. If `NULL`, the
  current RNG state is used.

- quiet:

  Logical; if `TRUE`, suppress progress output.

## Value

An object of class `c("skmle_bc", "bc")` with components:

- `coefficients`, `var`: the extrapolated estimate and its covariance,
  including the variance from choosing the bandwidth,

- `var_fixed`: the covariance treating the selected bandwidth and window
  as fixed,

- `h`: the selected bandwidth,

- `fit`: the
  [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) fit at
  `h`,

- `window`: the bandwidths used in the extrapolation,

- `weights`: the matrix of extrapolation weights \\a_j\\, one row per
  grid point and one column per coefficient,

- `grid`: a data frame with one row per bandwidth and coefficient giving
  the estimate, its standard error, the estimated bias and the estimated
  mean squared error,

- `cov_all`: the joint covariance of the estimates over the grid,
  ordered by bandwidth and then coefficient,

- `n_boot`: the number of simulated grids behind `var`,

- `degree`: the degree of the bias model and extrapolation,

- `call`: the matched call.

## The bias

At bandwidth \\h\\ the estimator evaluates the hazard at each time with
covariates observed up to \\h\\ earlier. The lag acts as measurement
error and attenuates the coefficients. For a covariate that changes by
jumps, the resulting bias is linear in \\h\\, and a small-sample term of
order \\1/(nh)\\ is added to it. The mean of \\\hat\beta_h\\ is
therefore modelled as \\\beta_0 + c h + e/h\\.

## The selected bandwidth

For each coefficient, \\c\\ and \\e\\ are estimated by a regression of
\\\hat\beta_h\\ on \\(1, h, 1/h)\\ across the grid, weighted by the
inverse sandwich variances. The selected bandwidth `h` is the grid point
that minimises the estimated mean squared error summed over
coefficients, \\\sum_k \\(\hat c_k h + \hat e_k / h)^2 +
\widehat{\mathrm{var}}(\hat\beta\_{k,h})\\\\. `fit` is the
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) fit at that
bandwidth.

## The extrapolated estimate

`coefficients` is the intercept of a weighted regression of
\\\hat\beta_h\\ on \\(1, h)\\ over the grid points in \\\[h/2, 2h\]\\.
This removes the linear term of the bias. The intercept is a linear
combination \\\sum_j a_j \hat\beta\_{h_j}\\ of the grid estimates, so
with the bandwidth and window held fixed its covariance `var_fixed`
follows from the joint sandwich \$\$\mathrm{Cov}(\hat\beta\_{h_j},
\hat\beta\_{h_m}) = A\_{h_j}^{-1} \Bigl(\tfrac{1}{n}\textstyle\sum_i
S\_{i,h_j} S\_{i,h_m}^\top\Bigr) A\_{h_m}^{-1} / n,\$\$ where
\\S\_{i,h}\\ is subject \\i\\'s score contribution at bandwidth \\h\\.

Because the bandwidth and window are chosen from the same estimates, the
extrapolated estimate has additional variance that `var_fixed` omits.
`var` includes it. `n_boot` grids are drawn from a normal distribution
whose mean is the grid estimates and whose covariance is the joint
sandwich; the selection and extrapolation are rerun on each simulated
grid, and `var` is the covariance of the resulting estimates. No model
is refitted, so this step takes a few seconds. In simulation at \\s =
1\\, the data-driven selection added 6–9\\ estimate.

The extrapolated estimate has smaller bias and larger variance than the
coefficients of `fit`. Use it with `var` for confidence intervals.

A selected bandwidth at an end of the grid raises a warning. Grid points
whose fit has no usable variance are dropped, also with a warning.

## See also

[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md),
[`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)

## Examples

``` r
# \donttest{
library(survival)
set.seed(123)
dat <- sim_skmle_data(
  n = 150,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)
bc <- skmle_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, s = 0, n_h = 8, quiet = TRUE)
bc
#> Call:
#> skmle_bc(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, n_h = 8, quiet = TRUE)
#> 
#> Selected h = 0.2757 (plug-in MSE); extrapolated (degree 1) over h in [0.1794, 0.4236]
#> 
#> Bias-corrected coefficients:
#> covariates1 covariates2 
#>   1.8464607  -0.7798997 
summary(bc)
#> Call:
#> skmle_bc(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, n_h = 8, quiet = TRUE)
#> 
#> Selected h = 0.2757; extrapolated over h in [0.1794, 0.4236]
#> 
#>             estimate std.error      z  p.value estimate_at_h std.error_at_h
#> covariates1   1.8465    0.4855  3.803 0.000143        1.1706         0.2639
#> covariates2  -0.7799    0.4462 -1.748 0.080463       -0.4342         0.3124
# }
```

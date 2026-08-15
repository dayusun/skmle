# Select the Bandwidth by Cross-Validation

Perform K-fold cross-validation to select the kernel bandwidth for
[`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md).

## Usage

``` r
skmle_cv(
  formula,
  data,
  id,
  obs_times,
  s,
  K = 5,
  h_grid = NULL,
  n_h = 10,
  nknots = 3,
  norder = 3,
  lq_nodes = 64,
  maxeval = 10000,
  xtol_rel = 1e-06,
  seed = NULL,
  quiet = FALSE,
  one_sided = TRUE
)

# S3 method for class 'cv.skmle'
print(x, ...)
```

## Arguments

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

  Longitudinal observation times aligned row-wise with `data`.

- s:

  Box-Cox transformation parameter. `s = 0` gives the proportional
  hazards model and `s = 1` gives the additive hazards model.

- K:

  Number of folds.

- h_grid:

  Optional numeric vector of candidate bandwidth values. If `NULL`, a
  grid is generated automatically from the observed time gaps.

- n_h:

  Number of candidate bandwidths to generate when `h_grid` is `NULL`.

- nknots:

  Number of interior knots used in the sieve approximation of the
  baseline component.

- norder:

  Order parameter supplied to the high-level interface for the spline
  approximation.

- lq_nodes:

  Number of Legendre-Gauss quadrature nodes used in numerical
  integration.

- maxeval:

  Maximum number of optimizer evaluations.

- xtol_rel:

  Relative convergence tolerance passed to the optimizer.

- seed:

  Optional integer seed for the random subject-to-fold assignment. If
  `NULL`, the current RNG state is used and no explicit seed is set.

- quiet:

  Logical; if `TRUE`, suppress progress output.

- one_sided:

  Logical. `TRUE` (the default) uses a **half** kernel: only covariate
  observations strictly before the event or quadrature time inform that
  time, which is the risk-set restriction and the estimator as
  published. `FALSE` uses a **full**, two-sided kernel, smoothing the
  covariate path from both sides. The switch applies to the risk-set
  averages inside the C++ backend as well as to the row weights, so the
  two are always consistent.

- x:

  A `cv.skmle` object.

- ...:

  Ignored.

## Value

An object of class `cv.skmle` with components:

- `h_cv`: selected bandwidth,

- `fit`: `skmle` fit refit on the full data,

- `cv_results`: data frame of candidate bandwidths and CV losses,

- `h_grid`: bandwidth grid used in the search,

- `fold_id`: the subject-to-fold assignment vector (length `n`),

- `seed`: the value of `seed` supplied by the user, or `NULL`,

- `call`: the matched call.

## Details

`skmle_cv()` splits subjects, not rows, across folds. This is the
appropriate unit for cross-validation because multiple rows belong to
the same subject in the long-format data structure.

After choosing the bandwidth with the smallest average validation loss,
the function refits
[`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) on the
full data set using the selected value. Because the fold assignment is
random, pass `seed` (or
[`set.seed()`](https://rdrr.io/r/base/Random.html) before calling) to
make the grid selection reproducible.

## Kernel choice

`one_sided` is used inside the fold loop as well as being passed through
to the refit, so the bandwidth is selected under the same kernel the
final fit uses.

## Examples

``` r
# \donttest{
library(survival)

set.seed(123)
dat <- sim_skmle_data(
  n = 60,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)

cv_fit <- skmle_cv(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,
  K = 3,
  h_grid = c(0.3, 0.4, 0.5),
  seed = 2026,
  quiet = TRUE
)

cv_fit$h_cv
#> [1] 0.5
cv_fit$cv_results
#>     h    cvloss
#> 1 0.3 0.8171964
#> 2 0.4 0.7012594
#> 3 0.5 0.6076993
summary(cv_fit$fit)
#> Call:
#> skmle::skmle(formula = Surv(X, delta) ~ covariates, data = dat, 
#>     id = id, obs_times = obs_times, s = 0, h = 0.5)
#> 
#>   n= 60
#> 
#>             Estimate Std. Error z value Pr(>|z|)  
#> covariates1  0.84088    0.35747  2.3523  0.01866 *
#> covariates2 -0.34645    0.41266 -0.8396  0.40115  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -0.1247 
# }
```

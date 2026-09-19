# Select the Bandwidth by Cross-Validation

Perform K-fold cross-validation to select the kernel bandwidth for
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md).

## Usage

``` r
skmle_cv(
  formula,
  data,
  id,
  obs_times,
  s = 0,
  K = 5,
  h_grid = NULL,
  n_h = 10,
  nknots = 3,
  lq_nodes = 64,
  maxeval = 10000,
  xtol_rel = 1e-06,
  seed = NULL,
  quiet = FALSE,
  weight = NULL,
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

  Longitudinal observation times aligned row-wise with `data`. Times may
  be on any scale; the sieve basis and the cumulative-hazard quadrature
  are built on the observed follow-up, so there is no need to rescale to
  the unit interval first. `h` must be on the same scale.

- s:

  Box-Cox transformation parameter, defaulting to `0`. `s = 0` is the
  proportional hazards model, `s = 1` is the additive hazards model, and
  values in between interpolate. If you do not have a reason to choose
  otherwise, the default is the familiar Cox model.

- K:

  Number of folds.

- h_grid:

  Optional numeric vector of candidate bandwidth values. If `NULL`, a
  grid is generated automatically from the observed time gaps.

- n_h:

  Number of candidate bandwidths to generate when `h_grid` is `NULL`.

- nknots:

  Number of interior knots used in the sieve approximation of the
  baseline component. Knots are placed at `(1:nknots)/(nknots + 1)`. The
  basis is a natural cubic spline
  ([`splines::ns`](https://rdrr.io/r/splines/ns.html)); its order is
  fixed, which is why there is no `norder` argument.

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

- weight:

  Weight function of the standardised lag \\(t - r)/h\\, or `NULL` (the
  default) for the Epanechnikov kernel implied by `one_sided`. Any R
  function will do; see
  [kernel-weights](https://www.sundayu.me/skmle/reference/kernel-weights.md).

- one_sided:

  Logical. `TRUE` (the default) uses a half kernel: only covariate
  observations strictly before the event or quadrature time inform that
  time, which is the risk-set restriction and the estimator as
  published. `FALSE` uses a full, two-sided kernel, smoothing the
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

- `cv_results`: data frame of candidate bandwidths, their CV losses, and
  the standard error of each loss across the folds,

- `h_grid`: bandwidth grid used in the search,

- `fold_id`: the subject-to-fold assignment vector (length `n`),

- `seed`: the value of `seed` supplied by the user, or `NULL`,

- `call`: the matched call.

## Kernel choice

`one_sided` is used inside the fold loop as well as being passed through
to the refit, so the bandwidth is selected under the same kernel the
final fit uses.

## How the folds are formed

`skmle_cv()` splits subjects across folds. Several rows belong to the
same subject in long format, so splitting by row would put one subject
on both sides of the split.

After choosing the bandwidth with the smallest average validation loss,
the function refits
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) on the full
data set using the selected value. Because the fold assignment is
random, pass `seed` (or
[`set.seed()`](https://rdrr.io/r/base/Random.html) before calling) to
make the grid selection reproducible.

## What the held-out loss is

Each training fold is fitted at the candidate bandwidth, and the fit is
scored on the held-out subjects by an ordinary log-likelihood per
subject,

\$\$-\frac{1}{n\_{\mathrm{test}}} \sum\_{i \in \mathrm{test}} \left\[
\delta_i \log g\\\hat\alpha(X_i) + Z_i(X_i)^\top \hat\beta\\ -
\int_0^{X_i} g\\\hat\alpha(t) + Z_i(t)^\top \hat\beta\\\\dt
\right\],\$\$

where \\Z_i(t)\\ is the covariate carried forward from the last
observation at or before \\t\\ (and the first observation carried back,
before the first observation time). The integral is exact up to the
Legendre rule applied between consecutive observation times, where the
carried-forward path jumps.

No kernel and no bandwidth appear in that expression, and that is the
point. The kernel-weighted log-likelihood
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) maximises
**cannot** be compared across bandwidths: its weights are \\W(u/h)/h\\,
so the weight each subject contributes falls away as `h` grows, and the
criterion decreases monotonically in `h` whatever the fit is worth.
Scored that way the largest candidate wins every grid on every data set,
and coefficients further from the truth are preferred to coefficients
nearer it. Dividing by the admitted weight, or scoring at one bandwidth
held fixed across the grid, does not rescue it. The carried-forward
likelihood is one yardstick for every candidate, so its approximation
cancels out of the comparison and what is left is the quality of the
fit.

## The default grid

When `h_grid` is `NULL` the grid is log-spaced over \\\[\max\\\min_i
(X_i - T\_{ij})\_+,\\ \tau n^{-0.6}\\,\\ \min\\\max_i \max_j (X_i -
T\_{ij})\_+,\\ \tau n^{-0.3}\\\]\\, with \\\tau = \max_i X_i\\, so it
adapts to the scale of the times on its own.

Always look at `cv_results`. A minimum at an endpoint of the grid raises
a warning: the selected value is then the best of the values offered
rather than a minimum, and the grid should be widened. On the automatic
grid that warning is common, because \\n^{-0.3}\\ is the rate the
asymptotics assume while the finite-sample minimum of the loss
frequently lies above it. Widening `h_grid` by hand shows where the
curve turns.

## How sharp the selection is

Not very, and the `se` column says so: it is the standard error of each
loss across the folds, and over a wide middle range of `h` the losses
sit inside one standard error of each other. Read the curve, not only
`h_cv`.

The criterion scores prediction of the held-out hazard, in which the
baseline \\\hat\alpha\\ can absorb attenuation in \\\hat\beta\\, so it
leans towards more smoothing than the coefficients on their own would
want. Over 10 replicates at \\n = 200\\ on a grid spanning `0.05` to
`0.9`, the mean squared error of \\\hat\beta\\ at the selected bandwidth
was `0.155`, against `0.205` at the largest candidate and `0.065` at the
bandwidth an oracle would have picked. Resist the temptation to correct
the lean by taking the smallest bandwidth within one standard error of
the minimum: that lands in the noisy small-`h` end, and scored `0.205`
over the same replicates – no better than taking the largest candidate.
The rise at the left of the curve is real.

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
  seed = 2026,
  quiet = TRUE
)

cv_fit$h_cv
#> [1] 0.158428
# Read the whole table, not just the selection: a minimum on the edge of the
# grid is a boundary artefact and warns.
cv_fit$cv_results
#> # A tibble: 10 × 3
#>         h cvloss     se
#>     <dbl>  <dbl>  <dbl>
#>  1 0.0969  0.295 0.0577
#>  2 0.110   0.294 0.0644
#>  3 0.124   0.258 0.0976
#>  4 0.140   0.235 0.115 
#>  5 0.158   0.228 0.119 
#>  6 0.179   0.233 0.116 
#>  7 0.203   0.238 0.116 
#>  8 0.229   0.238 0.109 
#>  9 0.259   0.229 0.100 
#> 10 0.293   0.233 0.0928
summary(cv_fit$fit)
#> Call:
#> skmle::skmle(formula = Surv(X, delta) ~ covariates, data = dat, 
#>     id = id, obs_times = obs_times, s = 0, h = 0.158427999136039)
#> 
#>   n= 60
#> 
#>             Estimate Std. Error z value Pr(>|z|)   
#> covariates1  1.52299    0.51344  2.9663 0.003014 **
#> covariates2 -0.36860    0.46092 -0.7997 0.423880   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 0.4295 
# }
```

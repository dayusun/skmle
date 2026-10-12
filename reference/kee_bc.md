# Bias-Corrected Kernel Estimating Equations by Extrapolation in the Bandwidth

Applies the bandwidth selection and bias correction of
[`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md) to
the kernel estimating-equation estimators of Cao and coauthors:
[`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) for the
proportional hazards model with sparse longitudinal covariates, and
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md) for
regression of an asynchronous longitudinal response. Each function fits
its estimator over a grid of bandwidths, chooses a bandwidth by a
plug-in estimate of the mean squared error of the coefficients, and
extrapolates the coefficients to bandwidth zero. The standard error
accounts for the correlation between fits at different bandwidths and
for the data-driven choice of bandwidth. See
[`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md) for
details of the method.

## Usage

``` r
kee_cox_bc(
  formula,
  data,
  id,
  obs_times,
  h_grid = NULL,
  n_h = 12,
  degree = 1,
  n_boot = 2000,
  seed = NULL,
  quiet = FALSE,
  ...
)

kee_async_bc(
  data_y,
  data_x,
  formula,
  id,
  time,
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

- formula, data, id, obs_times:

  As in
  [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md).

- h_grid:

  Optional vector of candidate bandwidths. If `NULL`, `n_h` log-spaced
  values are used. `kee_cox_bc()` uses the default grid of
  [`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md).
  `kee_async_bc()` uses values from \\2\\\mathrm{IQR}\\n^{-0.6}\\ to the
  range of the pooled observation times, where IQR is the interquartile
  range of those times.

- n_h:

  Number of candidate bandwidths when `h_grid` is `NULL`.

- degree:

  Degree in \\h\\ of the bias model and of the extrapolation, as in
  [`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md).
  The default `1` removes the first-order bias of a one-sided kernel,
  and of a two-sided kernel when the covariate path is not smooth. When
  the bias of a two-sided kernel is of order \\h^2\\, linear
  extrapolation leaves a remainder of order \\\hat h^2\\, which is small
  relative to the standard error at the selected bandwidth. `2` removes
  this remainder but gives a larger variance.

- n_boot, seed, quiet:

  As in
  [`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md).

- ...:

  Further arguments passed to
  [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) or
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md),
  such as `one_sided`, `weight` or `link`.

- data_y, data_x, time:

  As in
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md).

## Value

An object of class `c("kee_bc", "bc")` with the components described in
[`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md);
`fit` is the
[`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) or
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md) fit
at the selected bandwidth.

## See also

[`skmle_bc()`](https://www.sundayu.me/skmle/reference/skmle_bc.md),
[`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md),
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md),
[`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)

# Asynchronous longitudinal regression, time-dependent coefficients

Fit \\E\\Y(t) \mid X(t)\\ = g\\X(t)^\top \beta(t)\\\\ from asynchronous
sparse longitudinal data, estimating the coefficient curve \\\beta(t)\\
pointwise. This is the time-dependent coefficient estimator of Cao, Zeng
and Fine (2015), the companion to
[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md).

At a target time \\t\\ the estimating equation weights a response
occasion and a covariate occasion by their separate distances to \\t\\,
\$\$U_n\\\beta(t)\\ = n^{-1} \sum_i \sum_j \sum_k W\_{h_1}(T\_{ij} - t)
W\_{h_2}(S\_{ik} - t) X_i(S\_{ik}) \[ Y_i(T\_{ij}) -
g\\X_i(S\_{ik})^\top \beta(t)\\ \] = 0.\$\$

## Usage

``` r
kee_async_td(
  y_id,
  y_time,
  y,
  x_id,
  x_time,
  X,
  times,
  h,
  h2 = h,
  one_sided = FALSE,
  link = c("identity", "log", "logistic"),
  intercept = TRUE,
  maxit = 50L,
  tol = 1e-08
)

# S3 method for class 'kee_td'
print(x, ...)
```

## Arguments

- y_id, y_time, y:

  Subject identifier, observation time and value for each response
  occasion. Vectors of equal length.

- x_id, x_time:

  Subject identifier and observation time for each covariate occasion.
  Vectors of equal length, matching the rows of `X`.

- X:

  Covariate matrix, one row per covariate occasion.

- times:

  Numeric vector of target times at which to estimate \\\beta(t)\\.

- h:

  Bandwidth for the response side, \\h_1\\.

- h2:

  Bandwidth for the covariate side, \\h_2\\. Defaults to `h`.

- one_sided:

  Logical. `FALSE` (the default) is the full two-sided kernel. `TRUE`
  restricts the **covariate** side to observations strictly before the
  target time, so \\\beta(t)\\ is estimated from covariate values that
  were already available at `t`. The response side is always two-sided:
  it is a local average around `t`, not a filtering step. Both lags are
  measured as "how far in the past", matching
  [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  and the survival estimators.

- link:

  One of `"identity"` (closed form), `"log"` or `"logistic"`
  (Newton-Raphson).

- intercept:

  Logical, prepend a column of ones to `X`.

- maxit, tol:

  Newton-Raphson controls, ignored when `link = "identity"`.

- x:

  A `kee_td` object.

- ...:

  Ignored.

## Value

An object of class `kee_td`:

- coefficients:

  `length(times)` by `p` matrix of \\\hat\beta(t)\\.

- se, z, p.value:

  Pointwise standard errors and Wald statistics, same shape as
  `coefficients`.

- var:

  List of \\p \times p\\ sandwich matrices, one per target time.

- times, h, h2, link, n, call:

  Fit metadata.

- convergence:

  Integer per target time: 0 converged, 1 singular, 2 hit `maxit`, 3 no
  data in the window.

`x`, invisibly.

## Details

Because both time arguments are smoothed, the estimator converges at the
**bivariate** smoothing rate \\(n h_1 h_2)^{1/2}\\, slower than the
\\(nh)^{1/2}\\ of the time-invariant fit and much slower than the usual
varying-coefficient rate available under synchronous sampling. Expect
wide bands and choose bandwidths accordingly.

The product weight factorises over the two occasion indices, so no pair
is ever enumerated: the response side collapses to two per-subject
scalars and the covariate side to one weight per row. Each target time
then costs \\O((n_y + n_x p^2))\\. Fits at successive target times are
warm-started from the previous solution when the link is nonlinear.

Standard errors are the pointwise sandwich at each target time; they are
not simultaneous bands, and no undersmoothing correction is applied, so
the intervals are centred on the smoothed curve rather than on
\\\beta(t)\\.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
for time-invariant coefficients,
[`plot.kee_td()`](https://dayusun.github.io/skmle/reference/plot.kee_td.md).

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
fit <- kee_async_td(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
  times = seq(0.2, 0.8, by = 0.1), h = 0.3
)
coef(fit)
#>      (Intercept)        x
#> [1,]   0.5936535 1.194805
#> [2,]   0.6010658 1.171446
#> [3,]   0.6144281 1.207259
#> [4,]   0.6512632 1.259015
#> [5,]   0.6653495 1.360668
#> [6,]   0.6794972 1.418637
#> [7,]   0.6088444 1.482852
```

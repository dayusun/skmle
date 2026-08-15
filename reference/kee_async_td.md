# Asynchronous longitudinal regression with time-dependent coefficients

Fit \\E\\Y(t) \mid X(t)\\ = g\\X(t)^\top \beta(t)\\\\ from asynchronous
sparse longitudinal data, estimating the coefficient curve \\\beta(t)\\
pointwise. This is the time-dependent coefficient estimator of Cao, Zeng
and Fine (2015), the companion to
[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md).

At a target time \\t\\ a response occasion and a covariate occasion are
weighted by their **separate** distances to \\t\\, \$\$U_n\\\beta(t)\\ =
n^{-1} \sum_i \sum_j \sum_k W\_{h_1}(t - T\_{ij}) W\_{h_2}(t - S\_{ik})
X_i(S\_{ik}) \[ Y_i(T\_{ij}) - g\\X_i(S\_{ik})^\top \beta(t)\\ \] =
0.\$\$ The lag between the two occasions never enters, only their
distances to \\t\\. That is the difference from
[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md),
where the lag is everything.

## Usage

``` r
kee_async_td(
  formula,
  data_y,
  data_x,
  id,
  time,
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

- formula:

  A two-sided formula, `y ~ x1 + x2`. The left-hand side is evaluated in
  `data_y`, the right-hand side in `data_x`.

- data_y:

  Data frame of response occasions: one row per time the response was
  recorded.

- data_x:

  Data frame of covariate occasions: one row per time the covariates
  were recorded. Its time grid need not overlap `data_y`'s at all.

- id:

  Subject identifier, unquoted or as a string. Must name a column
  present in **both** tables. Non-numeric identifiers are allowed.

- time:

  Observation time, unquoted or as a string. Must name a column present
  in **both** tables.

- times:

  Numeric vector of target times at which to estimate \\\beta(t)\\. Keep
  them inside the range of the observed times.

- h:

  Bandwidth for the response side, \\h_1\\.

- h2:

  Bandwidth for the covariate side, \\h_2\\. Defaults to `h`.

- one_sided:

  Logical. `FALSE` (default) is the full two-sided kernel. `TRUE`
  restricts the **covariate** side to observations strictly before the
  target time, so \\\beta(t)\\ uses only covariate values already
  available at `t`. The response side stays two-sided: it is a local
  average around `t`, not a filtering step. Both lags are measured as
  "how far in the past", matching
  [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  and the survival estimators.

- link:

  One of `"identity"` (closed form), `"log"` or `"logistic"`
  (Newton-Raphson).

- intercept:

  Logical, include an intercept.

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

  Pointwise standard errors and Wald statistics, the same shape as
  `coefficients`.

- var:

  List of \\p \times p\\ sandwich matrices, one per target time.

- times, h, h2, one_sided, link, n:

  Fit settings.

- nactive:

  Covariate rows carrying nonzero weight at each target time.

- convergence:

  Per target time: `0` converged, `1` singular, `2` hit `maxit`, `3` no
  data in the window.

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods are
available.

## A slower rate than you may expect

Because both time arguments are smoothed, the estimator converges at the
**bivariate** rate \\(n h_1 h_2)^{1/2}\\. That is slower than the
\\(nh)^{1/2}\\ of the time-invariant fit, and much slower than the usual
varying-coefficient rate available when response and covariate are
observed together. Expect wide bands. Do not read a wiggle in the curve
as structure without checking it against them.

A practical consequence: a sample that gives a comfortable
time-invariant fit can be far too small for a credible curve. If the
bands cover a horizontal line across the whole range, the honest reading
is that the data do not resolve time variation, not that \\\beta(t)\\ is
flat.

## Reading the output

`se` holds pointwise sandwich standard errors at each target time. They
are pointwise, not simultaneous, so a curve leaving the band at one or
two target times is not evidence of anything on its own.

No undersmoothing correction is applied, so the intervals are centred on
the smoothed curve \\E\hat\beta(t)\\ rather than on \\\beta(t)\\. With a
large bandwidth the curve is flattened toward a constant and the bands
do not widen to say so. Refit at a smaller `h` to see how much of the
shape is smoothing.

Target times within a bandwidth of the edge of the observed time range
draw on a one-sided window and are the least reliable part of the curve.
Keep `times` inside the data.

## Cost

The product weight factorises over the two occasion indices, so no pair
is ever enumerated: the response side collapses to two per-subject
scalars and the covariate side to one weight per row. Each target time
costs \\O(n_y + n_x p^2)\\. Fits at successive target times are
warm-started from the previous solution when the link is nonlinear.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md),
[`plot.kee_td()`](https://dayusun.github.io/skmle/reference/plot.kee_td.md)

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
fit <- kee_async_td(y ~ x,
  data_y = d$y, data_x = d$x, id = id, time = time,
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

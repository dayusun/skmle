# Asynchronous longitudinal regression with time-dependent coefficients

Fit \\E\\Y(t) \mid X(t)\\ = g\\X(t)^\top \beta(t)\\\\ from asynchronous
sparse longitudinal data, estimating the coefficient curve \\\beta(t)\\
pointwise. This is the time-dependent coefficient estimator of Cao, Zeng
and Fine (2015), the companion to
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md).

At a target time \\t\\ a response occasion and a covariate occasion are
weighted by their **separate** distances to \\t\\, \$\$U_n\\\beta(t)\\ =
n^{-1} \sum_i \sum_j \sum_k W\_{h_1}(t - T\_{ij}) W\_{h_2}(t - S\_{ik})
X_i(S\_{ik}) \[ Y_i(T\_{ij}) - g\\X_i(S\_{ik})^\top \beta(t)\\ \] =
0.\$\$ The lag between the two occasions never enters, only their
distances to \\t\\. That is the difference from
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md),
where the lag is everything.

## Usage

``` r
kee_async_td(
  data_y,
  data_x,
  formula,
  id,
  time,
  times = NULL,
  h = NULL,
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

- data_y:

  Data frame of response occasions: one row per time the response was
  recorded. First argument, so the function pipes.

- data_x:

  Data frame of covariate occasions: one row per time the covariates
  were recorded. Its time grid need not overlap `data_y`'s at all.

- formula:

  A two-sided formula, `y ~ x1 + x2`. **The left-hand side is looked up
  in `data_y` and the right-hand side in `data_x`**, which is the one
  unusual thing about this interface and follows from the data being in
  two tables: there is no single frame holding both.

- id:

  Subject identifier naming a column present in **both** tables. Give it
  unquoted (`id = subject`), as a string (`id = "subject"`), or embraced
  (`id = {{ col }}`) when calling from inside another function.
  Non-numeric identifiers are allowed.

- time:

  Observation time, naming a column present in **both** tables. Accepts
  the same three forms as `id`.

- times:

  Numeric vector of target times at which to estimate \\\beta(t)\\.
  Defaults to 25 points spanning the 10th to 90th percentile of the
  observed response times, which keeps them inside the data: target
  times near the ends of the range are fitted from a one-sided window
  and are the least reliable part of the curve.

- h:

  Bandwidth for the response side, \\h_1\\. If omitted, a rule-of-thumb
  value is used and reported in a message; see
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md).

- h2:

  Bandwidth for the covariate side, \\h_2\\. Defaults to `h`.

- one_sided:

  Logical. `FALSE` (default) is the full two-sided kernel. `TRUE`
  restricts the **covariate** side to observations strictly before the
  target time, so \\\beta(t)\\ uses only covariate values already
  available at `t`. The response side stays two-sided: it is a local
  average around `t`; restricting it as well would give a filtering
  estimate, which answers a different question. Both lags are measured
  as "how far in the past", matching
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
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
is that the data cannot resolve time variation. That is a different
statement from evidence that \\\beta(t)\\ is flat.

## Reading the output

`se` holds sandwich standard errors at each target time separately. Each
band covers its own target time; read across the whole curve and the
coverage is much worse than the nominal level, so a curve straying
outside at one or two times means little on its own.

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

[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md),
[`plot.kee_td()`](https://www.sundayu.me/skmle/reference/plot.kee_td.md)

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
fit <- d$y |>
  kee_async_td(d$x, y ~ x,
    id = id, time = time,
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
tidy(fit)
#> # A tibble: 14 × 6
#>     time term        estimate std.error statistic  p.value
#>    <dbl> <chr>          <dbl>     <dbl>     <dbl>    <dbl>
#>  1   0.2 (Intercept)    0.594    0.0833      7.12 1.05e-12
#>  2   0.3 (Intercept)    0.601    0.0804      7.48 7.67e-14
#>  3   0.4 (Intercept)    0.614    0.0818      7.51 5.86e-14
#>  4   0.5 (Intercept)    0.651    0.0940      6.93 4.19e-12
#>  5   0.6 (Intercept)    0.665    0.112       5.95 2.74e- 9
#>  6   0.7 (Intercept)    0.679    0.112       6.05 1.41e- 9
#>  7   0.8 (Intercept)    0.609    0.117       5.22 1.77e- 7
#>  8   0.2 x              1.19     0.0817     14.6  2.10e-48
#>  9   0.3 x              1.17     0.0743     15.8  6.11e-56
#> 10   0.4 x              1.21     0.0724     16.7  2.06e-62
#> 11   0.5 x              1.26     0.0871     14.5  2.26e-47
#> 12   0.6 x              1.36     0.105      13.0  1.64e-38
#> 13   0.7 x              1.42     0.103      13.8  3.76e-43
#> 14   0.8 x              1.48     0.0993     14.9  2.19e-50
```

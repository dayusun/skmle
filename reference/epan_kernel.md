# Row weights from a weight function

The one place a weight is turned into the numbers the estimators use, so
there is one description of any given kernel rather than one per call
site. `epan_kernel()` is kept because three fitting functions used to
write the Epanechnikov out inline and older code may still reach for it.

## Usage

``` r
epan_kernel(u)

kernel_weights(lag, h, weight = NULL, one_sided = TRUE)
```

## Arguments

- u:

  Numeric vector or matrix of standardised lags.

- lag:

  Raw time difference `t - r`, a vector or a matrix. Matrices keep their
  shape, which the sieve quadrature relies on.

- h:

  Positive bandwidth.

- weight:

  A weight function of the standardised lag, or `NULL` for the
  Epanechnikov that `one_sided` implies; see
  [kernel-weights](https://www.sundayu.me/skmle/reference/kernel-weights.md).
  `NULL` is the default for the same reason it is the default at every
  entry point: a literal
  [`w_epan_half()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)
  here silently turns `one_sided = FALSE` into a one-sided fit, because
  the support stays `[0, 1]` and the negative lags are zeroed by the
  weight rather than kept. The package's own test suite caught that,
  which is the argument for not having written it twice.

- one_sided:

  Logical. When `TRUE` rows with a non-positive lag receive zero weight,
  which is the risk-set restriction of a hazard model: only covariate
  observations before the time inform it. It is applied after the weight
  has been evaluated, because it is a modelling choice rather than a
  property of the kernel.

## Value

`epan_kernel()` returns \\0.75(1 - u^2)\_+\\, shape preserved.

`kernel_weights()` returns the scaled weights `W(lag/h)/h`.

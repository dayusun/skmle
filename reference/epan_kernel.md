# Epanechnikov kernel and the row weights built from it

The kernel was written out inline in
[`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md),
[`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md) and
[`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md).
Keeping one copy matters now that the half/full choice is a user-facing
argument: three inline copies are three places for the switch to be
forgotten.

## Usage

``` r
epan_kernel(u)

kernel_weights(lag, h, one_sided = TRUE)
```

## Arguments

- u:

  Numeric vector or matrix of standardised lags.

- lag:

  Raw time difference `t - r`, a vector or a matrix. Matrices keep their
  shape, which the sieve quadrature relies on.

- h:

  Positive bandwidth.

- one_sided:

  Logical. When `TRUE` (the default) rows with a non-positive lag
  receive zero weight, which is the risk-set restriction of a hazard
  model: only covariate observations before the time inform it. `FALSE`
  smooths from both sides.

## Value

`epan_kernel()` returns \\0.75(1 - u^2)\_+\\, shape preserved.

`kernel_weights()` returns the scaled weights `W(lag/h)/h`.

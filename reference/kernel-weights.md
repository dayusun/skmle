# Kernel weight functions

A weight function says how much a covariate observation at lag `u` from
a time contributes to the estimating equation at that time, where `u` is
the lag divided by the bandwidth. Every estimator in this package takes
one through its `weight` argument, and any R function of `u` will do:
the Epanechnikov kernels here are a starting point, not the only option.

## Usage

``` r
w_epan_half()

w_epan_full()
```

## Value

A function of one argument, the standardised lag \\u = (t - r)/h\\ as a
vector or a matrix, returning the weight with the shape of its argument
preserved and carrying a `support` attribute.

## Functions

- `w_epan_half()`: The Epanechnikov half kernel \\0.75(1 - u^2)\_+\\ on
  \\\[0, 1\]\\, the default for the survival estimators and the weight
  every published result in this package was computed with.

- `w_epan_full()`: The Epanechnikov full kernel \\0.75(1 - u^2)\_+\\ on
  \\\[-1, 1\]\\, the default for the asynchronous longitudinal
  estimators and the kernel of Cao, Zeng and Fine (2015). Symmetry kills
  the first moment, so it is order 2.

## Writing your own

A weight is a function of one argument, the standardised lag \\u = (t -
r)/h\\, returning a value of the same shape. It is called on vectors and
on matrices, and it must preserve the shape it is given: the sieve
quadrature evaluates it on a matrix of node-by-observation lags. Zero
outside its own support is the function's own responsibility.

    gauss_half <- function(u) exp(-u^2 / 2) * (u >= 0 & u <= 3)
    kee_cox(Surv(X, delta) ~ z, data, id = id, obs_times = r,
            h = 0.3, weight = gauss_half)

Nothing about the function has to be polynomial, nonnegative or
symmetric. A weight that takes negative values is supported throughout:
the estimating equations admit a row on whether its weight is nonzero,
not on whether it is positive.

`support` is the one attribute worth setting, as a length-2 numeric
giving the interval outside which the weight is zero. It is used to
enumerate pairs in the asynchronous estimators, where scanning every
pair would otherwise be quadratic. Without it a conservative default is
assumed and the fit is correct but slower.

## Scaling

The weight is used as \\W(u)/h\\. A weight that integrates to something
other than 1 is not wrong: the normalisation cancels wherever the weight
appears in a ratio, which is everywhere in the risk-set averages, and is
absorbed into `h` elsewhere. `w_epan_half()` integrates to `1/2` for
that reason and is the weight every published result in this package
used.

## The one-sided restriction is separate

`one_sided` is a modelling choice, not a property of the weight. It is
the risk-set restriction of a hazard model: only covariate observations
strictly before a time may inform it. It is applied on the sign of the
lag after the weight has been evaluated, so a two-sided weight under
`one_sided = TRUE` is truncated at zero rather than rejected.

## Examples

``` r
W <- w_epan_half()
W(c(-0.5, 0, 0.5, 1, 1.5))
#> [1] 0.0000 0.7500 0.5625 0.0000 0.0000
attr(W, "support")
#> [1] 0 1

# Any function of the standardised lag is a weight.
tri <- function(u) pmax(1 - abs(u), 0)
tri(c(-0.5, 0, 0.5))
#> [1] 0.5 1.0 0.5
```

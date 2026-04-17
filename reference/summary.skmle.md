# Summary for skmle object

Summary for skmle object

## Usage

``` r
# S3 method for class 'skmle'
summary(object, ...)
```

## Arguments

- object:

  An object of class `skmle`.

- ...:

  Further arguments passed to or from other methods.

## Value

An object of class `summary.skmle` containing the call, a coefficient
table with estimates, standard errors, z-statistics and p-values, the
log-likelihood, convergence status, and the sample size `n`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.5)
summary(fit)
} # }
```

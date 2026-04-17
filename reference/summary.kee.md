# Summary for kee object

Summary for kee object

## Usage

``` r
# S3 method for class 'kee'
summary(object, ...)
```

## Arguments

- object:

  An object of class `kee`.

- ...:

  Further arguments passed to or from other methods.

## Value

An object of class `summary.kee` containing the call, a coefficient
table with estimates, standard errors, z-statistics and p-values, the
convergence code, and the sample size `n`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, h = 0.5)
summary(fit)
} # }
```

# Number of subjects contributing to a fit

The asymptotics of every estimator in this package are in the number of
subjects. Rows within a subject are repeated measurements of one
trajectory and carry far less than one observation's worth of
information, so [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns
the subject count, which is what
[`summary()`](https://rdrr.io/r/base/summary.html) prints as `n`.

## Usage

``` r
# S3 method for class 'skmle'
nobs(object, ...)

# S3 method for class 'kee'
nobs(object, ...)

# S3 method for class 'kee_td'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `skmle`, `kee` or `kee_td` object.

- ...:

  Unused.

## Value

An integer, the number of subjects.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 120)
fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
nobs(fit)
#> [1] 120
```

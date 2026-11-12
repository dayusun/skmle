# Number of subjects contributing to a fit

Returns the number of subjects, which is what
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

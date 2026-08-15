# One-row summary of a fit

One-row summary of a fit

## Usage

``` r
# S3 method for class 'skmle'
glance(x, ...)

# S3 method for class 'kee'
glance(x, ...)

# S3 method for class 'kee_td'
glance(x, ...)
```

## Arguments

- x:

  A fitted `skmle`, `kee` or `kee_td` object.

- ...:

  Unused.

## Value

A one-row tibble. `nobs` is the number of **subjects**, which is the
unit the asymptotics are in. No `AIC` or `logLik` column is reported:
the criterion these estimators optimise is a kernel-weighted
pseudo-likelihood, so an information criterion built from it would not
mean what it usually means.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150)
fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
glance(fit)
#> # A tibble: 1 × 7
#>    nobs nterms     h npair link     one_sided convergence
#>   <int>  <int> <dbl> <int> <chr>    <lgl>           <int>
#> 1   150      2   0.3  1863 identity FALSE               0
```

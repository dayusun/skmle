# Print kee object

Print kee object

## Usage

``` r
# S3 method for class 'kee'
print(x, ...)
```

## Arguments

- x:

  An object of class `kee`.

- ...:

  Further arguments passed to or from other methods.

## Value

`x`, invisibly. Called for its side effect of printing a brief summary
of the call and estimated coefficients.

## Examples

``` r
# \donttest{
library(survival)

set.seed(123)
dat <- sim_skmle_data(
  n = 60,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)

fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, h = 0.5)
print(fit)
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#> Coefficients:
#> covariates1 covariates2 
#>   0.9822679  -0.4326534 
# }
```

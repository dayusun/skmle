# Print skmle object

Print skmle object

## Usage

``` r
# S3 method for class 'skmle'
print(x, ...)
```

## Arguments

- x:

  An object of class `skmle`.

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

fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
print(fit)
#> Call:
#> skmle(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#> Coefficients:
#> covariates1 covariates2 
#>    1.045449   -0.459170 
# }
```

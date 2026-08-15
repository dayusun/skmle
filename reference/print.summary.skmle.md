# Print summary of skmle object

Print summary of skmle object

## Usage

``` r
# S3 method for class 'summary.skmle'
print(x, ...)
```

## Arguments

- x:

  An object of class `summary.skmle`.

- ...:

  Further arguments passed to or from other methods.

## Value

`x`, invisibly. Called for its side effect of printing the summary
table.

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
print(summary(fit))
#> Call:
#> skmle(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#>   n= 60
#> 
#>             Estimate Std. Error z value Pr(>|z|)  
#> covariates1  1.04545    0.41413  2.5245  0.01159 *
#> covariates2 -0.45917    0.44770 -1.0256  0.30507  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -0.1861 
# }
```

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
summary(fit)
#> Call:
#> skmle(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#>   n= 60
#> 
#>             Estimate Std. Error z value Pr(>|z|)  
#> covariates1  0.84088    0.35747  2.3523  0.01866 *
#> covariates2 -0.34645    0.41266 -0.8396  0.40115  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -0.1247 
# }
```

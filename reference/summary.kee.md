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
summary(fit)
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#> Cox-type proportional hazards, kernel estimating equation (half kernel)
#>   n= 60 subjects   bandwidth h = 0.5
#> 
#>             Estimate Std. Error z value Pr(>|z|)  
#> covariates1  0.98227    0.39352  2.4961  0.01256 *
#> covariates2 -0.43265    0.44736 -0.9671  0.33348  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```

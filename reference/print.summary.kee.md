# Print summary of kee object

Print summary of kee object

## Usage

``` r
# S3 method for class 'summary.kee'
print(x, ...)
```

## Arguments

- x:

  An object of class `summary.kee`.

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

fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, h = 0.5)
print(summary(fit))
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#>   n= 60
#> 
#>             Estimate Std. Error z value Pr(>|z|)  
#> covariates1  0.76656    0.33898  2.2613  0.02374 *
#> covariates2 -0.30130    0.41342 -0.7288  0.46612  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```

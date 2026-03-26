# Fit a Cox-Type KEE Model

Fit the proportional hazards model for sparse longitudinal covariate
data using a kernel estimating-equation approach.

## Usage

``` r
kee_cox(formula, data, id, obs_times, h)
```

## Arguments

- formula:

  A model formula with a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  response.

- data:

  Data frame containing all variables used in the fit.

- id:

  Subject identifier aligned row-wise with `data`.

- obs_times:

  Longitudinal observation times aligned row-wise with `data`.

- h:

  Positive kernel bandwidth.

## Value

An object of class `kee` containing coefficient estimates, the estimated
variance-covariance matrix, the estimating-equation matrices,
convergence status, and the original function call.

## Details

`kee_cox()` targets the proportional hazards case without estimating a
nonparametric baseline component. It is therefore a useful specialized
alternative to
[`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) when the
scientific model is Cox-type and the main interest is in the regression
coefficients.

## References

Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets
Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
*Journal of the American Statistical Association* (2025): 1-12.

Cao, Hongyuan, et al. "Inference for Cox models with sparse longitudinal
covariates." *Biometrika* (2015).

## Examples

``` r
library(survival)

set.seed(123)
dat <- sim_skmle_data(
  n = 80,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)

fit_cox <- kee_cox(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

fit_cox
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#> Coefficients:
#> covariates1 covariates2 
#>   0.8583207  -0.4983319 
summary(fit_cox)
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#>   n= 80
#> 
#>             Estimate Std. Error z value Pr(>|z|)   
#> covariates1  0.85832    0.29253  2.9341 0.003345 **
#> covariates2 -0.49833    0.34661 -1.4377 0.150516   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

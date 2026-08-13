# Plot the estimated baseline function for skmle model

Plot the estimated baseline function for skmle model

## Usage

``` r
# S3 method for class 'skmle'
plot(x, t_seq = seq(0, 1, length.out = 100), ...)
```

## Arguments

- x:

  An object of class `skmle`.

- t_seq:

  A numeric vector of time points to evaluate the baseline function.
  Default is `seq(0, 1, length.out = 100)`.

- ...:

  Further arguments passed to or from other methods.

## Value

A `ggplot` object showing the estimated nonparametric baseline function
evaluated on `t_seq`.

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
plot(fit)

# }
```

# Plot a cross-validation curve

The held-out loss against the candidate bandwidths, with the selected
value marked. Worth a glance every time: if the curve is still falling
at an endpoint of the grid, the "selected" bandwidth is where you
stopped looking rather than a minimum, and the grid should be widened.

## Usage

``` r
# S3 method for class 'cv.skmle'
plot(x, ...)
```

## Arguments

- x:

  A `cv.skmle` object.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly. Called for the plot.

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
cv <- skmle_cv(Surv(X, delta) ~ covariates, data = dat, id = id,
  obs_times = obs_times, s = 0, K = 3, seed = 2026, quiet = TRUE)
plot(cv)

# }
```

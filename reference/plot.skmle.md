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
if (FALSE) { # \dontrun{
fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.5)
plot(fit)
} # }
```

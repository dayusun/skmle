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
if (FALSE) { # \dontrun{
fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.5)
print(fit)
} # }
```

# Print kee object

Print kee object

## Usage

``` r
# S3 method for class 'kee'
print(x, ...)
```

## Arguments

- x:

  An object of class `kee`.

- ...:

  Further arguments passed to or from other methods.

## Value

`x`, invisibly. Called for its side effect of printing a brief summary
of the call and estimated coefficients.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, h = 0.5)
print(fit)
} # }
```

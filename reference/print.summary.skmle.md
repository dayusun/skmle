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
if (FALSE) { # \dontrun{
fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.5)
print(summary(fit))
} # }
```

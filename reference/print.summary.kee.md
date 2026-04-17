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
if (FALSE) { # \dontrun{
fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, h = 0.5)
print(summary(fit))
} # }
```

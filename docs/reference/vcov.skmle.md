# Extract the covariance matrix of a fitted model

Extract the covariance matrix of a fitted model

## Usage

``` r
# S3 method for class 'skmle'
vcov(object, ...)

# S3 method for class 'kee'
vcov(object, ...)

# S3 method for class 'kee_td'
vcov(object, time = NULL, ...)
```

## Arguments

- object:

  A fitted `skmle`, `kee` or `kee_td` object.

- ...:

  Unused.

- time:

  For a `kee_td` fit, the target time whose covariance matrix is wanted.
  Matched to the nearest fitted target time. When `NULL` the full list
  of matrices is returned, one per target time.

## Value

For `skmle` and `kee`, the estimated covariance matrix of the regression
coefficients. For `kee_td`, see `time`.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 120)
fit <- kee_async(y ~ x, data_y = d$y, data_x = d$x, id = id, time = time, h = 0.3)
vcov(fit)
#>               (Intercept)             x
#> (Intercept)  0.0099252555 -0.0002930651
#> x           -0.0002930651  0.0082048377
```

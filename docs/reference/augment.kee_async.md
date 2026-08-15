# Add fitted means to the covariate table

Fitted means belong to covariate occasions. The estimating equation
evaluates the link at each observed covariate vector, and the kernel
weight is what ties that vector to a response occasion. So
[`augment()`](https://generics.r-lib.org/reference/augment.html) returns
`data_x` with a `.fitted` column, \\g(X_i(S\_{ik})^\top \hat\beta)\\.

## Usage

``` r
# S3 method for class 'kee_async'
augment(x, data_x, ...)
```

## Arguments

- x:

  A `kee_async` fit.

- data_x:

  The covariate table the fit was made from.

- ...:

  Unused.

## Value

`data_x` as a tibble with a `.fitted` column added.

## Details

There is deliberately no `.resid`. A residual needs a response value at
\\S\_{ik}\\, and asynchronous data has none; any residual reported here
would have to be invented.

The fit does not retain its data, so `data_x` must be supplied.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150)
fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
augment(fit, d$x)
#> # A tibble: 742 × 4
#>       id   time       x .fitted
#>    <int>  <dbl>   <dbl>   <dbl>
#>  1     1 0.0618 -0.864   -0.595
#>  2     1 0.629   0.345    1.12 
#>  3     1 0.661   0.437    1.25 
#>  4     1 0.945   0.253    0.985
#>  5     2 0.0233  0.292    1.04 
#>  6     2 0.477   0.716    1.64 
#>  7     2 0.530   0.827    1.80 
#>  8     2 0.553   0.974    2.00 
#>  9     2 0.732   0.0543   0.704
#> 10     2 0.789  -0.101    0.484
#> # ℹ 732 more rows
```

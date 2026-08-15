# Summarise a fit as a tibble

One row per coefficient for `skmle` and `kee` fits, and one row per
(target time, coefficient) for `kee_td`, so a coefficient curve goes
straight into `ggplot2` without reshaping.

## Usage

``` r
# S3 method for class 'skmle'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'kee'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'kee_td'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted `skmle`, `kee` or `kee_td` object.

- conf.int:

  Logical, add `conf.low` and `conf.high` columns.

- conf.level:

  Confidence level for those columns.

- ...:

  Unused.

## Value

A tibble with columns `term`, `estimate`, `std.error`, `statistic` and
`p.value`, preceded by `time` for a `kee_td` fit, and followed by
`conf.low` and `conf.high` when `conf.int = TRUE`.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150)
fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
tidy(fit, conf.int = TRUE)
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.627    0.0972      6.45 1.13e-10    0.436     0.817
#> 2 x              1.41     0.1000     14.1  2.17e-45    1.22      1.61 
```

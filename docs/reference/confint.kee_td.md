# Wald confidence intervals for a coefficient curve

Pointwise Wald intervals at each target time. They are pointwise, not
simultaneous, and are not corrected for smoothing bias, so they cover
\\E\hat\beta(t)\\ rather than \\\beta(t)\\; see
[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md).

## Usage

``` r
# S3 method for class 'kee_td'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A `kee_td` object.

- parm:

  Coefficients to report, by name or position. Defaults to all.

- level:

  Confidence level.

- ...:

  Unused.

## Value

A data frame with one row per (target time, coefficient) and columns
`time`, `term`, `estimate`, `se`, and the two interval endpoints.

## Details

For `skmle` and `kee` fits no method is needed:
[`stats::confint.default()`](https://rdrr.io/r/stats/confint.html) works
once [`vcov()`](https://rdrr.io/r/stats/vcov.html) is available.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
fit <- kee_async_td(d$y, d$x,
  y ~ x, id = id, time = time,
  times = c(0.3, 0.5, 0.7), h = 0.3
)
confint(fit)
#>   time        term  estimate         se     2.5 %    97.5 %
#> 1  0.3 (Intercept) 0.6010658 0.08040101 0.4434828 0.7586489
#> 2  0.5 (Intercept) 0.6512632 0.09396858 0.4670882 0.8354383
#> 3  0.7 (Intercept) 0.6794972 0.11223274 0.4595251 0.8994694
#> 4  0.3           x 1.1714456 0.07434268 1.0257366 1.3171546
#> 5  0.5           x 1.2590145 0.08708608 1.0883289 1.4297001
#> 6  0.7           x 1.4186374 0.10300965 1.2167422 1.6205326
```

# Plot estimated coefficient curves

Draws \\\hat\beta_j(t)\\ against `t` with pointwise Wald bands, one
panel per coefficient.

## Usage

``` r
# S3 method for class 'kee_td'
plot(x, which = seq_len(ncol(x$coefficients)), level = 0.95, ...)
```

## Arguments

- x:

  A `kee_td` object from
  [`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md).

- which:

  Integer or character vector selecting coefficients to draw. Defaults
  to all.

- level:

  Confidence level for the pointwise bands.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly. Called for the plot.

## Details

The bands are pointwise, not simultaneous: reading them as a confidence
region for the whole curve overstates the evidence. They are also not
corrected for smoothing bias, so they cover \\E\hat\beta(t)\\ rather
than \\\beta(t)\\; a large bandwidth flattens the curve without widening
the band to say so.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
fit <- d$y |>
  kee_async_td(d$x, y ~ x,
    id = id, time = time,
    times = seq(0.2, 0.8, by = 0.05), h = 0.3
  )
plot(fit)
```

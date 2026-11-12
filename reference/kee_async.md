# Asynchronous longitudinal regression with time-invariant coefficients

Regress a sparsely observed longitudinal response on a sparsely observed
longitudinal covariate when the two are measured on **different** time
grids and are never seen together.

Each response occasion \\T\_{ij}\\ is paired with every covariate
occasion \\S\_{ik}\\ of the same subject, and each pair contributes in
proportion to its time separation, so no pair is discarded and no value
is carried forward. The estimator solves the kernel-weighted estimating
equation of Cao, Zeng and Fine (2015), \$\$U_n(\beta) = n^{-1} \sum_i
\sum_j \sum_k W_h(T\_{ij} - S\_{ik}) X_i(S\_{ik}) \[ Y_i(T\_{ij}) -
g\\X_i(S\_{ik})^\top \beta\\ \] = 0,\$\$ with \\W_h(u) = W(u/h)/h\\ and
\\W\\ the Epanechnikov kernel.

## Usage

``` r
kee_async(
  data_y,
  data_x,
  formula,
  id,
  time,
  h = NULL,
  one_sided = FALSE,
  link = c("identity", "log", "logistic"),
  intercept = TRUE,
  maxit = 50L,
  tol = 1e-08
)
```

## Arguments

- data_y:

  Data frame of response occasions: one row per time the response was
  recorded. First argument, so the function pipes.

- data_x:

  Data frame of covariate occasions: one row per time the covariates
  were recorded. Its time grid need not overlap `data_y`'s at all.

- formula:

  A two-sided formula, `y ~ x1 + x2`. **The left-hand side is looked up
  in `data_y` and the right-hand side in `data_x`**, which is the one
  unusual thing about this interface and follows from the data being in
  two tables: there is no single frame holding both.

- id:

  Subject identifier naming a column present in **both** tables. Give it
  unquoted (`id = subject`), as a string (`id = "subject"`), or embraced
  (`id = {{ col }}`) when calling from inside another function.
  Non-numeric identifiers are allowed.

- time:

  Observation time, naming a column present in **both** tables. Accepts
  the same three forms as `id`.

- h:

  Positive bandwidth, on the same scale as `time`. If omitted, a value
  is read off the observation times as a rule of thumb and reported in a
  message.

- one_sided:

  Logical. `FALSE` (default) is the full two-sided kernel of the paper;
  `TRUE` is the half kernel, admitting only covariate observations
  strictly before the response.

- link:

  One of `"identity"` (closed form), `"log"` or `"logistic"`
  (Newton-Raphson).

- intercept:

  Logical, include an intercept.

- maxit, tol:

  Newton-Raphson controls, ignored when `link = "identity"`.

## Value

An object of class `kee` with components

- coefficients:

  Named coefficient vector.

- var:

  Sandwich covariance matrix \\A^{-1} B A^{-1}\\.

- A, Sigma:

  The bread and the meat of the sandwich.

- scores:

  Per-subject score contributions, one row per subject; `Sigma` is
  `crossprod(scores)`.

- n:

  Number of subjects, the unit the asymptotics are in.

- npair:

  Number of weighted (response, covariate) pairs contributing.

- h, one_sided, link:

  Fit settings.

- convergence:

  `0` converged, `1` singular, `2` hit `maxit`.

- call:

  The matched call.

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods are
available.

## Bandwidth

`h` is on the same scale as the observation times. Small `h` admits few
pairs; large `h` admits pairs far apart in time. Fit at a few values of
`h` and compare the estimates. A warning is issued when fewer than 5% of
response occasions have any covariate observation in their window, which
usually means `h` is on the wrong scale.

## Half kernel or full kernel

`one_sided = FALSE`, the default, is the full two-sided kernel of the
paper. `one_sided = TRUE` admits only covariate observations strictly
**before** the response. It uses roughly half the pairs, so it needs a
larger `h` for the same effective sample size.

Observation times are assumed independent of the response and covariate
processes.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
for a coefficient curve \\\beta(t)\\,
[`sim_async_data()`](https://www.sundayu.me/skmle/reference/sim_async_data.md)
to simulate.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150, beta = c(0.5, 1.5))

# The two tables are on different time grids and share only `id`.
head(d$y, 3)
#> # A tibble: 3 × 3
#>      id  time      y
#>   <int> <dbl>  <dbl>
#> 1     1 0.202 -2.95 
#> 2     1 0.573 -0.819
#> 3     1 0.898  0.211
head(d$x, 3)
#> # A tibble: 3 × 3
#>      id   time      x
#>   <int>  <dbl>  <dbl>
#> 1     1 0.0618 -0.864
#> 2     1 0.629   0.345
#> 3     1 0.661   0.437

fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.25)
summary(fit)
#> Call:
#> kee_async(data_y = d$y, data_x = d$x, formula = y ~ x, id = id, 
#>     time = time, h = 0.25)
#> 
#> Asynchronous longitudinal regression (identity link, full kernel)
#>   n= 150 subjects   bandwidth h = 0.25
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.628762   0.098039  6.4134 1.423e-10 ***
#> x           1.438840   0.099817 14.4148 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Optimization status: 0 
confint(fit)
#>                 2.5 %   97.5 %
#> (Intercept) 0.4366094 0.820915
#> x           1.2432022 1.634478
tidy(fit, conf.int = TRUE)
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.629    0.0980      6.41 1.42e-10    0.437     0.821
#> 2 x              1.44     0.0998     14.4  4.18e-47    1.24      1.63 

# Data comes first, so it pipes.
d$y |>
  kee_async(d$x, y ~ x, id = id, time = time, h = 0.25) |>
  glance()
#> # A tibble: 1 × 7
#>    nobs nterms     h npair link     one_sided convergence
#>   <int>  <int> <dbl> <int> <chr>    <lgl>           <int>
#> 1   150      2  0.25  1569 identity FALSE               0

# Half kernel: only covariate observations preceding the response.
kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.25, one_sided = TRUE)
#> Call:
#> kee_async(data_y = d$y, data_x = d$x, formula = y ~ x, id = id, 
#>     time = time, h = 0.25, one_sided = TRUE)
#> 
#> Asynchronous longitudinal regression (identity link, half kernel)
#> subjects: 150   bandwidth h = 0.25
#> 
#> Coefficients:
#> (Intercept)           x 
#>   0.6233627   1.4148202 
```

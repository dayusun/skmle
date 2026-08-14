# Asynchronous longitudinal regression, time-invariant coefficients

Fit a regression of a sparsely observed longitudinal response on a
sparsely observed longitudinal covariate when the two are measured on
**different** time grids. Responses at \\T\_{ij}\\ are paired with the
same subject's covariate observations at \\S\_{ik}\\, and each pair is
weighted by its time separation, so no pair has to be discarded and no
value has to be carried forward.

This is the time-invariant coefficient estimator of Cao, Zeng and Fine
(2015), which solves \$\$U_n(\beta) = n^{-1} \sum_i \sum_j \sum_k
W_h(T\_{ij} - S\_{ik}) X_i(S\_{ik}) \[ Y_i(T\_{ij}) -
g\\X_i(S\_{ik})^\top \beta\\ \] = 0,\$\$ with \\W_h(t) = W(t/h)/h\\. See
[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
for the companion estimator with time-dependent coefficients
\\\beta(t)\\.

## Usage

``` r
kee_async(
  y_id,
  y_time,
  y,
  x_id,
  x_time,
  X,
  h,
  one_sided = FALSE,
  link = c("identity", "log", "logistic"),
  intercept = TRUE,
  maxit = 50L,
  tol = 1e-08
)
```

## Arguments

- y_id, y_time, y:

  Subject identifier, observation time and value for each response
  occasion. Vectors of equal length.

- x_id, x_time:

  Subject identifier and observation time for each covariate occasion.
  Vectors of equal length, matching the rows of `X`.

- X:

  Covariate matrix, one row per covariate occasion.

- h:

  Positive bandwidth.

- one_sided:

  Logical. `FALSE` (the default) is the full two-sided Epanechnikov
  kernel of the paper. `TRUE` is the half kernel: only covariate
  observations strictly before the response contribute.

- link:

  One of `"identity"` (closed form), `"log"` or `"logistic"`
  (Newton-Raphson).

- intercept:

  Logical, prepend a column of ones to `X`.

- maxit, tol:

  Newton-Raphson controls, ignored when `link = "identity"`.

## Value

An object of class `kee` with `coefficients`, the sandwich `var`, the
bread `A`, the meat `Sigma`, `n` (subjects), `npair` (weighted pairs
contributing) and the call.

## Details

The link is evaluated at each **observed** covariate vector and the
weight multiplies the resulting contribution. The covariate is never
smoothed and then substituted, which would be regression calibration and
is not consistent under sparse sampling. The estimator converges at the
univariate smoothing rate \\(nh)^{1/2}\\, not at \\\sqrt n\\.

**Half kernel or full kernel.** `one_sided = FALSE`, the default, is the
full two-sided Epanechnikov kernel of the paper: a response is paired
with covariate observations on either side of it. `one_sided = TRUE`
admits only covariate observations that strictly precede the response,
which is what a causal reading of the covariate path requires and what
the survival estimators in this package do by default.

Cao, Zeng and Fine assume the covariance of the covariate process is
twice differentiable and obtain a bias of order \\h^2\\. Their Section 6
notes that relaxing this to one-sided differentiability, which admits
processes with independent increments such as the Ornstein-Uhlenbeck
process, leaves a bias of order \\h\\ instead. Watch for an estimate
that drifts as `h` grows.

Pairs are enumerated in C++ and then collapsed onto covariate rows
before the coefficients are touched, so the full pair design is never
materialised and the Newton iterations cost \\O(n_x p^2)\\ rather than
\\O(n\_{pair} p^2)\\.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
for time-dependent coefficients.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 100, beta = c(0.5, 1.5))
fit <- kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x, h = 0.25)
summary(fit)
#> Call:
#> kee_async(y_id = d$y$id, y_time = d$y$time, y = d$y$y, x_id = d$x$id, 
#>     x_time = d$x$time, X = d$x$x, h = 0.25)
#> 
#>   n= 100
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 0.572213   0.106594  5.3682 7.955e-08 ***
#> x           1.340047   0.095312 14.0596 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Optimization status: 0 

# Half kernel: only covariate observations that precede the response.
kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
  h = 0.25, one_sided = TRUE
)
#> Call:
#> kee_async(y_id = d$y$id, y_time = d$y$time, y = d$y$y, x_id = d$x$id, 
#>     x_time = d$x$time, X = d$x$x, h = 0.25, one_sided = TRUE)
#> 
#> Coefficients:
#> (Intercept)           x 
#>   0.6221119   1.3433354 
```

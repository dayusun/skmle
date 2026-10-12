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
  tol = 1e-08,
  correct = c("none", "pairs", "gaussian")
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
  message. Use
  [`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
  to choose it from the data.

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

  Newton-Raphson controls, ignored when `link = "identity"` and
  `correct = "none"`.

- correct:

  `"none"` (default): the kernel equation of Cao, Zeng and Fine.
  `"pairs"`: the pair-moment correction, with no model for the
  covariates. `"gaussian"`: the correction for one stationary Gaussian
  covariate. See Details.

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

- h, one_sided, link, correct:

  Fit settings.

- variogram:

  With `correct = "gaussian"`, a list holding the estimated variogram
  `gamma` (a function of the lag), the covariate mean `m` and variance
  `s2`; otherwise `NULL`.

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

## Why not something simpler

Two obvious shortcuts are inconsistent here, which is the reason this
function exists.

*Last value carried forward* replaces \\X_i(T\_{ij})\\ by the most
recent observed value. Under sparse sampling the gap back to that value
does not shrink as \\n\\ grows, so the substituted covariate keeps a
non-vanishing error and the coefficient is attenuated toward zero.

*Regression calibration* smooths the covariate path first and
substitutes the smoothed value. Under sparse sampling the smoothing
window never empties, so the same attenuation appears; for a nonlinear
link there is a Jensen bias as well, because the estimating equation
needs the average of \\g(\cdot)\\ rather than \\g\\ of the average.

The kernel-weighted equation instead evaluates the link at each
**observed** covariate vector, and lets the weight express how much that
observation says about the response occasion.

## Rate and bandwidth

The estimator is consistent and asymptotically normal at the smoothing
rate \\(nh)^{1/2}\\. Standard errors shrink more slowly than in a
parametric fit, and the bandwidth becomes a real modelling choice: too
small and few pairs contribute, too large and the bias grows.
[`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
picks it from the data.

`h` is on the **same scale as the observation times**. Unlike the
survival `h` must be on whatever scale the times are on; nothing else is
assumed about that scale. A warning is issued when fewer than 5% of
response occasions have any covariate observation in their window, which
is the signature of a units mistake.

## Half kernel or full kernel

`one_sided = FALSE`, the default, is the full two-sided kernel of the
paper. `one_sided = TRUE` admits only covariate observations strictly
**before** the response, which is what a causal reading of the covariate
path requires and what the survival estimators in this package do by
default. It roughly halves the number of contributing pairs, so it needs
a larger `h` to reach the same effective sample size.

## Assumptions worth checking

Cao, Zeng and Fine assume the covariance function of the covariate
process is twice differentiable, which gives a bias of order \\h^2\\.
Their Section 6 notes that relaxing this to one-sided differentiability,
which admits processes with independent increments such as the
Ornstein-Uhlenbeck process, leaves a bias of order \\h\\ instead. The
practical symptom is an estimate that drifts steadily as `h` grows, so
fit at several bandwidths and look before trusting one number.

Observation times are assumed independent of the response and covariate
processes. Informative observation times need a different estimator.

## Asynchrony correction

`correct` removes the asynchrony bias instead of shrinking it with
\\h\\.

`correct = "pairs"` uses each subject's own covariate pairs in place of
the unobserved covariate at the response time. For two covariate
occasions \\S\_{ik} \ne S\_{ik'}\\ of the same subject, \\(X_i(S\_{ik}),
X_i(S\_{ik'}))\\ has the distribution of \\(X(S), X(T))\\ at the same
two times, so the equation \$\$\frac{\sum W_h(T\_{ij} - S\_{ik})
X_i(S\_{ik}) Y_i(T\_{ij})} {\sum W_h(T\_{ij} - S\_{ik})} -
\frac{\sum\_{k \ne k'} W_h(S\_{ik'} - S\_{ik}) r(S\_{ik'}) X_i(S\_{ik})
g\\X_i(S\_{ik'})^\top \beta\\} {\sum\_{k \ne k'} W_h(S\_{ik'} - S\_{ik})
r(S\_{ik'})} = 0\$\$ has mean zero at the true \\\beta\\ for every
\\h\\. Here \\r = \lambda_T / \lambda_S\\ is the ratio of the response
and covariate visit intensities, estimated from the observation times
alone by a logistic regression of the occasion type on a natural spline
of time (4 degrees of freedom). No model is assumed for the covariate
process: it need not be Gaussian or stationary, and any number of
covariates, including ones constant within a subject, is allowed. At
least one covariate must vary within a subject. The fit needs subjects
with two or more covariate occasions and assumes visit processes whose
second-order intensity factorises (for example Poisson visits).

`correct = "gaussian"` replaces \\g\\X(S)^\top\beta\\\\ by
\\E\[g\\X(T)^\top\beta\\ \mid x(S)\]\\ for a single covariate that is
second-order stationary with variogram \\\gamma(u) = E\\x(t) - x(t -
u)\\^2\\, mean \\m\\ and variance \\s^2\\, computed as if \\x(T) \mid
x(S)\\ were normal with mean \\m + \rho(u)\\x(S) - m\\\\ and variance
\\s^2\\1 - \rho(u)^2\\\\, \\\rho(u) = 1 - \gamma(u)/(2s^2)\\. The
variogram is estimated by isotonic regression of \\\\x(S\_{ik}) -
x(S\_{ik'})\\^2\\ on \\\|S\_{ik} - S\_{ik'}\|\\. With the identity link
it needs only stationarity; with the log and logistic links it also
needs normality. It has a smaller variance than `"pairs"` when these
assumptions hold and a bias that grows with \\h\\ when they do not.

Both corrected equations are unbiased at every fixed \\h\\ under their
assumptions, so \\h\\ affects only the variance. The standard errors are
the sandwich with the intensity ratio or the variogram treated as known.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
to choose `h`,
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

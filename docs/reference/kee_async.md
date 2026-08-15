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
  formula,
  data_y,
  data_x,
  id,
  time,
  h,
  one_sided = FALSE,
  link = c("identity", "log", "logistic"),
  intercept = TRUE,
  maxit = 50L,
  tol = 1e-08
)
```

## Arguments

- formula:

  A two-sided formula, `y ~ x1 + x2`. The left-hand side is evaluated in
  `data_y`, the right-hand side in `data_x`.

- data_y:

  Data frame of response occasions: one row per time the response was
  recorded.

- data_x:

  Data frame of covariate occasions: one row per time the covariates
  were recorded. Its time grid need not overlap `data_y`'s at all.

- id:

  Subject identifier, unquoted or as a string. Must name a column
  present in **both** tables. Non-numeric identifiers are allowed.

- time:

  Observation time, unquoted or as a string. Must name a column present
  in **both** tables.

- h:

  Positive bandwidth, on the same scale as `time`.

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
rate \\(nh)^{1/2}\\, not at \\\sqrt n\\. The bandwidth therefore
matters: too small and few pairs contribute, too large and the bias
grows. Use
[`kee_async_cv()`](https://dayusun.github.io/skmle/reference/kee_async_cv.md)
rather than guessing.

`h` is on the **same scale as the observation times**. Unlike the
survival estimators there is no requirement that times lie on \\\[0,
1\]\\, but `h` must match whatever scale is used. A warning is issued
when fewer than 5% of response occasions have any covariate observation
in their window, which is the signature of a units mistake.

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

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async_cv()`](https://dayusun.github.io/skmle/reference/kee_async_cv.md)
to choose `h`,
[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
for a coefficient curve \\\beta(t)\\,
[`sim_async_data()`](https://dayusun.github.io/skmle/reference/sim_async_data.md)
to simulate.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150, beta = c(0.5, 1.5))

# The two tables are on different time grids and share only `id`.
head(d$y, 3)
#>   id      time          y
#> 1  1 0.2016819 -2.9454025
#> 2  1 0.5728534 -0.8186868
#> 3  1 0.8983897  0.2109702
head(d$x, 3)
#>   id       time          x
#> 1  1 0.06178627 -0.8642242
#> 2  1 0.62911404  0.3452775
#> 3  1 0.66079779  0.4373871

fit <- kee_async(y ~ x,
  data_y = d$y, data_x = d$x,
  id = id, time = time, h = 0.25
)
summary(fit)
#> Call:
#> kee_async(formula = y ~ x, data_y = d$y, data_x = d$x, id = id, 
#>     time = time, h = 0.25)
#> 
#>   n= 150
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

# Half kernel: only covariate observations preceding the response.
kee_async(y ~ x,
  data_y = d$y, data_x = d$x,
  id = id, time = time, h = 0.25, one_sided = TRUE
)
#> Call:
#> kee_async(formula = y ~ x, data_y = d$y, data_x = d$x, id = id, 
#>     time = time, h = 0.25, one_sided = TRUE)
#> 
#> Coefficients:
#> (Intercept)           x 
#>   0.6233627   1.4148202 
```

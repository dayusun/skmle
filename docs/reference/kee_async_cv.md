# Choose the bandwidth for an asynchronous longitudinal fit

K-fold cross-validation over **subjects** for
[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md).
For each candidate bandwidth the estimator is fitted on the training
subjects and scored on the held-out subjects by the kernel-weighted
squared error of the fitted mean, \$\$ \frac{\sum\_{i \in \mathrm{test}}
\sum_j \sum_k W_h(T\_{ij} - S\_{ik}) \[Y_i(T\_{ij}) -
g\\X_i(S\_{ik})^\top \hat\beta\\\]^2} {\sum\_{i \in \mathrm{test}}
\sum_j \sum_k W_h(T\_{ij} - S\_{ik})}.\$\$

## Usage

``` r
kee_async_cv(
  data_y,
  data_x,
  formula,
  id,
  time,
  h_grid = NULL,
  n_h = 10L,
  K = 5L,
  one_sided = FALSE,
  link = c("identity", "log", "logistic"),
  intercept = TRUE,
  maxit = 50L,
  tol = 1e-08,
  seed = NULL,
  quiet = FALSE
)

# S3 method for class 'cv.kee_async'
print(x, ...)
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

- h_grid:

  Optional numeric vector of candidate bandwidths. When `NULL` a grid is
  generated from the observation times; see Details.

- n_h:

  Number of candidates to generate when `h_grid` is `NULL`.

- K:

  Number of folds. Reduced to the number of subjects if larger.

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

- seed:

  Optional integer seed for the subject-to-fold assignment. Pass it (or
  call [`set.seed()`](https://rdrr.io/r/base/Random.html) first) to make
  the selection reproducible.

- quiet:

  Logical; suppress progress output.

- x:

  A `cv.kee_async` object.

- ...:

  Ignored.

## Value

An object of class `cv.kee_async`:

- h_cv:

  The selected bandwidth.

- fit:

  The
  [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  fit refitted on all subjects at `h_cv`.

- cv_results:

  Data frame of `h`, the mean held-out loss `cvloss`, and `nfold_used`,
  the number of folds that produced a usable fit.

- h_grid, fold_id, seed, call:

  Selection metadata.

## How the folds are formed

Folds are formed over subjects, never over rows. Rows within a subject
come from one trajectory, so splitting by row would put a subject on
both sides of the split and leak the answer into the held-out score.

## Why the loss is normalised

Dividing by the total weight is what makes bandwidths comparable. A
larger `h` admits more pairs and so accumulates a larger unnormalised
error whatever the quality of the fit; without the denominator the
criterion would select the smallest candidate every time.

## Why squared error and not a likelihood

These estimators solve an estimating equation. There is no likelihood,
so nothing to evaluate as a held-out log-likelihood. The same weighted
squared error on the response scale is used for all three links, applied
after `link`.

## The default grid

When `h_grid` is `NULL` the grid is log-spaced over \\\[2 (Q_3 - Q_1)
n^{-0.7},\\ 2 (Q_3 - Q_1) n^{-0.3}\]\\, where \\Q_1, Q_3\\ are quartiles
of the pooled observation times of both tables and \\n\\ is the number
of subjects. This is the rule used by the reference implementation of
Cao, Zeng and Fine's method, and it adapts to the scale of `time` on its
own.

Always look at `cv_results`. If the minimum sits at an endpoint of the
grid a warning is issued, and the grid should be widened by hand: the
selected value is then a boundary artefact rather than a minimum.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150, beta = c(0.5, 1.5))
cv <- d$y |>
  kee_async_cv(d$x, y ~ x,
    id = id, time = time,
    h_grid = c(0.15, 0.25, 0.40), K = 3, seed = 1, quiet = TRUE
  )
#> Warning: the selected bandwidth is at an endpoint of 'h_grid'; widen the grid to check that the minimum is interior
cv
#> Call:
#> kee_async_cv(data_y = d$y, data_x = d$x, formula = y ~ x, id = id, 
#>     time = time, h_grid = c(0.15, 0.25, 0.4), K = 3, seed = 1, 
#>     quiet = TRUE)
#> 
#> 3-fold subject-level cross-validation
#> 
#> # A tibble: 3 × 3
#>       h cvloss nfold_used
#>   <dbl>  <dbl>      <dbl>
#> 1  0.15   1.35          3
#> 2  0.25   1.51          3
#> 3  0.4    1.67          3
#> 
#> Selected h = 0.15
#> 
#> Coefficients at the refit:
#> (Intercept)           x 
#>   0.6103921   1.5031392 
cv$cv_results
#> # A tibble: 3 × 3
#>       h cvloss nfold_used
#>   <dbl>  <dbl>      <dbl>
#> 1  0.15   1.35          3
#> 2  0.25   1.51          3
#> 3  0.4    1.67          3
```

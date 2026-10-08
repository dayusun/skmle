# Choose the bandwidth for an asynchronous longitudinal fit

Select the bandwidth for
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md) and
return the fit at it. The default, `method = "variance"`, selects the
bandwidth that minimises the estimated variance of an
asynchrony-corrected fit (`kee_async(correct = "pairs")` by default).
`method = "sync"` selects it for the same fit by K-fold cross-validation
over subjects, and `method = "plain"` is ordinary cross-validation of
the uncorrected fit.

## Usage

``` r
kee_async_cv(
  data_y,
  data_x,
  formula,
  id,
  time,
  h_grid = NULL,
  n_h = 100L,
  K = 5L,
  repeats = 1L,
  one_sided = FALSE,
  link = c("identity", "log", "logistic"),
  intercept = TRUE,
  maxit = 50L,
  tol = 1e-08,
  seed = NULL,
  quiet = FALSE,
  method = c("variance", "sync", "plain"),
  correct = c("pairs", "gaussian")
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
  built from the data; see Details.

- n_h:

  Number of candidates in the default grid.

- K:

  Number of folds. Reduced to the number of subjects if larger.
  Cross-validation methods only.

- repeats:

  Number of independent random K-fold splits; the held-out loss is
  averaged over them. Cross-validation methods only.

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

- seed:

  Optional integer seed for the subject-to-fold assignment. Pass it (or
  call [`set.seed()`](https://rdrr.io/r/base/Random.html) first) to make
  the selection reproducible. Cross-validation methods only.

- quiet:

  Logical; suppress progress output.

- method:

  `"variance"` (default): the bandwidth minimising the estimated
  variance of the asynchrony-corrected fit. `"sync"`: synchronous-loss
  cross-validation of the corrected fit. `"plain"`: kernel-weighted
  cross-validation of the uncorrected fit. See Details.

- correct:

  The correction used with `method = "variance"` or `"sync"`: `"pairs"`
  (default) or `"gaussian"`; see
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md).

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
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
  fit on all subjects at `h_cv`, with the chosen `correct` unless
  `method = "plain"`.

- cv_results:

  Data frame of `h` and the criterion. For `"variance"`, `trSV`,
  \\\mathrm{tr}\\\hat\Sigma_X \hat V(h)\\\\. For the cross-validation
  methods, `cvloss` (the intercept for `"sync"`, the weighted held-out
  loss for `"plain"`) and `nfold_used`, the number of folds over all
  repeats that produced a usable loss.

- h_grid, fold_id, repeats, method, correct, seed, call:

  Selection metadata; `fold_id` and `repeats` are `NULL` for
  `"variance"`.

## Minimum-variance bandwidth (`method = "variance"`)

The corrected equation has mean zero at the true coefficients for every
\\h\\ (see
[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)),
so the bandwidth affects only the variance and the mean squared error of
\\\hat\beta(h)\\ is its variance. The corrected equation is solved on
all subjects at every candidate, and the selected \\h\\ minimises
\$\$\mathrm{tr}\\\hat\Sigma_X \hat V(h)\\,\$\$ with \\\hat V(h)\\ the
sandwich variance of \\\hat\beta(h)\\ and \\\hat\Sigma_X\\ the second
moment of \\X\\ at the response times: the covariate rows weighted by
the estimated intensity ratio \\\lambda_T / \lambda_S\\ for
`correct = "pairs"`, unweighted for `correct = "gaussian"`. The
criterion is \\n\\ times the expected squared prediction error
\\E\[\\X(T)^\top(\hat\beta - \beta_0)\\^2\]\\, and it does not depend on
the scale of the covariates.

\\\hat V(h)\\ converges at every point of the finite grid, so when \\V\\
has a unique minimiser \\h^\*\\ on the grid the selected bandwidth
equals \\h^\*\\ with probability tending to one. The estimate at the
selected bandwidth is then asymptotically normal with variance
\\V(h^\*)\\, and the Wald interval has its nominal level asymptotically.
No assumption on the distribution of the covariates is involved beyond
those of the corrected equation.

The default grid is `n_h` log-spaced bandwidths up to the largest
within-subject lag \\\|T\_{ij} - S\_{ik}\|\\, starting at the smallest
\\h\\ at which at least \\\lceil n^{1/2} \rceil\\ subjects have a pair
inside the kernel (for `correct = "pairs"`, both a response-covariate
pair and a covariate pair). `K`, `repeats` and `seed` are not used.

## Synchronous-loss cross-validation (`method = "sync"`)

For each candidate \\h\\ and fold, the corrected equation is fitted on
the training subjects. Every within-subject pair \\(T, S)\\ of the
held-out subjects gets the unweighted loss of \\Y(T)\\ against the
fitted mean at \\X(S)\\: squared error (identity), Poisson deviance
(log) or Bernoulli log loss (logistic). The losses are regressed on an
intercept and the entries of the estimated variogram matrix
\\\hat\Gamma(\|T - S\|)\\, \\\Gamma(u) = E\[\\x(t) - x(t-u)\\\\x(t) -
x(t-u)\\^\top\]\\, of the covariates that vary within subjects. The
intercept, the held-out loss at zero asynchrony, is retained. Each entry
of \\\hat\Gamma\\ comes from isotonic regressions of squared increments
of \\x_j\\ and \\x_j + x_k\\ on the lag. For the identity link and a
stationary covariate the intercept's expectation is \\\sigma^2 +
\delta_h^\top \Sigma_X \delta_h\\, with \\\delta_h\\ the error of the
training fit, because the asynchrony part of the loss is exactly linear
in \\\Gamma(u)\\. Intercepts are averaged over folds (and over `repeats`
splits) and the minimiser is selected.

A candidate is retained only when every fold produces a fit whose
intercept is at least the mean saturated loss of the held-out pairs
(\\\min\_\eta\\ of the loss: 0 for squared error and log loss, \\y - y
\log y\\ for the Poisson loss). An intercept below it comes from a
degenerate fold fit.

The default grid is `n_h` log-spaced bandwidths up to the largest
within-subject lag \\\|T\_{ij} - S\_{ik}\|\\, starting at the smallest
\\h\\ at which every training fold has at least \\\lceil
n\_{\mathrm{train}}^{1/2} \rceil\\ subjects with a pair inside the
kernel (for `correct = "pairs"`, both a response-covariate pair and a
covariate pair). Below that, a fold fit uses a handful of pairs and its
error does not shrink with \\n\\.

The guarantee is consistency: the selected \\h\\ lies on a finite grid
on which every corrected fit is consistent. The criterion's differences
between candidates are of the same order as its noise, so the selected
\\h\\ does not converge. Because the same subjects choose \\h\\ and fit
\\\beta\\, the estimate at the selected \\h\\ can carry a selection bias
of the order of its standard error, and the Wald interval has no
coverage guarantee.

## Plain cross-validation (`method = "plain"`)

Each candidate is fitted without correction and scored on the held-out
subjects by the kernel-weighted squared error of the fitted mean at that
bandwidth, \$\$ \frac{\sum\_{i \in \mathrm{test}} \sum_j \sum_k
W_h(T\_{ij} - S\_{ik}) \[Y_i(T\_{ij}) - g\\X_i(S\_{ik})^\top
\hat\beta\\\]^2} {\sum\_{i \in \mathrm{test}} \sum_j \sum_k
W_h(T\_{ij} - S\_{ik})}.\$\$ The default grid runs from the smallest to
the largest within-subject lag (positive lags \\T\_{ij} - S\_{ik}\\ only
when `one_sided = TRUE`; exact zeros excluded). A candidate is dropped
only when it fails on every fold. This criterion is minimised at a
bandwidth much smaller than the MSE-optimal one, and the uncorrected fit
keeps a bias of the order of its standard error there.

## How the folds are formed

Folds are formed over subjects, never over rows. Rows within a subject
come from one trajectory, so splitting by row would put a subject on
both sides of the split and leak the answer into the held-out score.

Always look at `cv_results`. If the minimum sits at an endpoint of the
grid a warning is issued.

## See also

[`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 150, beta = c(0.5, 1.5))
cv <- d$y |>
  kee_async_cv(d$x, y ~ x,
    id = id, time = time,
    n_h = 20, K = 3, seed = 1, quiet = TRUE
  )
#> Warning: the selected bandwidth is at an endpoint of 'h_grid'; widen the grid to check that the minimum is interior
cv
#> Call:
#> kee_async_cv(data_y = d$y, data_x = d$x, formula = y ~ x, id = id, 
#>     time = time, n_h = 20, K = 3, seed = 1, quiet = TRUE)
#> 
#> Minimum-variance bandwidth: tr{Sigma_X V(h)} over the grid (pairs-corrected fits on all subjects)
#> 
#> # A tibble: 20 × 2
#>          h   trSV
#>      <dbl>  <dbl>
#>  1 0.00374 1.19  
#>  2 0.00502 0.560 
#>  3 0.00673 0.316 
#>  4 0.00903 0.220 
#>  5 0.0121  0.179 
#>  6 0.0162  0.171 
#>  7 0.0218  0.171 
#>  8 0.0292  0.138 
#>  9 0.0392  0.109 
#> 10 0.0526  0.0881
#> 11 0.0706  0.0872
#> 12 0.0946  0.0834
#> 13 0.127   0.0952
#> 14 0.170   0.0953
#> 15 0.228   0.0943
#> 16 0.306   0.0925
#> 17 0.411   0.0893
#> 18 0.551   0.0794
#> 19 0.739   0.0730
#> 20 0.992   0.0696
#> 
#> Selected h = 0.991538
#> 
#> Coefficients at the selected h:
#> (Intercept)           x 
#>   0.6313399   1.7001935 
cv$cv_results
#> # A tibble: 20 × 2
#>          h   trSV
#>      <dbl>  <dbl>
#>  1 0.00374 1.19  
#>  2 0.00502 0.560 
#>  3 0.00673 0.316 
#>  4 0.00903 0.220 
#>  5 0.0121  0.179 
#>  6 0.0162  0.171 
#>  7 0.0218  0.171 
#>  8 0.0292  0.138 
#>  9 0.0392  0.109 
#> 10 0.0526  0.0881
#> 11 0.0706  0.0872
#> 12 0.0946  0.0834
#> 13 0.127   0.0952
#> 14 0.170   0.0953
#> 15 0.228   0.0943
#> 16 0.306   0.0925
#> 17 0.411   0.0893
#> 18 0.551   0.0794
#> 19 0.739   0.0730
#> 20 0.992   0.0696
```

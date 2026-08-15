# skmle

Longitudinal covariates are rarely measured when you need them. A
biomarker is drawn at a handful of clinic visits; the event you are
modelling happens somewhere between two of them, and the outcome you are
modelling is recorded on a schedule of its own. Carrying the last value
forward answers a different question, and smoothing the covariate first
and substituting it is regression calibration, which is not consistent
when the sampling is sparse.

`skmle` handles the mismatch directly, by weighting each observation
according to how far its measurement time sits from the time being
modelled. It covers two settings:

For survival outcomes, the transformed hazards family, with covariates
observed sparsely and intermittently. The Box-Cox parameter `s` indexes
the family: `s = 0` is proportional hazards, `s = 1` is additive
hazards, and other values interpolate. This is the Sieve Maximum
Kernel-weighted Log-likelihood Estimator (SMKLE) of Sun, Sun, Zhao and
Cao (2025), plus faster specialised estimating equations for the two
named cases.

For longitudinal outcomes, generalised linear models where the response
and the covariate are measured on *different* time grids and never
observed together. This is the asynchronous case. This is the
kernel-weighted estimating equations of Cao, Zeng and Fine (2015), with
either time-invariant coefficients or a coefficient curve `β(t)`.

Both settings share one idea and one implementation: a kernel weight
bridging mismatched time grids, with the heavy numerical work in C++ via
`Rcpp` and `RcppArmadillo`, behind an ordinary R modelling interface.

## Which function do I need?

Two questions.

First, what is the outcome? A *time to an event* (death, relapse,
failure), possibly censored? Use the survival estimators; your data is
one long table. Or a *repeatedly measured quantity* (a score, a lab
value) recorded on its own schedule? Use the asynchronous estimators;
your data is two tables.

Then:

| Outcome | Situation | Function |
|:---|:---|:---|
| Survival | Start here; Cox model | [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md) |
| Survival | Additive hazards instead | [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md) |
| Survival | Want the baseline hazard, or a model between the two | [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) |
| Survival | Choose the bandwidth properly | [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md) |
| Longitudinal | Start here; one constant effect | [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md) |
| Longitudinal | Choose the bandwidth properly | [`kee_async_cv()`](https://dayusun.github.io/skmle/reference/kee_async_cv.md) |
| Longitudinal | The effect may change over time | [`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md) |

Every fitting function picks a bandwidth for you if you do not supply
one, and says in a message what it chose. That is enough for a first
answer; the `_cv()` functions choose it from the data, which is what to
report.

## The functions

**Survival outcomes**

- [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) fits
  the general transformed hazards model, by sieve maximum
  kernel-weighted likelihood.
- [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md)
  fits proportional hazards, by kernel estimating equation.
- [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
  fits additive hazards, in closed form.
- [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
  selects the bandwidth by subject-level cross-validation, with a refit
  on the full data.
- [`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md)
  simulates sparse longitudinal survival data.

**Asynchronous longitudinal outcomes**

- [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  fits one constant coefficient per covariate.
- [`kee_async_cv()`](https://dayusun.github.io/skmle/reference/kee_async_cv.md)
  selects the bandwidth by subject-level cross-validation.
- [`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
  estimates a coefficient curve `β(t)`, pointwise.
- [`sim_async_data()`](https://dayusun.github.io/skmle/reference/sim_async_data.md)
  simulates a response and a covariate on independent observation-time
  streams.

Fitted objects support [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and the broom
generics [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`augment()`](https://generics.r-lib.org/reference/augment.html), all
returning tibbles. The two data frames come first in the asynchronous
estimators so they compose with the native pipe.

Every estimator takes `one_sided`, choosing between a half kernel (only
measurements *before* the modelled time contribute, the causal reading
and the risk-set restriction of a hazard model) and a full, two-sided
kernel.

## Installation

Install the development version from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("dayusun/skmle", build_vignettes = TRUE)
```

`build_vignettes = TRUE` matters if you want to read the articles from R
with [`vignette()`](https://rdrr.io/r/utils/vignette.html). Without it
the code installs fine but the vignettes are skipped. You can also read
them on the website, linked above, without installing anything.

Because the package compiles C++ code, you need a working toolchain such
as `Rtools` on Windows or the Xcode command line tools on macOS.

## Survival outcomes

``` r

library(skmle)
library(survival)

set.seed(123)

dat <- sim_skmle_data(
  n = 200,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,   # Box-Cox transformation parameter; s = 0 is proportional hazards
  cen = 0.7
)

fit <- skmle(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,
  h = 0.5,
  nknots = 3
)

summary(fit)
plot(fit)   # estimated baseline component
```

The returned objects follow standard R model conventions:
[`print()`](https://rdrr.io/r/base/print.html) shows the call and
coefficients, [`summary()`](https://rdrr.io/r/base/summary.html) reports
estimates with standard errors, z statistics and p-values, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
estimated baseline.

If the scientific model is proportional or additive hazards
specifically, the specialised estimating equations are considerably
faster than the full sieve fit:

``` r

fit_cox <- kee_cox(
  Surv(X, delta) ~ covariates,
  data = dat, id = id, obs_times = obs_times, h = 0.5
)
```

Bandwidth selection:

``` r

cv_fit <- skmle_cv(
  Surv(X, delta) ~ covariates,
  data = dat, id = id, obs_times = obs_times,
  s = 0, K = 3, h_grid = c(0.3, 0.4, 0.5), quiet = TRUE
)

cv_fit$h_cv
summary(cv_fit$fit)
```

## Asynchronous longitudinal outcomes

Here the outcome is itself a sparsely observed process, and its
measurement times do not line up with the covariate’s. The two are
supplied as separate tables, because there is no single data frame
holding both without inventing rows:

``` r

set.seed(202)
d <- sim_async_data(n = 300, beta = c(0.5, 1.5))

head(d$y)   # id, time, y   -- response occasions
head(d$x)   # id, time, x   -- covariate occasions, on a different grid

# Data comes first, so it pipes. The formula spans both tables:
# `y` is looked up in the first, `x` in the second.
fit_a <- d$y |>
  kee_async(d$x, y ~ x, id = id, time = time, h = 0.25)

summary(fit_a)
tidy(fit_a, conf.int = TRUE)   # a tibble, ready for ggplot2 or dplyr
```

`h` is the one consequential choice;
[`kee_async_cv()`](https://dayusun.github.io/skmle/reference/kee_async_cv.md)
makes it by cross-validation over subjects:

``` r

cv <- d$y |>
  kee_async_cv(d$x, y ~ x, id = id, time = time, K = 5, seed = 1)

cv$h_cv
```

Each (response, covariate) pair within a subject contributes in
proportion to its time separation, so nothing is discarded and nothing
is carried forward. The estimator converges at the smoothing rate
`(nh)^(1/2)`, so standard errors shrink more slowly than in a parametric
fit.

If the coefficients vary with time,
[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
estimates the curve pointwise. Both time arguments are smoothed there,
so it converges at the bivariate rate `(n h1 h2)^(1/2)` and the bands
are correspondingly wide:

``` r

fit_td <- d$y |>
  kee_async_td(d$x, y ~ x,
               id = id, time = time,
               times = seq(0.2, 0.8, by = 0.05), h = 0.25)

plot(fit_td)          # beta_j(t) with pointwise 95% bands
tidy(fit_td)          # one row per (time, term)
```

There is a full article on this setting: why last-value-carried-forward
and regression calibration fail, how to read the bandwidth diagnostics,
and what the half kernel changes.

``` r

vignette("asynchronous", package = "skmle")
```

or read it at <https://www.sundayu.me/skmle/articles/asynchronous.html>.

Identity, log and logistic links are supported in both.

## Benchmarks against SurvSparse

The package includes a benchmark vignette comparing `skmle` with
`SurvSparse` on matched sparse longitudinal survival-data settings.

In the current benchmark sweep:

| Scenario | `SurvSparse` median | `skmle` median | Relative result |
|:---|---:|---:|:---|
| Additive hazards, `n = 200` | ~288 ms | `kee_additive`: ~87 ms | `skmle` faster |
| Additive hazards, `n = 500` | ~835 ms | `kee_additive`: ~116 ms | `skmle` faster |
| Transformed hazards, `n = 100` | ~1096 ms | `skmle(s = 0)`: ~137 ms | `skmle` much faster |
| Transformed hazards, `n = 200` | ~2238 ms | `skmle(s = 0)`: ~216 ms | `skmle` much faster |

The additive comparison also includes the general spline-based
`skmle(s = 1)` fit, which remains competitive but is slower than
[`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
because it solves the broader joint optimization problem rather than a
specialized estimating equation.

For the full code and benchmark setup, see the vignette:

``` r

vignette("benchmark_survsparse", package = "skmle")
```

For a package tutorial covering both settings, see:

``` r

vignette("tutorial", package = "skmle")
```

All three are also on the website:
<https://www.sundayu.me/skmle/articles/>.

## References

Sun, D., Sun, Z., Zhao, X., & Cao, H. (2025). Kernel Meets Sieve:
Transformed Hazards Models with Sparse Longitudinal Covariates. *Journal
of the American Statistical Association, 120*(552), 2580-2591.
<https://doi.org/10.1080/01621459.2025.2476781>

Cao, H., Zeng, D., & Fine, J. P. (2015). Regression Analysis of Sparse
Asynchronous Longitudinal Data. *Journal of the Royal Statistical
Society: Series B, 77*(4), 755-776. <https://doi.org/10.1111/rssb.12086>

# Fit a Cox-Type KEE Model

Fit the proportional hazards model for sparse longitudinal covariate
data using a kernel estimating-equation approach.

## Usage

``` r
kee_cox(
  formula,
  data,
  id,
  obs_times,
  h = NULL,
  weight = NULL,
  one_sided = TRUE
)
```

## Arguments

- formula:

  A model formula with a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  response.

- data:

  Data frame containing all variables used in the fit.

- id:

  Subject identifier aligned row-wise with `data`.

- obs_times:

  Longitudinal observation times aligned row-wise with `data`. Times may
  be on any scale; the sieve basis and the cumulative-hazard quadrature
  are built on the observed follow-up, so there is no need to rescale to
  the unit interval first. `h` must be on the same scale. Observation
  times may be negative; see “Negative observation times” in
  `kee_cox()`.

- h:

  Positive kernel bandwidth. If omitted, one is read off the observation
  times as a rule of thumb and reported in a message. Use
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md) to
  choose it from the data.

- weight:

  Weight function of the standardised lag \\(t - r)/h\\, or `NULL` (the
  default) for the Epanechnikov kernel implied by `one_sided`. Any R
  function will do; see
  [kernel-weights](https://www.sundayu.me/skmle/reference/kernel-weights.md).

- one_sided:

  Logical. `TRUE` (the default) uses a half kernel: only covariate
  observations strictly before the event or quadrature time inform that
  time, which is the risk-set restriction and the estimator as
  published. `FALSE` uses a full, two-sided kernel, smoothing the
  covariate path from both sides. The switch applies to the risk-set
  averages inside the C++ backend as well as to the row weights, so the
  two are always consistent.

## Value

An object of class `kee` containing coefficient estimates, the estimated
variance-covariance matrix, the estimating-equation matrices,
convergence status, and the original function call.

## Details

`kee_cox()` targets the proportional hazards case without estimating a
nonparametric baseline component. It is therefore a useful specialized
alternative to
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) when the
scientific model is Cox-type and the main interest is in the regression
coefficients.

## Negative observation times

`obs_times` may be negative. Event times may not: the sieve basis and
the cumulative-hazard quadrature are built on `[0, max(X)]`, so an event
time outside that span is extrapolation. Observation times carry no such
constraint. They enter the score only through the lag `X - obs_times`,
and the half kernel admits a row on the sign of that lag rather than on
the sign of the time, so a covariate observed at `-180` and an event at
`30` sit a lag of `210` apart and are handled like any other pair.

Reach still decides whether such a row counts for anything. A pre-origin
observation contributes only where it falls within `h` of an event time,
so one further back than a bandwidth is silently outside every window
and the fit is identical to the fit with that row deleted. Admitting
negative times widens what is representable, not what is reachable.

**The permission is numerical; the interpretation is not.** A negative
observation time means a covariate measured before the time origin, that
is, before the participant entered follow-up. Shifting the origin so
that every time is positive is a different operation: it moves `X` too,
leaves every lag and every risk set unchanged, and buys nothing. Holding
`X` anchored at zero and admitting negative observation times adds
information the fit would otherwise discard, and that is a claim about
the data:

- **The pre-entry measurement must be commensurable with the on-study
  ones.** A biomarker assayed on a different platform, in a different
  laboratory, or under a different protocol before enrolment is not the
  same variable, and the estimator cannot detect that it is not. Carry a
  provenance indicator and check that it does not predict the outcome.

- **Entry is often not independent of the covariate.** If a participant
  was enrolled because of the value measured beforehand, conditioning on
  entry has already selected on the covariate, and the pre-entry row
  reintroduces that selection into the risk-set average.

- **Immortal time.** A participant contributes a pre-entry covariate
  only by having survived to enrol. Nothing in the weighting corrects
  for this.

None of these is a defect in the estimator, and none of them is checked.
If the pre-entry observations are ordinary measurements on the same
scale, for example a run-in visit or a screening draw, using them is the
point. If they come from elsewhere, the fit will be quiet and wrong.

## References

Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets
Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
*Journal of the American Statistical Association* (2025): 1-12.

Cao, Hongyuan, et al. "Inference for Cox models with sparse longitudinal
covariates." *Biometrika* (2015).

## Examples

``` r
# \donttest{
library(survival)

set.seed(123)
dat <- sim_skmle_data(
  n = 80,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)

fit_cox <- kee_cox(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

summary(fit_cox)
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#> Cox-type proportional hazards, kernel estimating equation (half kernel)
#>   n= 80 subjects   bandwidth h = 0.5
#> 
#>             Estimate Std. Error z value Pr(>|z|)   
#> covariates1  1.02085    0.36037  2.8328 0.004615 **
#> covariates2 -0.36751    0.29950 -1.2271 0.219797   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```

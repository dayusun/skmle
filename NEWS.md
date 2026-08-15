# skmle 0.1.0

## Scope

The package covers two settings, not one. Alongside the transformed hazards
models for survival outcomes it now fits generalised linear models for
asynchronous longitudinal outcomes, where the response and the covariate are
recorded on different time grids. The `Title` drops its "for Survival Models"
restriction, and the `Description`, README, package help page and tutorial have
been rewritten to present the two settings on equal footing rather than
treating the second as an extra.

## Interface

* The asynchronous estimators take the two data frames **first**, so they
  compose with the native pipe:
  `data_y |> kee_async(data_x, y ~ x, id = id, time = time, h = 0.25)`.
  The formula still spans both tables -- its left-hand side is looked up in
  `data_y`, its right-hand side in `data_x` -- which is unavoidable when the
  data genuinely lives in two frames, and is now stated in the argument docs
  rather than left to be discovered.
* `id` and `time` accept a bare name, a string, or `{{ col }}`, so the
  estimators can be wrapped in other functions. They previously used
  `substitute()`, which deparsed `{{ col }}` literally.
* `tidy()`, `glance()` and `augment()` methods. `augment()` attaches `.fitted`
  to the covariate table, which is where a fitted value lives here; there is no
  `.resid`, because a residual would need a response value at the covariate
  time and that is exactly what asynchronous data lacks.
* `sim_async_data()`, `confint()`, and both `cv_results` tables return tibbles.
* Requires R (>= 4.1) for the native pipe used throughout the examples.

## Usability

* `kee_async()` and `kee_async_td()` take a formula spanning the two tables --
  `kee_async(y ~ x, data_y, data_x, id, time, h)` -- instead of seven
  positional vectors. Swapping the response and covariate tables used to return
  plausible numbers with no complaint.
* `kee_async_cv()` selects the bandwidth by subject-level cross-validation,
  scoring candidates by kernel-weighted squared error on held-out subjects. The
  asynchronous estimators previously offered no guidance on `h` at all.
* `kee_cox()` and `kee_additive()` now check that times lie on `[0, 1]`, which
  `skmle()` already did. `kee_additive()` genuinely requires it -- its
  quadrature is built on `[0, 1]` -- so times in other units were silently
  wrong rather than merely unusual.
* `kee_async()` warns when fewer than 5% of response occasions have a covariate
  observation in their window, which is the signature of a bandwidth on the
  wrong scale. The asynchronous estimators are scale-free, so this cannot be
  checked by a range test.
* `vcov()`, `nobs()` and (for `kee_td`) `confint()` methods. `confint()`
  previously failed on every fitted object in the package, because there was no
  `vcov()` for the default method to call.
* `skmle()` and `skmle_cv()` no longer take `norder`. It was validated and
  documented but never used: the sieve basis is a natural cubic spline, whose
  order is fixed. Existing calls that pass it will now error, and should drop
  the argument.
* `sim_async_data()` returns the covariates as plain columns (`x`, or `x1`,
  `x2`, ...) rather than a matrix column, so they can be named in a formula.
* New article, `vignette("asynchronous")`: why last-value-carried-forward and
  regression calibration are inconsistent here, how to read the bandwidth
  sensitivity plot, what the half kernel changes, and how to get the units
  right.

## Asynchronous longitudinal data

* `kee_async()` and `kee_async_td()` implement the kernel-weighted estimating
  equations of Cao, Zeng and Fine (2015) for a longitudinal response and a
  longitudinal covariate observed on **different** time grids, with
  time-invariant and time-dependent coefficients respectively. Identity, log
  and logistic links are supported.
* `sim_async_data()` simulates from their Section 4 design.
* `plot()` and `print()` methods for the `kee_td` coefficient curves.
* Both estimators are backed by C++: pairs are enumerated once and collapsed
  onto covariate rows, so the pair design never enters the Newton loop. The
  time-dependent weight factorises over the two occasion indices, so no pair is
  enumerated there at all.

## Half and full kernels

* `skmle()`, `kee_cox()`, `kee_additive()`, `skmle_cv()`, `kee_async()` and
  `kee_async_td()` take a `one_sided` argument. `TRUE` is the half kernel:
  only covariate observations preceding a time inform it. `FALSE` is the full,
  two-sided kernel. The survival estimators default to `TRUE`, the estimator as
  published; the asynchronous ones default to `FALSE`, the kernel of Cao, Zeng
  and Fine.
* For the survival estimators the switch reaches the risk-set averages in the
  C++ backend as well as the row weights, so the two halves of an estimator
  cannot disagree. `skmle_cv()` also uses it inside the fold loop, which
  previously had the half kernel hardcoded, so the bandwidth is now selected
  under the kernel the refit uses.
* The Epanechnikov kernel had been written out inline in three fitting
  functions; it now lives in one internal helper, so the half/full switch
  cannot be applied to two of the three.

## Initial release

* Initial CRAN release.
* `skmle()` fits transformed hazards models by sieve maximum kernel-weighted
  log-likelihood estimation (SMKLE).
* `kee_cox()` and `kee_additive()` fit Cox and additive hazards models via
  kernel-weighted estimating equations.
* `skmle_cv()` selects the kernel bandwidth by K-fold cross-validation.
* `sim_skmle_data()` simulates survival data with sparse, intermittently
  observed longitudinal covariates.
* `plot()`, `print()` and `summary()` methods for fitted objects.

# Changelog

## skmle 0.1.0

### Any weight function, not only the ones supplied

- **[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md),
  [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md),
  [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
  and [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
  take a `weight` argument.** Any R function of the standardised lag
  `u = (t - r)/h` may be passed; the fitting code privileges no
  particular kernel.
  [`w_epan_half()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)
  and
  [`w_epan_full()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)
  build the two Epanechnikov kernels and are exported.
- A weight takes one argument and returns the shape it was given. The
  sieve quadrature evaluates it on a matrix of node-by-observation lags,
  not only on a vector, so a weight built from
  [`apply()`](https://rdrr.io/r/base/apply.html) or a scalar `if()`
  would be recycled into the wrong cells. Both the length and the matrix
  shape are checked when the argument is resolved, rather than failing
  later inside the fit.
- A weight may take negative values and need not integrate to 1. The
  normalisation cancels in the risk-set ratios and is absorbed into `h`
  elsewhere, which is why
  [`w_epan_half()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)
  integrates to 1/2.
- `weight` defaults to `NULL`, meaning the Epanechnikov that `one_sided`
  implies. A literal
  [`w_epan_half()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)
  default would have been wrong, and was wrong for one commit: its
  support is `[0, 1]`, so `one_sided = FALSE` with that default zeroes
  every negative lag through the support and returns a one-sided fit
  under a two-sided name. The package’s own suite caught it and now pins
  it.
- `one_sided` stays a separate argument, because it is a property of the
  model and not of the kernel. It is the risk-set restriction of a
  hazard model and applies to the sign of the lag after the weight is
  evaluated, so a two-sided weight under `one_sided = TRUE` is truncated
  at zero rather than rejected.
- [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md),
  [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  and
  [`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
  do not take `weight`; they use the full Epanechnikov of Cao, Zeng and
  Fine.

### The cross-validation loop runs in C++

- **The fold loop and the held-out score moved back into C++**, as
  `skmle_cv_cpp()`. The reason is not speed, which is dominated by the
  fit either way. Scoring in R required a second implementation of the
  model’s transformation: `R/utils.R` carried a `trans_link()` that
  mirrored `trans_fun()` in `src/skmle_cpp.cpp` by hand, floor and all,
  and because `trans_fun` is not exported nothing could compare them and
  no test did. The criterion and the objective now call the same
  function, which is a guarantee rather than a test someone has to
  remember to write. `trans_link()` is gone.
- The optimiser set-up is factored into one internal `fit_core()`,
  shared by `skmle_cpp_fit()` and the loop, so the fit the
  cross-validation scores cannot drift from the fit
  [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) returns.
- `locf_score()` stays in R. It is set-up rather than a loop, it needs
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html) and
  [`findInterval()`](https://rdrr.io/r/base/findInterval.html), and it
  evaluates no part of the model: it returns a basis, quadrature weights
  and row indices.
- The criterion is unchanged, and that is checked rather than asserted.
  The losses from the R loop were recorded before the port, and the C++
  loop reproduces them to 4.4e-16 with identical coefficients.
  `test-cv-reference.R` keeps the comparison in the suite.
- The fold loop evaluates the caller’s weight function itself, as an
  `Rcpp::Function`, so the kernel keeps one description. An earlier
  draft of this change gave C++ a polynomial description of the weight
  alongside the R closure and a test comparing the two; handing the
  closure across removes the second description rather than testing it.

### Covariate observation times may be negative

- **`obs_times` is no longer required to be non-negative.** Event times
  still are, because the sieve basis and the cumulative-hazard
  quadrature are built on `[0, max(X)]` and anything outside that span
  is extrapolation. An observation time is never an argument to a basis:
  it reaches the score only through the lag `X - obs_times`, and the
  one-sided rule tests the sign of that lag, not the sign of the time.
  The old check applied a constraint belonging to one argument to an
  argument that does not carry it.
- What this buys is a covariate measured before the participant entered
  follow-up, a screening draw or a run-in visit, without re-origining
  the time axis. Shifting the origin moves `X` too and leaves every lag
  and every risk set unchanged, so it was never a way to recover that
  information.
- **The permission is numerical and the interpretation is not.** A
  pre-entry measurement has to be commensurable with the on-study ones,
  entry may not be independent of the covariate that preceded it, and a
  participant contributes a pre-entry row only by having survived to
  enrol. None of this is checked.
  [`?kee_cox`](https://www.sundayu.me/skmle/reference/kee_cox.md) states
  it under “Negative observation times”.
- Reach is unchanged: a row further back than one bandwidth from every
  event time still contributes nothing, and the fit equals the fit with
  that row deleted.

### Weights that can take negative values

- **The estimating equations and the sieve likelihood guard on a nonzero
  row weight rather than a positive one.** A weight function supplied by
  a caller is not required to be nonnegative, and the old `> 0` test
  dropped negative-weight rows silently instead of failing, which made
  the score wrong with no diagnostic. The objective, its gradient, and
  the bread and meat of the sandwich now sum over the same rows;
  previously a negative weight would have put them on different row
  sets.
- The internal risk-set weights keep their `> 0` test, which is a
  skip-zero optimisation on a kernel that is nonnegative by
  construction, not a sign assumption about the caller’s weight.
- The kernels shipped with the package are nonnegative, so every result
  the package has produced is unchanged to the last bit; the difference
  appears only under a caller-supplied weight.

### Scope

The package covers two settings, not one. Alongside the transformed
hazards models for survival outcomes it now fits generalised linear
models for asynchronous longitudinal outcomes, where the response and
the covariate are recorded on different time grids. The `Title` drops
its “for Survival Models” restriction, and the `Description`, README,
package help page and tutorial have been rewritten to present the two
settings on equal footing rather than treating the second as an extra.

### Getting a first answer

The package assumed you already knew the method. Three things a newcomer
could not supply now have defaults, each announced rather than hidden.

- **`h` is optional in every fitting function.** When omitted, a
  rule-of-thumb bandwidth is read off the observation times – the
  geometric midpoint of the grid the matching `_cv()` function searches
  – and a message reports the value, says it is a starting point rather
  than a tuned choice, and names the function that tunes it. Supplying
  `h` keeps the message quiet.
- **`s` defaults to `0` in
  [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) and
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)**,
  the proportional hazards model, so the Box-Cox parameter no longer has
  to be understood before the first fit.
- **`times` defaults in
  [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)**
  to 25 points spanning the 10th to 90th percentile of the observed
  response times, which keeps them away from the edges where a one-sided
  window makes the curve unreliable.

Together these mean
`kee_cox(Surv(X, delta) ~ z, data, id = id, obs_times = obs_times)` and
`data_y |> kee_async(data_x, y ~ x, id = id, time = time)` are complete
calls.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) now name the model
  in words and report the bandwidth and kernel, so the output says what
  was fitted rather than assuming you remember.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `cv.skmle` or `cv.kee_async` object draws the held-out loss against
  bandwidth with the selection marked, which is the quickest way to see
  a minimum sitting on the edge of the grid.
- A “Which function do I need?” table in `?skmle-package`, the README
  and the tutorial.

### Interface

- The asynchronous estimators take the two data frames **first**, so
  they compose with the native pipe:
  `data_y |> kee_async(data_x, y ~ x, id = id, time = time, h = 0.25)`.
  The formula still spans both tables – its left-hand side is looked up
  in `data_y`, its right-hand side in `data_x` – which is unavoidable
  when the data genuinely lives in two frames, and is now stated in the
  argument docs rather than left to be discovered.
- `id` and `time` accept a bare name, a string, or `{{ col }}`, so the
  estimators can be wrapped in other functions. They previously used
  [`substitute()`](https://rdrr.io/r/base/substitute.html), which
  deparsed `{{ col }}` literally.
- [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html)
  methods.
  [`augment()`](https://generics.r-lib.org/reference/augment.html)
  attaches `.fitted` to the covariate table, which is where a fitted
  value lives here; there is no `.resid`, because a residual would need
  a response value at the covariate time and that is exactly what
  asynchronous data lacks.
- [`sim_async_data()`](https://www.sundayu.me/skmle/reference/sim_async_data.md),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and both
  `cv_results` tables return tibbles.
- Requires R (\>= 4.1) for the native pipe used throughout the examples.

### Usability

- [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
  and
  [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  take a formula spanning the two tables –
  `kee_async(y ~ x, data_y, data_x, id, time, h)` – instead of seven
  positional vectors. Swapping the response and covariate tables used to
  return plausible numbers with no complaint.
- [`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
  selects the bandwidth by subject-level cross-validation, scoring
  candidates by kernel-weighted squared error on held-out subjects. The
  asynchronous estimators previously offered no guidance on `h` at all.
- [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) and
  [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
  now check that times lie on `[0, 1]`, which
  [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) already
  did.
  [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
  genuinely requires it – its quadrature is built on `[0, 1]` – so times
  in other units were silently wrong rather than merely unusual.
- [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
  warns when fewer than 5% of response occasions have a covariate
  observation in their window, which is the signature of a bandwidth on
  the wrong scale. The asynchronous estimators are scale-free, so this
  cannot be checked by a range test.
- [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and (for `kee_td`)
  [`confint()`](https://rdrr.io/r/stats/confint.html) methods.
  [`confint()`](https://rdrr.io/r/stats/confint.html) previously failed
  on every fitted object in the package, because there was no
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) for the default method
  to call.
- [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) and
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md) no
  longer take `norder`. It was validated and documented but never used:
  the sieve basis is a natural cubic spline, whose order is fixed.
  Existing calls that pass it will now error, and should drop the
  argument.
- [`sim_async_data()`](https://www.sundayu.me/skmle/reference/sim_async_data.md)
  returns the covariates as plain columns (`x`, or `x1`, `x2`, …) rather
  than a matrix column, so they can be named in a formula.
- New article,
  [`vignette("asynchronous")`](https://www.sundayu.me/skmle/articles/asynchronous.md):
  why last-value-carried-forward and regression calibration are
  inconsistent here, how to read the bandwidth sensitivity plot, what
  the half kernel changes, and how to get the units right.

### Asynchronous longitudinal data

- [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
  and
  [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  implement the kernel-weighted estimating equations of Cao, Zeng and
  Fine (2015) for a longitudinal response and a longitudinal covariate
  observed on **different** time grids, with time-invariant and
  time-dependent coefficients respectively. Identity, log and logistic
  links are supported.
- [`sim_async_data()`](https://www.sundayu.me/skmle/reference/sim_async_data.md)
  simulates from their Section 4 design.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods for the
  `kee_td` coefficient curves.
- Both estimators are backed by C++: pairs are enumerated once and
  collapsed onto covariate rows, so the pair design never enters the
  Newton loop. The time-dependent weight factorises over the two
  occasion indices, so no pair is enumerated there at all.

### Bandwidth selection in `skmle_cv()`

- **The held-out loss no longer falls away with the bandwidth.** It was
  the kernel-weighted log-likelihood evaluated on the held-out fold, and
  the kernel weight is `W(u/h)/h`, so the weight a subject contributes
  shrinks as `h` grows and the criterion decreased monotonically
  whatever the fit was worth. The largest candidate in `h_grid` won
  every grid on every data set – including the automatic grid – and
  coefficients further from the truth were preferred to coefficients
  nearer it. Reported by a reader of the tutorial, whose CV curve was
  monotone across `n = 200/1000`, censoring `0.2/0.7` and five seeds.
- The score is now an ordinary log-likelihood evaluated on the held-out
  subjects, with the covariate path carried forward from the last
  observation. No kernel and no bandwidth enter it, so it is the same
  yardstick for every candidate and the comparison is about the fit. See
  [`?skmle_cv`](https://www.sundayu.me/skmle/reference/skmle_cv.md).
- **A minimum at an endpoint of `h_grid` now warns**, as
  [`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
  already did: the selected value is then the best of the values offered
  rather than a minimum. On the automatic grid this fires often, because
  the grid stops at the rate-based `tau * n^-0.3` while the
  finite-sample minimum of the loss commonly lies above it.
- **`cv_results` gains an `se` column**, the standard error of each loss
  across the folds, because the criterion is flat: at the sample sizes
  in the tutorial every candidate is within one standard error of the
  minimum. It also scores prediction of the held-out hazard, where the
  baseline can absorb attenuation in `beta`, so it leans towards more
  smoothing than the coefficients on their own would want. Over 10
  replicates at `n = 200` on a grid spanning 0.05 to 0.9, the squared
  error of `beta-hat` at the selected bandwidth averaged 0.155, against
  0.205 at the largest candidate – what the old criterion always
  returned – and 0.065 at the oracle bandwidth. A one-standard-error
  rule does *not* help: it lands in the noisy small-`h` end and scored
  0.205, and
  [`?skmle_cv`](https://www.sundayu.me/skmle/reference/skmle_cv.md) says
  so.
- **`h_grid = NULL` worked again.** The automatic grid referred to the
  end of follow-up before it was computed, so the documented default
  errored with `object 'tau' not found` for any data set.
- The fold loop moved from C++ to R. It only subsets matrices and calls
  the optimiser, which R does in a dozen lines, and the held-out score
  needs a spline basis at subject-specific quadrature nodes.

### Half and full kernels

- [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md),
  [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md),
  [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md),
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md),
  [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
  and
  [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  take a `one_sided` argument. `TRUE` is the half kernel: only covariate
  observations preceding a time inform it. `FALSE` is the full,
  two-sided kernel. The survival estimators default to `TRUE`, the
  estimator as published; the asynchronous ones default to `FALSE`, the
  kernel of Cao, Zeng and Fine.
- For the survival estimators the switch reaches the risk-set averages
  in the C++ backend as well as the row weights, so the two halves of an
  estimator cannot disagree.
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
  also uses it inside the fold loop, which previously had the half
  kernel hardcoded, so the bandwidth is now selected under the kernel
  the refit uses.
- The Epanechnikov kernel had been written out inline in three fitting
  functions; it now lives in one internal helper, so the half/full
  switch cannot be applied to two of the three.

### Initial release

- Initial CRAN release.
- [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) fits
  transformed hazards models by sieve maximum kernel-weighted
  log-likelihood estimation (SMKLE).
- [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) and
  [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
  fit Cox and additive hazards models via kernel-weighted estimating
  equations.
- [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
  selects the kernel bandwidth by K-fold cross-validation.
- [`sim_skmle_data()`](https://www.sundayu.me/skmle/reference/sim_skmle_data.md)
  simulates survival data with sparse, intermittently observed
  longitudinal covariates.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for fitted
  objects.

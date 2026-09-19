## Negative covariate observation times.
##
## `obs_times` may be negative; `X` may not.  The asymmetry is not an oversight.
## The sieve basis and the cumulative-hazard quadrature are built on
## `[0, max(X)]`, so an event time outside that span is extrapolation, while an
## observation time is never an argument to a basis at all: it reaches the score
## only through the lag `X - obs_times`, and the one-sided rule tests the sign
## of that lag rather than the sign of the time.
##
## The sharp test is invariance under a common shift.  Subtracting a constant
## from both `X` and `obs_times` leaves every lag and every risk-set ordering
## untouched, so the fit must not move even once the shift has driven a large
## share of the observation times below zero.  Times are in days here, both to
## give the shift room to bite and because the unit interval is no longer a
## requirement.

make_day_sim <- function(n = 120, seed = 7) {
  set.seed(seed)
  x_i <- stats::runif(n, 400, 2000)
  d_i <- stats::rbinom(n, 1, 0.6)
  rows <- lapply(seq_len(n), function(i) {
    r <- sort(stats::runif(sample(2:5, 1), 0, x_i[i]))
    data.frame(
      id = i, obs_times = r, X = x_i[i], delta = d_i[i],
      z1 = stats::rnorm(length(r)), z2 = stats::rnorm(length(r))
    )
  })
  do.call(rbind, rows)
}

fit_day <- function(dat, h = 365, ...) {
  kee_cox(survival::Surv(X, delta) ~ z1 + z2,
    data = dat, id = id, obs_times = obs_times, h = h, ...
  )
}

test_that("check_time_scale accepts negative observation times", {
  expect_silent(check_time_scale(c(1, 2, 3), c(-5, 0, 2)))
  expect_silent(check_time_scale(c(1, 2, 3), c(-5, -4, -3)))
})

test_that("check_time_scale still rejects negative event times", {
  expect_error(
    check_time_scale(c(-1, 2, 3), c(0, 1, 2)),
    "event times must be non-negative"
  )
})

test_that("check_time_scale still rejects non-finite times", {
  expect_error(check_time_scale(c(1, NA, 3), c(0, 1, 2)), "finite")
  expect_error(check_time_scale(c(1, 2, 3), c(0, Inf, 2)), "finite")
})

test_that("kee_cox is invariant to a common shift that makes obs_times negative", {
  skip_on_cran()
  dat <- make_day_sim()
  by <- 350

  shifted <- dat
  shifted$X <- shifted$X - by
  shifted$obs_times <- shifted$obs_times - by

  # The shift must actually exercise the relaxation and must leave the event
  # times legal, or the test would pass for the wrong reason.
  expect_gt(mean(shifted$obs_times < 0), 0.10)
  expect_gte(min(shifted$X), 0)

  fit0 <- fit_day(dat)
  fit1 <- fit_day(shifted)

  expect_equal(unname(fit1$coefficients), unname(fit0$coefficients),
    tolerance = 1e-10
  )
  expect_equal(unname(diag(fit1$var)), unname(diag(fit0$var)),
    tolerance = 1e-10
  )
})

test_that("the two-sided kernel is also invariant to the shift", {
  skip_on_cran()
  dat <- make_day_sim()
  shifted <- dat
  shifted$X <- shifted$X - 350
  shifted$obs_times <- shifted$obs_times - 350

  fit0 <- fit_day(dat, one_sided = FALSE)
  fit1 <- fit_day(shifted, one_sided = FALSE)

  expect_equal(unname(fit1$coefficients), unname(fit0$coefficients),
    tolerance = 1e-10
  )
})

test_that("a pre-origin covariate row within one bandwidth is used", {
  skip_on_cran()
  dat <- make_day_sim()

  # Move one subject's earliest row before the time origin.  The comparison is
  # against DELETING that row, which is what a non-negativity rule forces a
  # caller to do.  If the two agree, the row is being ignored.  `h = 730`
  # because the earliest event here is at day 400: at `h = 365` a row at day
  # -120 is outside every window and correctly contributes nothing, which is
  # the point of the companion test below.
  target <- which(dat$id == 1L)[1]
  moved <- dat
  moved$obs_times[target] <- -120
  dropped <- dat[-target, , drop = FALSE]

  fit_moved <- fit_day(moved, h = 730)
  fit_dropped <- fit_day(dropped, h = 730)

  expect_true(all(is.finite(fit_moved$coefficients)))
  expect_false(isTRUE(all.equal(
    unname(fit_moved$coefficients), unname(fit_dropped$coefficients),
    tolerance = 1e-8
  )))
})

test_that("a pre-origin row further back than the bandwidth contributes nothing", {
  skip_on_cran()
  dat <- make_day_sim()

  # Reach, not sign, is what decides whether a row counts.  The earliest event
  # is at day 400, so with `h = 365` no event time can look back as far as day
  # -120 and the fit must equal the fit with that row deleted.  Admitting
  # negative times does not smuggle unreachable observations into the score.
  target <- which(dat$id == 1L)[1]
  moved <- dat
  moved$obs_times[target] <- -120
  dropped <- dat[-target, , drop = FALSE]

  expect_equal(
    unname(fit_day(moved)$coefficients),
    unname(fit_day(dropped)$coefficients),
    tolerance = 1e-10
  )
})

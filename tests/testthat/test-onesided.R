## Half kernel vs full kernel for the survival estimators.
##
## `one_sided = TRUE` (the default) is the risk-set restriction: a covariate
## observation informs a time only if it PRECEDES it.  The restriction lives in
## two places -- the row weights built in R, and the risk-set averages built in
## C++ -- and the flag has to reach both, so the test perturbs covariate rows
## that no event time can ever look back at and checks the fit does not move.

make_onesided_sim <- function(n, s = 0, seed = 101) {
    set.seed(seed)
    sim_skmle_data(
        n = n,
        mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
        mu_bar = 8,
        alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
        beta = c(1, -0.5),
        s = s,
        cen = 0.7
    )
}

test_that("the half kernel ignores covariate rows no event time can look back at", {
    skip_on_cran()
    dat <- make_onesided_sim(80)

    # Rows whose covariate was observed at or after their own subject's exit
    # time.  Such a row is invisible to a half kernel in both places at once:
    # its own event-term weight is zero, and in any risk set it belongs to
    # (which requires X[i] <= X[k']) the lag X[i] - obs_times[k] is <= 0.
    late <- dat$obs_times >= dat$X
    skip_if_not(sum(late) >= 3, "no late covariate rows in this draw")

    bump <- dat
    bump$covariates[late, ] <- bump$covariates[late, ] + 5

    f0 <- kee_cox(survival::Surv(X, delta) ~ covariates,
        data = dat, id = id, obs_times = obs_times, h = 0.5
    )
    f1 <- kee_cox(survival::Surv(X, delta) ~ covariates,
        data = bump, id = id, obs_times = obs_times, h = 0.5
    )
    expect_equal(coef(f0), coef(f1), tolerance = 1e-10)

    # The full kernel does see them, so the same perturbation must move it.
    g0 <- kee_cox(survival::Surv(X, delta) ~ covariates,
        data = dat, id = id, obs_times = obs_times, h = 0.5, one_sided = FALSE
    )
    g1 <- kee_cox(survival::Surv(X, delta) ~ covariates,
        data = bump, id = id, obs_times = obs_times, h = 0.5, one_sided = FALSE
    )
    expect_false(isTRUE(all.equal(coef(g0), coef(g1))))
})

test_that("every survival estimator accepts a full kernel and stays finite", {
    skip_on_cran()
    dat <- make_onesided_sim(80)

    fc <- kee_cox(survival::Surv(X, delta) ~ covariates,
        data = dat, id = id, obs_times = obs_times, h = 0.5, one_sided = FALSE
    )
    expect_true(all(is.finite(coef(fc))))
    expect_true(all(diag(fc$var) >= 0))

    fa <- kee_additive(survival::Surv(X, delta) ~ covariates,
        data = make_onesided_sim(80, s = 1), id = id, obs_times = obs_times,
        h = 0.5, one_sided = FALSE
    )
    expect_true(all(is.finite(coef(fa))))

    fs <- skmle(survival::Surv(X, delta) ~ covariates,
        data = dat, id = id, obs_times = obs_times, s = 0, h = 0.5,
        nknots = 3, one_sided = FALSE
    )
    expect_true(all(is.finite(coef(fs))))
    expect_true(all(diag(fs$var) >= 0))
})

test_that("cross-validation tunes under the kernel the refit will use", {
    skip_on_cran()
    dat <- make_onesided_sim(60)
    cv <- function(...) {
        skmle_cv(survival::Surv(X, delta) ~ covariates,
            data = dat, id = id, obs_times = obs_times, s = 0,
            K = 3, h_grid = c(0.4, 0.6), seed = 1, quiet = TRUE, ...
        )
    }
    half <- cv()
    full <- cv(one_sided = FALSE)

    # Different kernel, different held-out loss: the fold loop rebuilds the
    # weights per bandwidth in C++ and had the half kernel hardcoded there.
    expect_false(isTRUE(all.equal(half$cv_results$cvloss, full$cv_results$cvloss)))
    # And the choice reaches the refit rather than being dropped from the call.
    expect_identical(full$fit$call$one_sided, FALSE)
})

test_that("the shared kernel helper matches the definition it replaced", {
    lag <- seq(-1, 1, by = 0.1)
    h <- 0.4
    kerfun <- function(xx) pmax((1 - xx^2) * 0.75, 0)
    expect_equal(
        kernel_weights(lag, h, one_sided = TRUE),
        kerfun(lag / h) / h * as.numeric(lag > 0)
    )
    expect_equal(kernel_weights(lag, h, one_sided = FALSE), kerfun(lag / h) / h)

    # Matrix shape survives: the sieve quadrature passes a node-by-row matrix.
    m <- matrix(lag[1:20], 4, 5)
    expect_equal(dim(kernel_weights(m, h)), c(4L, 5L))
})

test_that("the fit does not depend on the units time is measured in", {
  skip_on_cran()
  dat <- make_onesided_sim(120)
  cc <- 365
  days <- dat
  days$X <- dat$X * cc
  days$obs_times <- dat$obs_times * cc

  # The sieve basis and the cumulative-hazard quadrature are built on
  # [0, max(X)], so rescaling time and the bandwidth together must leave the
  # coefficients alone.  When both were pinned to [0, 1] this silently
  # integrated over the wrong interval for any other scale.
  fit <- function(d, h) {
    coef(skmle(survival::Surv(X, delta) ~ covariates,
      data = d, id = id, obs_times = obs_times, s = 0, h = h, nknots = 3
    ))
  }
  expect_equal(fit(dat, 0.3), fit(days, 0.3 * cc), tolerance = 1e-5)

  # kee_cox sees time only through comparisons and lags, so it is invariant
  # for the same reason and needs no rescaling at all.
  cox <- function(d, h) {
    coef(kee_cox(survival::Surv(X, delta) ~ covariates,
      data = d, id = id, obs_times = obs_times, h = h
    ))
  }
  expect_equal(cox(dat, 0.3), cox(days, 0.3 * cc), tolerance = 1e-4)

  # An additive hazard is a rate per unit time, so its coefficients carry
  # units of 1/time and scale by cc rather than staying put.
  add_dat <- make_onesided_sim(120, s = 1, seed = 202)
  add_days <- add_dat
  add_days$X <- add_dat$X * cc
  add_days$obs_times <- add_dat$obs_times * cc
  add <- function(d, h) {
    coef(kee_additive(survival::Surv(X, delta) ~ covariates,
      data = d, id = id, obs_times = obs_times, h = h
    ))
  }
  expect_equal(add(add_dat, 0.3), add(add_days, 0.3 * cc) * cc, tolerance = 1e-4)
})

test_that("times are still required to be finite and non-negative", {
  dat <- make_onesided_sim(40)
  bad <- dat
  bad$X[1] <- -1
  expect_error(
    kee_cox(survival::Surv(X, delta) ~ covariates,
      data = bad, id = id, obs_times = obs_times, h = 0.5
    ),
    "non-negative"
  )
})

test_that("skmle no longer advertises an argument it ignores", {
  expect_false("norder" %in% names(formals(skmle)))
  expect_false("norder" %in% names(formals(skmle_cv)))
})

test_that("vcov and confint work on the survival fits", {
  skip_on_cran()
  dat <- make_onesided_sim(60)
  f <- kee_cox(survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times, h = 0.5
  )
  expect_equal(vcov(f), f$var)
  ci <- confint(f)
  expect_equal(dim(ci), c(2L, 2L))
  expect_identical(nobs(f), 60L)

  s <- skmle(survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times, s = 0, h = 0.5, nknots = 3
  )
  expect_equal(vcov(s), s$var)
  expect_equal(dim(confint(s)), c(2L, 2L))
  expect_identical(nobs(s), 60L)
})

## ----------------------------------------------------- defaults for novices

test_that("a first fit needs no bandwidth and says what it picked", {
  skip_on_cran()
  dat <- make_onesided_sim(60)

  expect_message(
    f <- kee_cox(survival::Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times
    ),
    "rule of thumb"
  )
  expect_true(all(is.finite(coef(f))))
  expect_true(f$h > 0)

  # The chosen value must sit inside the grid skmle_cv() would search, since
  # that is the claim the message makes.
  n <- length(unique(dat$id))
  pos <- (dat$X - dat$obs_times)[dat$X - dat$obs_times > 0]
  expect_gte(f$h, max(min(pos), n^(-0.6)) - 1e-12)
  expect_lte(f$h, min(max(pos), n^(-0.3)) + 1e-12)

  # s defaults to the Cox model, so skmle() needs neither s nor h.
  expect_message(
    s0 <- skmle(survival::Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times, nknots = 3
    ),
    "rule of thumb"
  )
  expect_equal(s0$s, 0)
  expect_true(all(is.finite(coef(s0))))
})

test_that("supplying h keeps the message quiet", {
  skip_on_cran()
  dat <- make_onesided_sim(60)
  expect_no_message(
    kee_cox(survival::Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times, h = 0.5
    )
  )
})

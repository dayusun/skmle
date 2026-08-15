test_that("skmle fits s=0 (Cox) and returns finite estimates", {
  skip_on_cran()
  dat <- make_sim(n = 80, s = 0)
  fit <- skmle(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, h = 0.5, nknots = 3
  )
  expect_s3_class(fit, "skmle")
  expect_length(fit$coefficients, 2L)
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(is.finite(diag(fit$var))))
  expect_true(all(diag(fit$var) >= 0))
  expect_equal(fit$n, 80L)
  expect_equal(fit$s, 0)
  expect_equal(fit$h, 0.5)
})

test_that("skmle fits s!=0 (transformed hazards) and respects Box-Cox constraint", {
  skip_on_cran()
  dat <- make_sim(n = 60, s = 0.5)
  fit <- skmle(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0.5, h = 0.5, nknots = 3
  )
  expect_s3_class(fit, "skmle")
  expect_length(fit$coefficients, 2L)
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("skmle rejects intercept-only formula and bad inputs", {
  dat <- make_sim(n = 20)
  expect_error(
    skmle(survival::Surv(X, delta) ~ 1, data = dat, id = id,
          obs_times = obs_times, s = 0, h = 0.5),
    "at least one covariate"
  )
  expect_error(
    skmle(survival::Surv(X, delta) ~ covariates, data = dat, id = id,
          obs_times = obs_times, s = 0, h = -0.1),
    "positive"
  )
})

test_that("skmle rejects times outside [0, 1]", {
  dat <- make_sim(n = 20)
  dat$X <- dat$X * 10
  expect_error(
    skmle(survival::Surv(X, delta) ~ covariates, data = dat, id = id,
          obs_times = obs_times, s = 0, h = 0.5),
    "Event and observation times"
  )
})

test_that("skmle accepts non-numeric (character/factor) id", {
  skip_on_cran()
  dat <- make_sim(n = 30)
  dat$id_chr <- paste0("subj_", dat$id)
  fit <- skmle(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id_chr, obs_times = obs_times,
    s = 0, h = 0.5, nknots = 3
  )
  expect_equal(fit$n, 30L)
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("print/summary/plot methods work on skmle objects", {
  skip_on_cran()
  dat <- make_sim(n = 40)
  fit <- skmle(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, h = 0.5, nknots = 3
  )
  expect_output(print(fit), "Coefficients")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.skmle")
  expect_output(print(sm), "n=")
  expect_s3_class(plot(fit), "ggplot")
})

test_that("sim_skmle_data never returns a negative event time", {
  # The failure time solves Lambda(t) = -log(u).  Handing that to a general
  # root finder from a starting value let it converge to a spurious negative
  # root on roughly one draw in twenty at s = 1, which produced survival data
  # with X < 0.  Nothing rejected it until kee_additive() started checking the
  # time scale, and then it surfaced as a cross-platform test failure.
  for (s in c(0, 1)) {
    for (sd in 1:25) {
      set.seed(sd)
      d <- sim_skmle_data(
        n = 15,
        mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
        mu_bar = 8,
        alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
        beta = c(1, -0.5), s = s, cen = 0.7
      )
      expect_true(all(d$X > 0), label = paste0("s=", s, " seed=", sd, " min X"))
      expect_true(all(d$X <= 1))
      expect_true(all(d$obs_times >= 0 & d$obs_times <= 1))
      # Censored subjects must be exactly at their censoring time.
      expect_equal(d$X[!d$delta], d$censoring[!d$delta])
    }
  }
})

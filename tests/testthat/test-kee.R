test_that("kee_cox fits and returns finite estimates", {
  skip_on_cran()
  dat <- make_sim(n = 80, s = 0)
  fit <- kee_cox(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    h = 0.5
  )
  expect_s3_class(fit, "kee")
  expect_length(fit$coefficients, 2L)
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(diag(fit$var) >= 0))
  expect_equal(fit$n, 80L)
})

test_that("kee_additive fits and returns finite estimates", {
  skip_on_cran()
  dat <- make_sim(n = 80, s = 1)
  fit <- kee_additive(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    h = 0.5
  )
  expect_s3_class(fit, "kee")
  expect_length(fit$coefficients, 2L)
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(diag(fit$var) >= 0))
})

test_that("kee_* input validation triggers", {
  dat <- make_sim(n = 20)
  expect_error(kee_cox(survival::Surv(X, delta) ~ covariates, data = dat,
                       id = id, obs_times = obs_times, h = 0),
               "positive")
  expect_error(kee_additive(survival::Surv(X, delta) ~ covariates, data = dat,
                            id = id, obs_times = obs_times, h = -1),
               "positive")
})

test_that("print/summary on kee objects work", {
  skip_on_cran()
  dat <- make_sim(n = 40, s = 0)
  fit <- kee_cox(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times, h = 0.5
  )
  expect_output(print(fit), "Coefficients")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.kee")
  expect_output(print(sm), "n=")
})

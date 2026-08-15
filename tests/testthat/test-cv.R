test_that("skmle_cv selects an h on the supplied grid", {
  skip_on_cran()
  set.seed(2026)
  dat <- make_sim(n = 60, s = 0)
  h_grid <- c(0.3, 0.5, 0.7)
  cv <- skmle_cv(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, K = 3, h_grid = h_grid, nknots = 3, quiet = TRUE
  )
  expect_s3_class(cv, "cv.skmle")
  expect_true(cv$h_cv %in% h_grid)
  expect_equal(nrow(cv$cv_results), length(h_grid))
  expect_s3_class(cv$fit, "skmle")
  expect_equal(cv$fit$h, cv$h_cv)
})

test_that("skmle_cv clamps K when it exceeds n and warns", {
  skip_on_cran()
  set.seed(2026)
  dat <- make_sim(n = 15, s = 0)
  expect_message(
    skmle_cv(
      survival::Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times,
      s = 0, K = 50, h_grid = c(0.4, 0.6), nknots = 3, quiet = FALSE
    ),
    "exceeds"
  )
})

test_that("skmle_cv is reproducible under `seed`", {
  skip_on_cran()
  dat <- make_sim(n = 40, s = 0)
  cv1 <- skmle_cv(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, K = 3, h_grid = c(0.3, 0.5, 0.7), seed = 42, quiet = TRUE
  )
  cv2 <- skmle_cv(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, K = 3, h_grid = c(0.3, 0.5, 0.7), seed = 42, quiet = TRUE
  )
  expect_equal(cv1$h_cv, cv2$h_cv)
  expect_equal(cv1$fold_id, cv2$fold_id)
  expect_equal(cv1$cv_results$cvloss, cv2$cv_results$cvloss)
  expect_equal(cv1$seed, 42)
})

test_that("the cross-validation refit shows a call, not the whole data set", {
  skip_on_cran()
  set.seed(4)
  dat <- sim_skmle_data(
    n = 60,
    mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
    mu_bar = 8,
    alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
    beta = c(1, -0.5), s = 0, cen = 0.7
  )
  cv <- skmle_cv(survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, K = 3, h_grid = c(0.3, 0.5), seed = 1, quiet = TRUE
  )

  # skmle_cv() freezes `data` into the refit call so that lazy re-evaluation
  # cannot pick up a different frame than cross-validation used.  That is
  # correct, but the frozen frame must not survive into the printed call.
  txt <- paste(deparse(cv$fit$call), collapse = " ")
  expect_lt(nchar(txt), 300)
  expect_match(txt, "data = dat", fixed = TRUE)
  expect_false(grepl("structure(", txt, fixed = TRUE))

  # print() and summary() are the things that were unreadable.
  expect_lt(length(capture.output(print(cv$fit))), 20)
  expect_lt(length(capture.output(print(summary(cv$fit)))), 25)

  # The selected bandwidth still reaches the refit.
  expect_equal(cv$fit$call$h, cv$h_cv)
})

test_that("cv.skmle prints a summary rather than its internals", {
  skip_on_cran()
  set.seed(6)
  dat <- sim_skmle_data(
    n = 50,
    mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
    mu_bar = 8,
    alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
    beta = c(1, -0.5), s = 0, cen = 0.7
  )
  cv <- skmle_cv(survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, K = 2, h_grid = c(0.3, 0.5), seed = 1, quiet = TRUE
  )
  out <- capture.output(print(cv))
  expect_lt(length(out), 20)
  expect_true(any(grepl("cross-validation", out)))
  expect_true(any(grepl("Selected h", out)))
  # No raw list internals.
  expect_false(any(grepl("^\\$", out)))
})

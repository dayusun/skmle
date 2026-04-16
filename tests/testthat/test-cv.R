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

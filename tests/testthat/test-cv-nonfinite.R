# A bandwidth whose training fit fails or whose held-out score overflows on any
# fold must lose the grid as +Inf, never NaN, and a grid with no finite
# candidate must stop with a message rather than select its first entry.

cv_dat <- function() {
  set.seed(123)
  sim_skmle_data(
    n = 60, mu = function(tt) 8 * (0.75 + (0.5 - tt)^2), mu_bar = 8,
    alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
    beta = c(1, -0.5), s = 0, cen = 0.7
  )
}

test_that("a failed or overflowing candidate scores +Inf, not NaN", {
  skip_on_cran()
  dat <- cv_dat()
  r <- suppressWarnings(skmle_cv(survival::Surv(X, delta) ~ covariates, data = dat,
    id = id, obs_times = obs_times, s = 0, K = 3,
    h_grid = c(0.004, 0.01, 0.05, 0.2, 0.4), seed = 1, quiet = TRUE))
  expect_false(anyNA(r$cv_results$cvloss))
  expect_true(all(is.infinite(r$cv_results$cvloss[1:2])))
  expect_true(all(is.finite(r$cv_results$cvloss[3:5])))
  expect_true(is.finite(r$h_cv))
  expect_equal(r$h_cv, r$cv_results$h[which.min(r$cv_results$cvloss)])
})

test_that("a grid with no finite candidate stops instead of selecting", {
  skip_on_cran()
  dat <- cv_dat()
  expect_error(
    skmle_cv(survival::Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, K = 3,
             h_grid = c(0.002, 0.004, 0.006), seed = 1, quiet = TRUE),
    "no bandwidth in 'h_grid' gave a finite"
  )
})

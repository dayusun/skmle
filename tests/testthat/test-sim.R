test_that("sim_skmle_data returns a long-format data frame with expected columns", {
  dat <- make_sim(n = 30)
  expect_s3_class(dat, "data.frame")
  expect_true(all(c("id", "X", "delta", "covariates", "obs_times") %in% names(dat)))
  expect_true(is.matrix(dat$covariates))
  expect_equal(ncol(dat$covariates), 2L)
  expect_equal(length(unique(dat$id)), 30L)
  expect_true(all(dat$X > 0))
})

test_that("sim_skmle_data rejects invalid beta length", {
  expect_error(
    sim_skmle_data(
      n = 10,
      mu = function(tt) 8,
      mu_bar = 8,
      alpha = function(tt) 0.5,
      beta = c(1, -0.5, 0.2),
      s = 0,
      cen = 0.7
    ),
    "two covariates"
  )
})

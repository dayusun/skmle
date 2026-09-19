## The cross-validation reference.
##
## `ref-cv-losses.rds` holds the held-out losses produced by the R fold loop at
## commit c2286f4, the commit that replaced the bandwidth-dependent criterion.
## The loop is being moved back into C++; this file is what the port has to
## reproduce.  Without it the port would be checked against nothing, which is
## the situation that turns a rewrite into a silent change of answer.
##
## A difference here is not automatically a regression, but it is never a
## detail.  Either the port changed the arithmetic, or the reference is stale
## and was regenerated deliberately.  Say which in the commit message.

test_that("skmle_cv reproduces the recorded reference losses", {
  skip_on_cran()
  ref <- readRDS(test_path("ref-cv-losses.rds"))

  set.seed(4242)
  dat <- make_sim(n = 60, s = 0)
  got <- suppressWarnings(skmle_cv(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times,
    s = 0, K = 3, h_grid = c(0.15, 0.25, 0.35, 0.5),
    seed = 2026, quiet = TRUE
  ))

  expect_equal(got$cv_results$h, ref$cv_results$h)
  # Tolerance, not identity: the port may sum the held-out terms in a different
  # order.  1e-8 is far below the spacing between candidates, which is what the
  # criterion is for, and far above floating-point reassociation.
  expect_equal(got$cv_results$cvloss, ref$cv_results$cvloss, tolerance = 1e-8)
  expect_equal(got$cv_results$se, ref$cv_results$se, tolerance = 1e-8)
  expect_equal(got$h_cv, ref$h_cv)
  expect_equal(unname(got$fit$coefficients), ref$coef, tolerance = 1e-8)
})

test_that("the recorded reference is the one the port targets", {
  ref <- readRDS(test_path("ref-cv-losses.rds"))
  expect_identical(ref$generated, "R fold loop, commit c2286f4")
  expect_equal(nrow(ref$cv_results), 4L)
  expect_true(all(is.finite(ref$cv_results$cvloss)))
})

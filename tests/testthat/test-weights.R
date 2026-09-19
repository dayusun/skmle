## The weight argument.
##
## `weight` takes any R function of the standardised lag. The two Epanechnikov
## kernels shipped here are a starting point, not the only option, so these
## tests cover a weight that is not polynomial and one that takes negative
## values as well as the shipped pair.

fit_cv <- function(...) {
  set.seed(4242)
  dat <- make_sim(n = 60, s = 0)
  suppressWarnings(skmle_cv(
    survival::Surv(X, delta) ~ covariates,
    data = dat, id = id, obs_times = obs_times, s = 0, K = 3,
    h_grid = c(0.15, 0.25, 0.35, 0.5), seed = 2026, quiet = TRUE, ...
  ))
}

test_that("the shipped kernels are the documented functions", {
  half <- w_epan_half()
  full <- w_epan_full()

  u <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  expect_equal(half(u), c(0, 0, 0, 0.75, 0.5625, 0, 0))
  expect_equal(full(u), c(0, 0, 0.5625, 0.75, 0.5625, 0, 0))

  expect_equal(attr(half, "support"), c(0, 1))
  expect_equal(attr(full, "support"), c(-1, 1))

  # Shape preservation, which the sieve quadrature depends on.
  m <- matrix(u[1:6], nrow = 2)
  expect_identical(dim(half(m)), dim(m))
  expect_identical(dim(full(m)), dim(m))
})

test_that("weight = NULL is the kernel one_sided implies", {
  expect_equal(
    kernel_weights(c(-0.3, 0.2), 1, resolve_weight(NULL, TRUE), TRUE),
    kernel_weights(c(-0.3, 0.2), 1, w_epan_half(), TRUE)
  )
  expect_equal(
    kernel_weights(c(-0.3, 0.2), 1, resolve_weight(NULL, FALSE), FALSE),
    kernel_weights(c(-0.3, 0.2), 1, w_epan_full(), FALSE)
  )
})

test_that("naming the default kernel changes nothing", {
  skip_on_cran()
  # A literal default of w_epan_half() would have broken one_sided = FALSE,
  # because the support would still be [0, 1]. NULL is what preserves both.
  expect_equal(
    fit_cv(weight = w_epan_half())$cv_results$cvloss,
    fit_cv()$cv_results$cvloss
  )
})

test_that("a non-polynomial weight is accepted and used", {
  skip_on_cran()
  gauss_half <- function(u) exp(-u^2 / 2) * (u >= 0 & u <= 3)
  got <- fit_cv(weight = gauss_half)

  expect_true(all(is.finite(got$cv_results$cvloss)))
  expect_true(got$h_cv %in% c(0.15, 0.25, 0.35, 0.5))
  # It must actually be used, not quietly ignored in favour of the default.
  expect_false(isTRUE(all.equal(
    got$cv_results$cvloss, fit_cv()$cv_results$cvloss,
    tolerance = 1e-6
  )))
})

test_that("a signed weight is used rather than clamped away", {
  skip_on_cran()
  # W_2(u) = 4 - 6u is negative above u = 2/3. A `> 0` guard anywhere in the
  # score would drop those rows silently and the fit would differ from the
  # same weight with its negative part removed.
  signed_w <- function(u) (4 - 6 * u) * (u >= 0 & u <= 1)
  clamped <- function(u) pmax(4 - 6 * u, 0) * (u >= 0 & u <= 1)

  a <- fit_cv(weight = signed_w)$cv_results$cvloss
  b <- fit_cv(weight = clamped)$cv_results$cvloss
  expect_false(isTRUE(all.equal(a, b, tolerance = 1e-6)))
})

test_that("a weight that is not a function, or misshapen, is refused", {
  expect_error(resolve_weight(1:3, TRUE), "must be a function")
  expect_error(
    resolve_weight(function(u) "nope", TRUE),
    "one finite numeric value per element"
  )
  expect_error(
    resolve_weight(function(u) sum(u), TRUE),
    "one finite numeric value per element"
  )
  # Returns the right length but drops the matrix shape: this is the failure
  # that would otherwise recycle weights into the wrong quadrature cells.
  expect_error(
    resolve_weight(function(u) as.numeric(u) * 0 + 1, TRUE),
    "must preserve the shape"
  )
})

test_that("weight_support falls back rather than failing", {
  expect_equal(weight_support(w_epan_half()), c(0, 1))
  expect_equal(weight_support(function(u) u), c(-1, 1))
})

test_that("kernel_weights defaults to the kernel one_sided implies", {
  # Regression test. kernel_weights() briefly defaulted to a literal
  # w_epan_half(), which meant a caller passing only one_sided = FALSE got a
  # HALF kernel: the negative lags were zeroed by the weight's support before
  # the one-sided rule ever ran. test-onesided.R caught it. The fix is the
  # same NULL resolution the entry points use, and this pins it at the helper.
  lag <- seq(-1.2, 1.2, by = 0.2)
  h <- 0.8
  expect_equal(
    kernel_weights(lag, h, one_sided = FALSE),
    kernel_weights(lag, h, weight = w_epan_full(), one_sided = FALSE)
  )
  expect_true(any(kernel_weights(lag, h, one_sided = FALSE)[lag < 0] > 0))
  expect_true(all(kernel_weights(lag, h, one_sided = TRUE)[lag <= 0] == 0))
})

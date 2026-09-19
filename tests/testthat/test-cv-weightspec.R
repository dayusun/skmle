## The kernel is described twice, so the two descriptions are held together.
##
## `kernel_weights()` evaluates the Epanechnikov in R, for the fit. The
## cross-validation loop cannot call it, because it rebuilds the weights for
## every candidate bandwidth inside C++, so the same kernel also travels as the
## four fields `epan_weight_spec()` returns and `calc_kerfun()` reads back.
##
## Two descriptions of one kernel with nothing comparing them is the shape of
## defect this port was done to remove, so it is not allowed to reappear in the
## weight. The comparison is against `kernel_weights()`, the description the fit
## actually uses.

test_that("the weight spec reproduces kernel_weights on both kernels", {
  # calc_kerfun() in C++, transcribed: coefficients, support, mirror, and the
  # one-sided rule on the sign of the lag rather than of the time.
  from_spec <- function(lag, h, spec, one_sided) {
    z <- lag / h
    v <- if (spec$mirror) abs(z) else z
    acc <- rowSums(outer(v, seq_along(spec$coef) - 1L, "^") *
      rep(spec$coef, each = length(v)))
    acc[z < spec$a | z > spec$b] <- 0
    if (one_sided) acc[lag <= 0] <- 0
    acc / h
  }

  set.seed(99)
  lag <- c(runif(3000, -3, 3), 0, 1, -1)

  for (os in c(TRUE, FALSE)) {
    for (h in c(0.05, 0.3, 1, 7.5)) {
      spec <- epan_weight_spec(os)
      expect_equal(
        from_spec(lag, h, spec, os),
        kernel_weights(lag, h, os),
        tolerance = 1e-12,
        info = sprintf("one_sided = %s, h = %g", os, h)
      )
    }
  }
})

test_that("the spec is the plain Epanechnikov and nothing else", {
  # If this changes, the public build has acquired a weight family and the
  # release decision that kept those private needs revisiting.
  expect_equal(epan_weight_spec(TRUE)$coef, c(0.75, 0, -0.75))
  expect_equal(epan_weight_spec(TRUE)$a, 0)
  expect_equal(epan_weight_spec(FALSE)$a, -1)
  expect_equal(epan_weight_spec(TRUE)$b, 1)
  expect_false(epan_weight_spec(TRUE)$mirror)
})

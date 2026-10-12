## The interior-point optimiser.
##
## optimizer = "barrier" maximises the kernel-weighted log-likelihood with
## 1 + s * eta > 0 at every point the objective evaluates. At s = 0 there are
## no constraints and it must agree with SLSQP; at s != 0 its solution must be
## strictly feasible everywhere and no worse than SLSQP given the same
## constraint set.

barrier_fit <- function(s, n = 120, seed = 5, ...) {
  d <- make_sim(n = n, s = s, seed = seed)
  f <- survival::Surv(X, delta) ~ covariates
  list(
    d = d,
    sq = skmle(f, d, id = id, obs_times = obs_times, s = s, h = 0.3, ...),
    br = skmle(f, d, id = id, obs_times = obs_times, s = s, h = 0.3,
               optimizer = "barrier", ...)
  )
}

test_that("the barrier fit is the SLSQP fit when there is nothing to constrain", {
  skip_on_cran()
  r <- barrier_fit(s = 0)
  expect_equal(coef(r$br), coef(r$sq), tolerance = 1e-6)
  expect_equal(r$br$loglik, r$sq$loglik, tolerance = 1e-8)
  expect_identical(r$br$optimizer, "barrier")
})

test_that("the barrier fit is strictly feasible at every evaluated point", {
  skip_on_cran()
  r <- barrier_fit(s = 1)
  d <- r$d
  fit <- r$br
  expect_gte(fit$convergence, 0)
  tau <- max(d$X)
  basis <- function(t) {
    as.matrix(splines::ns(t, knots = tau * (1:3) / 4, intercept = TRUE,
                          Boundary.knots = c(0, tau)))
  }
  zb <- as.vector(as.matrix(d$covariates) %*% coef(fit))
  kv <- kernel_weights(d$X - d$obs_times, 0.3)
  ev <- d$delta == 1 & kv != 0
  expect_gt(min(1 + as.vector(basis(d$X) %*% fit$gamma)[ev] + zb[ev]), 0)
  lq <- gaussquad::legendre.quadrature.rules(64)[[64]]
  tts <- 0.5 * tau * (lq$x + 1)
  used <- kernel_weights(outer(tts, d$obs_times, "-"), 0.3) != 0 &
    outer(tts, d$X, "<")
  eta <- outer(as.vector(basis(tts) %*% fit$gamma), zb, "+")
  expect_gt(min((1 + eta)[used]), 0)
})

test_that("the barrier fit is no worse than SLSQP on the fully constrained problem", {
  skip_on_cran()
  # SLSQP given every evaluated point as a constraint solves the barrier's
  # problem; the concave objective has one optimum, so the two must meet.
  full <- skmle
  b <- deparse(body(full))
  at <- grep("res <- if \\(optimizer == \"barrier\"\\)", b)
  b <- append(b, c(
    "ineqmat <- unique(rbind(ineqmat, ineq_rescue))",
    "ineq_rescue <- ineqmat[0, , drop = FALSE]"
  ), after = at - 1)
  body(full) <- parse(text = paste(b, collapse = "\n"))[[1]]
  environment(full) <- asNamespace("skmle")
  r <- barrier_fit(s = 0.5)
  ref <- suppressWarnings(full(survival::Surv(X, delta) ~ covariates, r$d,
    id = id, obs_times = obs_times, s = 0.5, h = 0.3))
  expect_gte(r$br$loglik, ref$loglik - 1e-7)
  expect_equal(coef(r$br), coef(ref), tolerance = 1e-4)
})

test_that("optimizer is validated", {
  d <- make_sim(n = 30, s = 0, seed = 1)
  expect_error(
    skmle(survival::Surv(X, delta) ~ covariates, d, id = id,
          obs_times = obs_times, h = 0.3, optimizer = "newton"),
    "should be one of"
  )
})

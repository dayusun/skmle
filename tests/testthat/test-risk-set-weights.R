## The risk-set averages use the caller's weight.
##
## Every estimating equation here has a main term weighted by W(lag/h)/h and a
## risk-set average S1/S0 inside it. Both must use the same weight and the same
## one-sided rule (Sun et al. 2025, eq. 4 and Remark 7), or the fit solves an
## equation that belongs to no estimator. These tests rebuild each equation in
## plain R from the weight alone and check the C++ fit against it, with a
## weight that reaches past h, where a hardcoded Epanechnikov would differ.

gauss_half <- function(u) exp(-u^2 / 2) * (u >= 0 & u <= 3)
attr(gauss_half, "support") <- c(0, 3)

rs_setup <- function(s = 0, n = 60, seed = 11) {
  d <- make_sim(n = n, s = s, seed = seed)
  list(
    d = d, X = d$X, r = d$obs_times, Z = as.matrix(d$covariates),
    dl = d$delta, n = length(unique(d$id))
  )
}

# The R reference for the risk set at time t: weights over rows k.
rs_w <- function(t, X, r, h, weight, one_sided) {
  (t <= X) * kernel_weights(t - r, h, weight, one_sided)
}

test_that("kee_cox solves the score built from its own weight", {
  for (cfg in list(list(w = gauss_half, os = TRUE, h = 0.1),
                   list(w = NULL, os = FALSE, h = 0.2))) {
    e <- rs_setup()
    fit <- kee_cox(survival::Surv(X, delta) ~ covariates, e$d, id = id,
                   obs_times = obs_times, h = cfg$h, one_sided = cfg$os,
                   weight = cfg$w)
    w <- resolve_weight(cfg$w, cfg$os)
    kv <- kernel_weights(e$X - e$r, cfg$h, w, cfg$os)
    eb <- exp(as.vector(e$Z %*% coef(fit)))
    U <- 0
    for (i in which(e$dl == 1 & kv != 0)) {
      a <- rs_w(e$X[i], e$X, e$r, cfg$h, w, cfg$os) * eb
      U <- U + kv[i] * (e$Z[i, ] - colSums(a * e$Z) / sum(a))
    }
    expect_lt(max(abs(U / e$n)), 1e-6)
  }
})

test_that("kee_additive is the closed form built from its own weight", {
  for (cfg in list(list(w = gauss_half, os = TRUE, h = 0.1),
                   list(w = NULL, os = FALSE, h = 0.2))) {
    e <- rs_setup(s = 1)
    fit <- kee_additive(survival::Surv(X, delta) ~ covariates, e$d, id = id,
                        obs_times = obs_times, h = cfg$h, lq_nodes = 32,
                        one_sided = cfg$os, weight = cfg$w)
    w <- resolve_weight(cfg$w, cfg$os)
    kv <- kernel_weights(e$X - e$r, cfg$h, w, cfg$os)
    tau <- max(e$X)
    lq <- gaussquad::legendre.quadrature.rules(32)[[32]]
    A <- matrix(0, 2, 2)
    for (q in seq_along(lq$x)) {
      tt <- 0.5 * tau * (lq$x[q] + 1)
      a <- rs_w(tt, e$X, e$r, cfg$h, w, cfg$os)
      if (sum(a) == 0) next
      zbar <- colSums(a * e$Z) / sum(a)
      A <- A + 0.5 * tau * lq$w[q] *
        crossprod(sweep(e$Z, 2, zbar) * a, e$Z)
    }
    B <- 0
    for (i in which(e$dl == 1 & kv != 0)) {
      a <- rs_w(e$X[i], e$X, e$r, cfg$h, w, cfg$os)
      B <- B + kv[i] * (e$Z[i, ] - colSums(a * e$Z) / sum(a))
    }
    expect_equal(unname(coef(fit)), as.vector(solve(A, B)), tolerance = 1e-8)
  }
})

test_that("the skmle bread uses the fit's weight and alpha at the risk-set time", {
  # s = 0 so that H''/H = exp(.) does not cancel: at s = 0.5 it is identically
  # 1 and a baseline evaluated at the wrong time would go unseen.
  e <- rs_setup(s = 0)
  h <- 0.1
  fit <- skmle(survival::Surv(X, delta) ~ covariates, e$d, id = id,
               obs_times = obs_times, s = 0, h = h, weight = gauss_half)
  tau <- max(e$X)
  bs <- as.matrix(splines::ns(e$X, knots = tau * (1:3) / 4, intercept = TRUE,
                              Boundary.knots = c(0, tau)))
  al <- as.vector(bs %*% fit$gamma)
  zb <- as.vector(e$Z %*% coef(fit))
  kv <- kernel_weights(e$X - e$r, h, gauss_half, TRUE)
  A <- matrix(0, 2, 2)
  for (i in which(e$dl == 1 & kv != 0)) {
    a <- rs_w(e$X[i], e$X, e$r, h, gauss_half, TRUE) * exp(al[i] + zb)
    S1 <- colSums(a * e$Z) / sum(a)
    S2 <- crossprod(e$Z * a, e$Z) / sum(a)
    A <- A + (S2 - tcrossprod(S1)) * kv[i]
  }
  expect_equal(unname(fit$A), A / e$n, tolerance = 1e-8)
})

test_that("the full kernel reaches the skmle meat as well as the fit", {
  # calc_B used to require a positive lag whatever one_sided said, so rows the
  # objective and the bread kept were missing from the meat.
  # The full scheme drops I(R <= X) and I(R <= t) everywhere (Remark 7), so
  # the meat is rebuilt here with the full kernel in every term. s = 0, where
  # H^(1) = 1 and H' = H''/H = exp.
  e <- rs_setup(s = 0)
  h <- 0.2
  fit <- skmle(survival::Surv(X, delta) ~ covariates, e$d, id = id,
               obs_times = obs_times, s = 0, h = h, lq_nodes = 32,
               one_sided = FALSE)
  w <- resolve_weight(NULL, FALSE)
  lag <- e$X - e$r
  kv <- kernel_weights(lag, h, w, FALSE)
  expect_true(any(e$dl == 1 & lag <= 0 & kv != 0))

  tau <- max(e$X)
  knots <- tau * (1:3) / 4
  basis <- function(t) {
    as.matrix(splines::ns(t, knots = knots, intercept = TRUE,
                          Boundary.knots = c(0, tau)))
  }
  al <- as.vector(basis(e$X) %*% fit$gamma)
  zb <- as.vector(e$Z %*% coef(fit))
  sid <- as.integer(factor(e$d$id))
  bb1 <- matrix(0, e$n, 2)
  for (i in which(e$dl == 1 & kv != 0)) {
    a <- rs_w(e$X[i], e$X, e$r, h, w, FALSE) * exp(al[i] + zb)
    bb1[sid[i], ] <- bb1[sid[i], ] + (colSums(a * e$Z) / sum(a) - e$Z[i, ]) * kv[i]
  }
  lq <- gaussquad::legendre.quadrature.rules(32)[[32]]
  tts <- 0.5 * tau * (lq$x + 1)
  al_tt <- as.vector(basis(tts) %*% fit$gamma)
  bb2 <- matrix(0, e$n, 2)
  for (q in seq_along(tts)) {
    k_tt <- rs_w(tts[q], e$X, e$r, h, w, FALSE)
    a <- k_tt * exp(al_tt[q] + zb)
    if (sum(a) == 0) next
    r1 <- a * sweep(-e$Z, 2, colSums(a * e$Z) / sum(a), "+")
    bb2 <- bb2 + 0.5 * tau * lq$w[q] * rowsum(r1, sid, reorder = TRUE)
  }
  expect_equal(unname(fit$B), crossprod(bb1 - bb2) / e$n, tolerance = 1e-8)
})

test_that("feasibility constraints follow the weight's reach", {
  lag <- c(-0.3, -0.1, 0, 0.05, 0.1, 0.2, 0.29)
  h <- 0.1
  # The Epanechnikov kernels reproduce the earlier |lag| <= h rule exactly.
  for (os in c(TRUE, FALSE)) {
    w <- resolve_weight(NULL, os)
    kv <- kernel_weights(lag, h, w, os)
    expect_identical(constraint_rows(lag, kv, h, w, os),
                     abs(lag) <= h & (!os | lag > 0))
  }
  # A weight reaching to 3h constrains every row it weights.
  kv <- kernel_weights(lag, h, gauss_half, TRUE)
  rows <- constraint_rows(lag, kv, h, gauss_half, TRUE)
  expect_true(all(rows[kv != 0]))
  expect_identical(rows, lag > 0 & lag <= 3 * h)
})

test_that("one_sided keeps its position from before weight existed", {
  e <- rs_setup()
  f <- survival::Surv(X, delta) ~ covariates
  expect_equal(
    coef(kee_cox(f, e$d, id, obs_times, 0.2, FALSE)),
    coef(kee_cox(f, e$d, id = id, obs_times = obs_times, h = 0.2,
                 one_sided = FALSE))
  )
})

test_that("an s = 1 fit that stalls is restarted, not returned", {
  skip_on_cran()
  # Replicate 108 of the BigRed200 meat study. SLSQP stalled on the floor of
  # the cumulative-hazard integrand (status -1, log-likelihood -1.057) with an
  # event row on the boundary and SEs of 1e-16, which skmle() returned without
  # comment. Restarted from the quadrature-constrained solution, the paper's
  # problem converges to -0.832.
  s <- 1
  n <- 200L
  set.seed(1e6 * s + 1e4 * n + 108)
  d <- sim_skmle_data(n = n, mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
    mu_bar = 8,
    alpha = function(tt) {
      (1 + s) / 2 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25))))
    },
    beta = c(1, -0.5), s = s, cen = 0.7)
  fit <- expect_no_warning(skmle(survival::Surv(X, delta) ~ covariates, d,
    id = id, obs_times = obs_times, s = s, h = n^(-0.4)))
  expect_gte(fit$convergence, 0)
  expect_gt(fit$loglik, -1)
  expect_true(all(sqrt(diag(fit$var)) > 0.05))
})

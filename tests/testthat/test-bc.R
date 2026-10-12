fake_fits <- function(h_grid, beta0, slope, p = length(beta0), n = 50) {
  set.seed(1)
  sc <- matrix(rnorm(n * p), n, p)
  lapply(h_grid, function(h) {
    A <- diag(p) * (1 + h)
    list(coefficients = setNames(beta0 + slope * h, paste0("z", seq_len(p))),
         A = A, scores = sc * sqrt(h), n = n, convergence = 0L,
         var = solve(A) %*% (crossprod(sc * sqrt(h)) / n) %*% solve(A) / n)
  })
}

test_that("bc_combine recovers beta0 exactly when the bias is linear in h", {
  h <- exp(seq(log(0.05), log(1), length.out = 8))
  fits <- fake_fits(h, beta0 = c(1, -0.5), slope = c(-0.4, 0.2))
  est <- bc_combine(fits, h)
  expect_equal(unname(est$coefficients), c(1, -0.5), tolerance = 1e-10)
})

test_that("extrapolation weights sum to one and annihilate h", {
  h <- exp(seq(log(0.05), log(1), length.out = 10))
  est <- bc_combine(fake_fits(h, c(0.3, 0.1), c(1, -1)), h)
  expect_equal(unname(colSums(est$weights)), c(1, 1), tolerance = 1e-10)
  expect_equal(unname(colSums(est$weights * h)), c(0, 0), tolerance = 1e-10)
  expect_true(all(est$weights[!h %in% est$window, ] == 0))
})

test_that("the joint covariance's diagonal blocks are each fit's variance", {
  dat <- make_sim(80, s = 0, seed = 5)
  # A five-point grid often selects an endpoint, which warns.
  bc <- suppressWarnings(
    skmle_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0,
             h_grid = c(0.1, 0.2, 0.3, 0.5, 0.8), quiet = TRUE)
  )
  f <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.3)
  expect_equal(unname(bc$cov_all[5:6, 5:6]), unname(f$var), tolerance = 1e-8)
  expect_s3_class(bc, "skmle_bc")
  expect_true(all(is.finite(bc$var)))
  expect_true(bc$h %in% c(0.1, 0.2, 0.3, 0.5, 0.8))
})

test_that("skmle's B is the mean outer product of its scores", {
  dat <- make_sim(60, s = 0, seed = 2)
  f <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, h = 0.3)
  expect_equal(dim(f$scores), c(f$n, 2L))
  expect_equal(f$B, crossprod(f$scores) / f$n)
})

test_that("skmle_bc rejects a grid with fewer than four bandwidths", {
  dat <- make_sim(40, s = 0, seed = 1)
  expect_error(
    skmle_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, h_grid = c(0.1, 0.2, 0.3), quiet = TRUE),
    "at least 4"
  )
})

test_that("a fit whose variance has collapsed is not used", {
  h <- exp(seq(log(0.05), log(1), length.out = 8))
  fits <- fake_fits(h, c(1, -0.5), c(-0.4, 0.2))
  fits[[3]]$var <- fits[[3]]$var * 1e-30
  fits[[5]]$convergence <- -1
  for (j in setdiff(seq_along(fits), 5)) fits[[j]]$convergence <- 4
  ok <- bc_usable(fits)
  expect_equal(which(!ok), c(3L, 5L))
})

test_that("n_boot = 0 reports the fixed-window variance", {
  h <- exp(seq(log(0.05), log(1), length.out = 8))
  est <- bc_combine(fake_fits(h, c(1, -0.5), c(-0.4, 0.2)), h, n_boot = 0)
  expect_equal(est$var, est$var_fixed)
})

test_that("the selection-adjusted variance is reproducible and not smaller", {
  dat <- make_sim(120, s = 0, seed = 8)
  run <- function() suppressWarnings(
    skmle_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, s = 0, n_h = 8, n_boot = 4000,
             seed = 3, quiet = TRUE)
  )
  a <- run()
  b <- run()
  expect_equal(a$var, b$var)
  # On this data set the selection adds spread; allow Monte Carlo error.
  expect_true(all(diag(a$var) >= 0.95 * diag(a$var_fixed)))
})

test_that("skmle_bc rejects n_boot = 1", {
  dat <- make_sim(40, s = 0, seed = 1)
  expect_error(
    skmle_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
             obs_times = obs_times, n_boot = 1, quiet = TRUE),
    "n_boot"
  )
})

test_that("degree 2 recovers beta0 exactly when the bias is quadratic in h", {
  h <- exp(seq(log(0.05), log(1), length.out = 10))
  fits <- fake_fits(h, beta0 = c(1, -0.5), slope = c(0, 0))
  for (j in seq_along(fits)) {
    fits[[j]]$coefficients <- fits[[j]]$coefficients +
      c(-0.4, 0.2) * h[j] + c(0.3, -0.6) * h[j]^2
  }
  est <- bc_combine(fits, h, n_boot = 0, degree = 2)
  expect_equal(unname(est$coefficients), c(1, -0.5), tolerance = 1e-10)
  expect_equal(unname(colSums(est$weights * h^2)), c(0, 0), tolerance = 1e-10)
})

test_that("each estimator's influence contributions reproduce its sandwich", {
  dat <- make_sim(80, s = 0, seed = 3)
  f1 <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
              obs_times = obs_times, s = 0, h = 0.3)
  f2 <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
                obs_times = obs_times, h = 0.3)
  d <- sim_async_data(n = 60)
  f3 <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
  for (f in list(f1, f2, f3)) {
    psi <- bc_psi(f)
    expect_equal(nrow(psi), f$n)
    expect_equal(unname(crossprod(psi)), unname(f$var), tolerance = 1e-8)
  }
})

test_that("kee_cox_bc and kee_async_bc run and pass their arguments on", {
  dat <- make_sim(100, s = 0, seed = 4)
  b1 <- suppressWarnings(
    kee_cox_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
               obs_times = obs_times, n_h = 8, n_boot = 200, seed = 1,
               quiet = TRUE)
  )
  expect_s3_class(b1, "kee_bc")
  expect_equal(b1$degree, 1L)
  expect_true(all(is.finite(b1$var)))
  expect_equal(unname(b1$fit$coefficients),
               unname(b1$grid$estimate[b1$grid$h == b1$h]))

  d <- sim_async_data(n = 100)
  b2 <- suppressWarnings(
    kee_async_bc(d$y, d$x, y ~ x, id = id, time = time, n_h = 8,
                 n_boot = 200, seed = 1, quiet = TRUE)
  )
  expect_s3_class(b2, "kee_bc")
  expect_equal(b2$degree, 1L)
  expect_true(all(is.finite(b2$var)))
  b3 <- suppressWarnings(
    kee_async_bc(d$y, d$x, y ~ x, id = id, time = time, n_h = 8, degree = 2,
                 n_boot = 0, one_sided = TRUE, quiet = TRUE)
  )
  expect_equal(b3$degree, 2L)
  expect_true(b3$fit$one_sided)
})

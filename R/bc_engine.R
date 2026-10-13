# Shared machinery for skmle_bc(), kee_cox_bc() and kee_async_bc() ----------
#
# Each estimator is refitted over a grid of bandwidths. Everything after that
# depends on the fits only through their coefficients and their per-subject
# influence contributions psi (n x p, beta_h - beta ~ sum_i psi_i), so one
# engine serves all three: the joint covariance of the grid estimates is
# crossprod(psi_j, psi_m), and the bias fit, the selection, the extrapolation
# and the selection-adjusted variance are functions of the grid estimates and
# that covariance alone.

# Influence contributions of one fit, rows aligned by subject code.
#   skmle:     var = A^-1 (S'S / n) A^-1 / n      -> psi = S A^-T / n
#   kee_cox:   var = W^-1 (U'U / n^2) W^-1        -> psi = U W^-T / n
#   kee_async: var = A^-1 (U'U) A^-1              -> psi = U A^-T
bc_psi <- function(fit) {
  if (!is.null(fit$psi)) return(fit$psi)
  out <- tryCatch({
    if (inherits(fit, "kee_async")) {
      fit$scores %*% t(solve(fit$A))
    } else if (inherits(fit, "kee")) {
      fit$scores %*% t(solve(fit$W)) / fit$n
    } else {
      fit$scores %*% t(solve(fit$A)) / fit$n
    }
  }, error = function(e) NULL)
  if (is.null(out) || !all(is.finite(out))) return(NULL)
  out
}

# Which fits can enter the combination. Besides failed fits, failed
# optimisations and non-finite variances, a fit can return a sandwich that has
# collapsed to rounding error (variances near 1e-30 were seen in simulation
# with status 4). Its inverse-variance weight then swamps the regressions, so
# a variance below 1e-8 times the grid's median is treated as unusable.
bc_usable <- function(fits) {
  ok <- vapply(fits, function(f) {
    !is.null(f) && all(is.finite(f$var)) && f$convergence >= 0 &&
      all(diag(f$var) > 0)
  }, logical(1))
  if (any(ok)) {
    v <- vapply(fits, function(f) {
      if (!is.null(f) && all(is.finite(f$var))) min(diag(f$var)) else NA_real_
    }, numeric(1))
    ok <- ok & !is.na(v) & v > 1e-8 * stats::median(v[ok])
  }
  ok
}

bc_check_n_boot <- function(n_boot) {
  if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) ||
      n_boot == 1 || n_boot < 0) {
    stop("'n_boot' must be 0 or at least 2", call. = FALSE)
  }
  invisible(TRUE)
}

bc_check_degree <- function(degree) {
  if (!is.numeric(degree) || length(degree) != 1 || !degree %in% c(1, 2)) {
    stop("'degree' must be 1 or 2", call. = FALSE)
  }
  as.integer(degree)
}

# Fit `fit_call` at every bandwidth in h_grid, drop unusable fits, and combine.
# A fit that stops with an error is dropped like one that fails to converge.
bc_fit_grid <- function(fit_call, env, h_grid, n_boot, seed, quiet, degree) {
  h_grid <- sort(unique(as.numeric(h_grid)))
  if (length(h_grid) < 4 || any(!is.finite(h_grid)) || any(h_grid <= 0)) {
    stop("'h_grid' must hold at least 4 distinct positive bandwidths",
         call. = FALSE)
  }
  fits <- lapply(seq_along(h_grid), function(j) {
    if (!quiet) cat(sprintf("fitting h = %.4g (%d/%d)\n", h_grid[j], j,
                            length(h_grid)))
    fit_call$h <- h_grid[j]
    f <- tryCatch(eval(fit_call, envir = env), error = function(e) NULL)
    if (!is.null(f)) f$psi <- bc_psi(f)
    if (!is.null(f) && is.null(f$psi)) f <- NULL
    f
  })

  usable <- bc_usable(fits)
  if (!all(usable)) {
    warning("dropping bandwidth(s) ", paste(format(h_grid[!usable], digits = 3),
                                           collapse = ", "),
            ": the fit failed or its variance is unusable", call. = FALSE)
    fits <- fits[usable]
    h_grid <- h_grid[usable]
  }
  if (length(h_grid) < 4) {
    stop("fewer than 4 bandwidths have usable fits", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  est <- bc_combine(fits, h_grid, n_boot = n_boot, degree = degree)
  if (est$h %in% range(h_grid)) {
    warning("the selected bandwidth is an endpoint of 'h_grid'; widen the grid",
            call. = FALSE)
  }
  fit_h <- fits[[match(est$h, h_grid)]]
  fit_h$psi <- NULL
  list(est = est, fit = fit_h)
}

# The arithmetic on a list of fits at bandwidths h_grid. Kept apart from the
# fitting so it can be checked on its own.
#
# degree = 1: bias c1 h + e / h, extrapolation linear in h. Removes the
#   first-order bias of a one-sided kernel, and of a two-sided kernel on a
#   rough covariate path (a jump process, or a covariance like exp(-|t - s|)).
#   On a smooth path a two-sided kernel's bias is c2 h^2; linear extrapolation
#   leaves O(h^2) at an h-hat of order n^(-1/3), which is o((n h)^(-1/2)). In
#   the kee_async simulation (dev/sim_kee) degree 1 covered at 0.93-0.95 under
#   both covariances with the smallest RMSE.
# degree = 2: bias c1 h + c2 h^2 + e / h, extrapolation quadratic in h.
#   Removes the h^2 term too, with a larger variance.
bc_combine <- function(fits, h_grid, n_boot = 2000L, degree = 1L) {
  H <- length(h_grid)
  p <- length(fits[[1]]$coefficients)
  cn <- names(fits[[1]]$coefficients)

  B <- t(vapply(fits, function(f) unname(f$coefficients), numeric(p)))
  B <- matrix(B, H, p)
  V <- matrix(t(vapply(fits, function(f) diag(f$var), numeric(p))), H, p)

  # Joint covariance, ordered (h_1: coef 1..p, h_2: coef 1..p, ...).
  psi <- lapply(fits, bc_psi)
  cov_all <- matrix(0, H * p, H * p)
  for (j in seq_len(H)) {
    for (k in j:H) {
      blk <- crossprod(psi[[j]], psi[[k]])
      rj <- (j - 1) * p + seq_len(p)
      rk <- (k - 1) * p + seq_len(p)
      cov_all[rj, rk] <- blk
      cov_all[rk, rj] <- t(blk)
    }
  }

  # With the variances held at their estimates, every step below is linear in
  # the grid estimates: the bias fit is a fixed projection, and for each
  # candidate bandwidth the extrapolation is a fixed set of weights. Building
  # those maps once lets the selection be rerun on simulated grids cheaply.
  # The regressions run on h / h_ref, which leaves every fitted value and
  # intercept unchanged but keeps the design well conditioned whatever the
  # unit of h: with h in days, (1, h, 1/h) spans about eight orders of
  # magnitude, and together with inverse-variance weights that differ by
  # several orders across the grid the normal equations became numerically
  # singular. For the same reason the weighted fits go through a QR
  # decomposition of the weighted design rather than its normal equations.
  h_ref <- exp(mean(log(h_grid)))
  powers <- function(h) outer(h / h_ref, seq_len(degree), `^`)
  # Rows of P give the weighted least-squares coefficients as linear maps of
  # the response: coef = P %*% y.
  wls_map <- function(X, W) {
    sw <- sqrt(W / max(W))
    # qr()'s default rank tolerance (1e-7) misreads a full-rank design whose
    # weights span many orders of magnitude as deficient (ADNI amy:plasma
    # GFAP_Q); only a design that is singular to rounding error is refused.
    qx <- qr(X * sw, tol = 1e-12)
    if (qx$rank < ncol(X)) stop("rank-deficient bias regression", call. = FALSE)
    qr.coef(qx, diag(sw, length(sw)))
  }

  # Plug-in: bias per coefficient, weighted by inverse variance.
  # bias_map[[k]] takes column k of the grid estimates to the fitted bias.
  Xb <- cbind(1, powers(h_grid), h_ref / h_grid)
  bias_map <- lapply(seq_len(p), function(k) {
    P <- wls_map(Xb, 1 / V[, k])
    Xb[, -1, drop = FALSE] %*% P[-1, , drop = FALSE]
  })
  select_h <- function(Bs) {
    # Bs: H x p x m array of grid estimates; returns the selected index per
    # slice.
    m <- dim(Bs)[3]
    mse <- matrix(rowSums(V), H, m)
    for (k in seq_len(p)) mse <- mse + (bias_map[[k]] %*% Bs[, k, ])^2
    apply(mse, 2, which.min)
  }

  # Extrapolation over [h/2, 2h]; at least degree + 1 grid points, taking the
  # nearest neighbours of h in log scale when the window holds fewer.
  # ext_w[[j]] holds the weights (H x p) used when h_grid[j] is selected.
  window_of <- function(j) {
    in_win <- h_grid >= h_grid[j] / 2 * (1 - 1e-8) &
      h_grid <= 2 * h_grid[j] * (1 + 1e-8)
    if (sum(in_win) < degree + 1) {
      in_win[order(abs(log(h_grid / h_grid[j])))[seq_len(degree + 1)]] <- TRUE
    }
    in_win
  }
  ext_w <- lapply(seq_len(H), function(j) {
    in_win <- window_of(j)
    X2 <- cbind(1, powers(h_grid[in_win]))
    a <- matrix(0, H, p)
    for (k in seq_len(p)) {
      M <- wls_map(X2, 1 / V[in_win, k])
      a[in_win, k] <- M[1, ]
    }
    a
  })

  bias <- vapply(seq_len(p), function(k) drop(bias_map[[k]] %*% B[, k]),
                 numeric(H))
  bias <- matrix(bias, H, p)
  mse <- rowSums(bias^2 + V)
  if (!any(is.finite(mse))) {
    stop("the plug-in mean squared error could not be computed on any ",
         "bandwidth", call. = FALSE)
  }
  j_sel <- which.min(mse)
  h_sel <- h_grid[j_sel]
  in_win <- window_of(j_sel)
  a <- ext_w[[j_sel]]
  beta_bc <- colSums(a * B)

  # Variance treating the selected bandwidth and window as fixed:
  # L cov_all t(L), with L picking the terms sum_j a_jk beta_{h_j, k}.
  L <- matrix(0, p, H * p)
  for (k in seq_len(p)) L[k, (seq_len(H) - 1) * p + k] <- a[, k]
  var_fixed <- L %*% cov_all %*% t(L)

  # Variance including the selection. Draw grids from N(B, cov_all), rerun the
  # selection and extrapolation on each, and take the covariance of the
  # results. In simulation the fixed-window formula understated the spread by
  # 6-9% at s = 1, because the selected bandwidth and window move with the
  # data.
  var_bc <- var_fixed
  if (n_boot > 0) {
    eg <- eigen(cov_all, symmetric = TRUE)
    root <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), nrow(cov_all))
    z <- matrix(stats::rnorm(n_boot * H * p), H * p, n_boot)
    draws <- as.vector(t(B)) + root %*% z
    Bs <- array(0, c(H, p, n_boot))
    for (k in seq_len(p)) Bs[, k, ] <- draws[(seq_len(H) - 1) * p + k, ]
    js <- select_h(Bs)
    bt <- t(vapply(seq_len(n_boot), function(b) {
      colSums(ext_w[[js[b]]] * Bs[, , b])
    }, numeric(p)))
    var_bc <- stats::cov(matrix(bt, n_boot, p))
  }

  names(beta_bc) <- cn
  dimnames(var_bc) <- dimnames(var_fixed) <- list(cn, cn)
  dimnames(a) <- list(format(h_grid, digits = 4), cn)
  lab <- paste(rep(format(h_grid, digits = 4), each = p), rep(cn, H), sep = ":")
  dimnames(cov_all) <- list(lab, lab)

  grid <- data.frame(
    h = rep(h_grid, each = p),
    term = rep(cn, H),
    estimate = as.vector(t(B)),
    std.error = sqrt(as.vector(t(V))),
    bias = as.vector(t(bias)),
    mse = rep(mse, each = p),
    stringsAsFactors = FALSE
  )

  list(coefficients = beta_bc, var = var_bc, var_fixed = var_fixed,
       h = h_sel, window = h_grid[in_win], weights = a, grid = grid,
       cov_all = cov_all, n_boot = n_boot, degree = degree)
}

#' @param x,object A bias-corrected fit from [skmle_bc()], [kee_cox_bc()] or
#'   [kee_async_bc()].
#' @param ... Ignored.
#' @rdname skmle_bc
#' @export
print.bc <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(sprintf(paste0("\nSelected h = %s (plug-in MSE); extrapolated (degree %d)",
                     " over h in [%s, %s]\n"),
              format(x$h, digits = 4), x$degree, format(min(x$window), digits = 4),
              format(max(x$window), digits = 4)))
  cat("\nBias-corrected coefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' @rdname skmle_bc
#' @export
summary.bc <- function(object, ...) {
  se_bc <- sqrt(diag(object$var))
  se_h <- sqrt(diag(object$fit$var))
  tab <- data.frame(
    estimate = object$coefficients,
    std.error = se_bc,
    z = object$coefficients / se_bc,
    p.value = 2 * stats::pnorm(-abs(object$coefficients / se_bc)),
    estimate_at_h = object$fit$coefficients,
    std.error_at_h = se_h
  )
  res <- list(call = object$call, h = object$h, window = object$window,
              coefficients = tab)
  class(res) <- "summary.bc"
  res
}

#' @rdname skmle_bc
#' @export
print.summary.bc <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(sprintf("\nSelected h = %s; extrapolated over h in [%s, %s]\n\n",
              format(x$h, digits = 4), format(min(x$window), digits = 4),
              format(max(x$window), digits = 4)))
  print(format(x$coefficients, digits = 4))
  invisible(x)
}

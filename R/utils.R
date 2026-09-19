`%||%` <- function(a, b) if (is.null(a)) b else a

#' Validate the time scale
#'
#' Times may be on any scale. The sieve basis and the cumulative-hazard
#' quadrature are built on `[0, max(X_time)]`, so nothing has to be rescaled to
#' the unit interval first; that was a convention of the original prototype
#' rather than a requirement of the method. What is required of `X_time` is
#' that it is finite and non-negative, and that follow-up has positive length.
#'
#' `obs_times_vec` is held to the weaker standard of being finite, and may be
#' negative. No basis and no quadrature rule is ever evaluated at a covariate
#' observation time: it reaches the score only through the lag
#' `X_time - obs_times_vec`, and the one-sided rule tests the sign of that lag,
#' not the sign of the time. The earlier non-negativity requirement applied a
#' constraint belonging to `X_time` to an argument that does not carry it. A
#' negative observation time means a covariate measured before the time origin,
#' which is a scientific statement about the data and not a numerical one; see
#' \dQuote{Negative observation times} in [kee_cox()].
#'
#' @param X_time Event or censoring times.
#' @param obs_times_vec Covariate observation times. May be negative.
#' @return `NULL`, invisibly. Called for the error.
#' @keywords internal
check_time_scale <- function(X_time, obs_times_vec) {
  if (anyNA(X_time) || anyNA(obs_times_vec) ||
    !all(is.finite(X_time)) || !all(is.finite(obs_times_vec))) {
    stop("event and observation times must be finite", call. = FALSE)
  }
  if (min(X_time) < 0) {
    stop(
      "event times must be non-negative; observed minimum ",
      format(min(X_time)),
      call. = FALSE
    )
  }
  if (max(X_time) <= 0) {
    stop("follow-up has zero length: all event times are 0", call. = FALSE)
  }
  invisible(NULL)
}

#' Build a sandwich covariance matrix with a reciprocal-condition guard
#'
#' Computes `A^{-1} B A^{-1} * scale`. Falls back to an `NA` matrix with an
#' informative warning if `A` is numerically singular; this surfaces the
#' underlying numerical issue rather than the generic `solve()` error.
#'
#' @param A Square "bread" matrix.
#' @param B Symmetric "meat" matrix of the same size as `A`.
#' @param scale Scalar multiplier (typically `1/n` or `1`).
#' @param what Short label used in the warning message.
#'
#' @return A `ncol(A)` by `ncol(A)` matrix: either the sandwich
#'   `A^{-1} B A^{-1} * scale` or a same-shaped matrix of `NA_real_`.
#' @noRd
safe_sandwich <- function(A, B, scale = 1, what = "variance") {
  rc <- tryCatch(rcond(A), error = function(e) 0)
  if (!is.finite(rc) || rc < .Machine$double.eps) {
    warning(sprintf(
      "%s estimation is unreliable: 'A' matrix is near-singular (reciprocal condition %.3g). Returning NA covariance.",
      what, rc
    ), call. = FALSE)
    return(matrix(NA_real_, ncol(A), ncol(A)))
  }
  A_inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(A_inv)) {
    warning(sprintf(
      "%s estimation failed: solve(A) raised an error despite finite reciprocal condition.",
      what
    ), call. = FALSE)
    return(matrix(NA_real_, ncol(A), ncol(A)))
  }
  A_inv %*% B %*% A_inv * scale
}


#' Epanechnikov kernel and the row weights built from it
#'
#' The kernel was written out inline in `skmle()`, `kee_cox()` and
#' `kee_additive()`. Keeping one copy matters now that the half/full choice is
#' a user-facing argument: three inline copies are three places for the switch
#' to be forgotten.
#'
#' @param u Numeric vector or matrix of standardised lags.
#' @return `epan_kernel()` returns \eqn{0.75(1 - u^2)_+}, shape preserved.
#' @keywords internal
epan_kernel <- function(u) {
    d <- dim(u)
    val <- pmax((1 - as.numeric(u)^2) * 0.75, 0)
    dim(val) <- d
    val
}

#' @param lag Raw time difference `t - r`, a vector or a matrix. Matrices keep
#'   their shape, which the sieve quadrature relies on.
#' @param h Positive bandwidth.
#' @param one_sided Logical. When `TRUE` (the default) rows with a non-positive
#'   lag receive zero weight, which is the risk-set restriction of a hazard
#'   model: only covariate observations before the time inform it. `FALSE`
#'   smooths from both sides.
#' @return `kernel_weights()` returns the scaled weights `W(lag/h)/h`.
#' @rdname epan_kernel
#' @keywords internal
kernel_weights <- function(lag, h, one_sided = TRUE) {
    kv <- epan_kernel(lag / h) / h
    if (one_sided) kv <- kv * (lag > 0)
    kv
}

#' The Epanechnikov kernel described as a polynomial
#'
#' The cross-validation loop rebuilds the weights for every candidate bandwidth
#' inside C++ and cannot call an R function to do it, so the weight travels as
#' four plain values instead: the coefficients of \eqn{\sum_j c_j u^j}, the
#' support, and whether the shape functions are \eqn{|u|^j} rather than
#' \eqn{u^j}. For \eqn{0.75(1 - u^2)} the coefficients are `c(0.75, 0, -0.75)`
#' and the support is `[0, 1]` for the half kernel, `[-1, 1]` for the full one.
#' `calc_kerfun()` in `src/skmle_cpp.cpp` reads them back.
#'
#' This is the second description of a kernel the package already has in
#' [kernel_weights()], so `test-cv-weightspec.R` holds the two against each
#' other. Two descriptions that nothing compares is how they drift apart.
#'
#' @param one_sided Logical; `TRUE` for the half kernel.
#' @return A list with `coef`, `a`, `b` and `mirror`.
#' @keywords internal
epan_weight_spec <- function(one_sided = TRUE) {
    list(
        coef = c(0.75, 0, -0.75),
        a = if (one_sided) 0 else -1,
        b = 1,
        mirror = FALSE
    )
}


# Data-driven default bandwidths.
#
# A student should not have to invent a number to get a first fit.  These are
# rules of thumb.  Cross-validation does the job properly, and the message that
# accompanies the default says so and names the function.

#' Default bandwidth for the survival estimators
#'
#' Geometric midpoint of the grid `skmle_cv()` searches, which is built from the
#' observed lags between covariate observation times and event times.
#'
#' @param X_time,obs_times Event/censoring and covariate observation times.
#' @param n Number of subjects.
#' @return A positive bandwidth.
#' @keywords internal
default_bandwidth_surv <- function(X_time, obs_times, n) {
  pos <- X_time - obs_times
  pos <- pos[pos > 0 & is.finite(pos)]
  if (!length(pos)) {
    stop(
      "no covariate observation precedes an event time, so no bandwidth can ",
      "be chosen automatically; supply 'h'",
      call. = FALSE
    )
  }
  tau <- max(X_time)
  lo <- max(min(pos), tau * n^(-0.6))
  hi <- min(max(pos), tau * n^(-0.3))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) lo else sqrt(lo * hi)
}

#' Default bandwidth for the asynchronous estimators
#'
#' Geometric midpoint of the grid `kee_async_cv()` searches,
#' \eqn{2 (Q_3 - Q_1) n^{-1/2}}, so it adapts to whatever units `time` is in.
#'
#' @param times Pooled observation times from both tables.
#' @param n Number of subjects.
#' @return A positive bandwidth.
#' @keywords internal
default_bandwidth_async <- function(times, n) {
  iqr <- diff(stats::quantile(times, c(0.25, 0.75), names = FALSE, na.rm = TRUE))
  if (!is.finite(iqr) || iqr <= 0) iqr <- diff(range(times, na.rm = TRUE))
  if (!is.finite(iqr) || iqr <= 0) {
    stop("observation times are constant; supply 'h'", call. = FALSE)
  }
  2 * iqr * n^(-0.5)
}

#' Tell the user which bandwidth was chosen for them
#'
#' @param h The chosen bandwidth.
#' @param cv_fun Name of the cross-validation function to point at.
#' @return `h`, invisibly.
#' @keywords internal
announce_bandwidth <- function(h, cv_fun) {
  message(sprintf(
    paste0(
      "'h' not supplied. Using h = %s, read off the observation times as a\n",
      "rule of thumb. See %s() to choose the bandwidth from the data."
    ),
    format(h, digits = 3), cv_fun
  ))
  invisible(h)
}

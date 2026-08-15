`%||%` <- function(a, b) if (is.null(a)) b else a

#' Validate that event and observation times are on the unit interval
#'
#' The sieve approximation of the baseline hazard uses a natural cubic
#' B-spline with `Boundary.knots = c(0, 1)`. Times outside this interval
#' would rely on extrapolation and produce a sieve that is not the one
#' implied by the paper. Users are expected to rescale their time axis
#' beforehand.
#'
#' @param X_time Numeric vector of event/censoring times.
#' @param obs_times_vec Numeric vector of covariate observation times.
#'
#' @return `NULL`, invisibly. Called for the side effect of raising an
#'   informative error when either vector extends outside `[0, 1]`.
#' @noRd
check_time_scale <- function(X_time, obs_times_vec) {
  rng_X <- range(X_time, na.rm = TRUE)
  rng_O <- range(obs_times_vec, na.rm = TRUE)
  if (rng_X[1] < 0 || rng_X[2] > 1 || rng_O[1] < 0 || rng_O[2] > 1) {
    stop(
      "Event and observation times must lie in [0, 1]. ",
      "Observed range for event/censoring times: [",
      format(rng_X[1]), ", ", format(rng_X[2]), "]; ",
      "for observation times: [",
      format(rng_O[1]), ", ", format(rng_O[2]), "]. ",
      "Rescale before calling (e.g. divide by max follow-up).",
      call. = FALSE
    )
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


# Data-driven default bandwidths.
#
# A student should not have to invent a number to get a first fit.  These are
# rules of thumb, not a substitute for cross-validation, and the message that
# accompanies them says so and names the function that does it properly.

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
  lo <- max(min(pos), n^(-0.6))
  hi <- min(max(pos), n^(-0.3))
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
      "'h' not supplied. Using h = %s, a rule-of-thumb value read off the\n",
      "observation times. It is a starting point, not a tuned choice --\n",
      "see %s() to select the bandwidth from the data."
    ),
    format(h, digits = 3), cv_fun
  ))
  invisible(h)
}

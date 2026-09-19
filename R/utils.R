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


#' Row weights from a weight function
#'
#' The one place a weight is turned into the numbers the estimators use, so
#' there is one description of any given kernel rather than one per call site.
#' `epan_kernel()` is kept because three fitting functions used to write the
#' Epanechnikov out inline and older code may still reach for it.
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
#' @param weight A weight function of the standardised lag, or `NULL` for the
#'   Epanechnikov that `one_sided` implies; see [kernel-weights]. `NULL` is the
#'   default for the same reason it is the default at every entry point: a
#'   literal `w_epan_half()` here silently turns `one_sided = FALSE` into a
#'   one-sided fit, because the support stays `[0, 1]` and the negative lags
#'   are zeroed by the weight rather than kept. The package's own test suite
#'   caught that, which is the argument for not having written it twice.
#' @param one_sided Logical. When `TRUE` rows with a non-positive lag receive
#'   zero weight, which is the risk-set restriction of a hazard model: only
#'   covariate observations before the time inform it. It is applied after the
#'   weight has been evaluated, because it is a modelling choice rather than a
#'   property of the kernel.
#' @return `kernel_weights()` returns the scaled weights `W(lag/h)/h`.
#' @rdname epan_kernel
#' @keywords internal
kernel_weights <- function(lag, h, weight = NULL, one_sided = TRUE) {
    weight <- resolve_weight(weight, one_sided)
    kv <- weight(lag / h) / h
    # A weight is allowed to ignore the shape of its argument only if it is
    # already scalar-shaped; resolve_weight() checks the matrix case up front,
    # so anything arriving here that lost its dim is a weight that slipped
    # through a path which did not resolve. Restore rather than fail silently.
    if (!identical(dim(kv), dim(lag))) dim(kv) <- dim(lag)
    if (one_sided) kv <- kv * (lag > 0)
    kv
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

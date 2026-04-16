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

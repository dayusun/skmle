# Standard model generics.
#
# `coef()` worked only by accident before this file existed: the default method
# happens to find a `$coefficients` element.  `vcov()` did not, so `confint()`
# failed on every fitted object in the package even though the covariance matrix
# was sitting in `$var`.

#' Extract the covariance matrix of a fitted model
#'
#' @param object A fitted `skmle`, `kee` or `kee_td` object.
#' @param ... Unused.
#' @return For `skmle` and `kee`, the estimated covariance matrix of the
#'   regression coefficients. For `kee_td`, see `time`.
#' @importFrom stats vcov nobs confint setNames
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 120)
#' fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
#' vcov(fit)
#' @name vcov.skmle
#' @export
vcov.skmle <- function(object, ...) object$var

#' @rdname vcov.skmle
#' @export
vcov.kee <- function(object, ...) object$var

#' @param time For a `kee_td` fit, the target time whose covariance matrix is
#'   wanted. Matched to the nearest fitted target time. When `NULL` the full
#'   list of matrices is returned, one per target time.
#' @rdname vcov.skmle
#' @export
vcov.kee_td <- function(object, time = NULL, ...) {
    if (is.null(time)) {
        return(stats::setNames(object$var, format(object$times)))
    }
    if (!is.numeric(time) || length(time) != 1L || is.na(time)) {
        stop("'time' must be a single number or NULL", call. = FALSE)
    }
    object$var[[which.min(abs(object$times - time))]]
}


#' Number of subjects contributing to a fit
#'
#' The asymptotics of every estimator in this package are in the number of
#' subjects. Rows within a subject are repeated measurements of one trajectory
#' and carry far less than one observation's worth of information, so `nobs()`
#' returns the subject count, which is what `summary()` prints as `n`.
#'
#' @param object A fitted `skmle`, `kee` or `kee_td` object.
#' @param ... Unused.
#' @return An integer, the number of subjects.
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 120)
#' fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
#' nobs(fit)
#' @name nobs.skmle
#' @export
nobs.skmle <- function(object, ...) as.integer(object$n)

#' @rdname nobs.skmle
#' @export
nobs.kee <- function(object, ...) as.integer(object$n)

#' @rdname nobs.skmle
#' @export
nobs.kee_td <- function(object, ...) as.integer(object$n)


#' Wald confidence intervals for a coefficient curve
#'
#' Pointwise Wald intervals at each target time. They are pointwise, not
#' simultaneous, and are not corrected for smoothing bias, so they cover
#' \eqn{E\hat\beta(t)} rather than \eqn{\beta(t)}; see [kee_async_td()].
#'
#' For `skmle` and `kee` fits no method is needed: [stats::confint.default()]
#' works once [vcov()] is available.
#'
#' @param object A `kee_td` object.
#' @param parm Coefficients to report, by name or position. Defaults to all.
#' @param level Confidence level.
#' @param ... Unused.
#' @return A data frame with one row per (target time, coefficient) and columns
#'   `time`, `term`, `estimate`, `se`, and the two interval endpoints.
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
#' fit <- kee_async_td(d$y, d$x,
#'   y ~ x, id = id, time = time,
#'   times = c(0.3, 0.5, 0.7), h = 0.3
#' )
#' confint(fit)
#' @importFrom stats qnorm
#' @export
confint.kee_td <- function(object, parm, level = 0.95, ...) {
    nms <- colnames(object$coefficients)
    if (missing(parm)) {
        parm <- nms
    } else if (is.numeric(parm)) parm <- nms[parm]
    parm <- parm[!is.na(match(parm, nms))]
    if (!length(parm)) stop("'parm' selects no coefficient", call. = FALSE)
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
        stop("'level' must be a single number in (0, 1)", call. = FALSE)
    }

    a <- (1 - level) / 2
    zq <- qnorm(1 - a)
    est <- object$coefficients[, parm, drop = FALSE]
    se <- object$se[, parm, drop = FALSE]
    out <- data.frame(
        time = rep(object$times, times = length(parm)),
        term = rep(parm, each = length(object$times)),
        estimate = as.numeric(est),
        se = as.numeric(se),
        lower = as.numeric(est - zq * se),
        upper = as.numeric(est + zq * se),
        stringsAsFactors = FALSE
    )
    names(out)[5:6] <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
    out
}


#' Plot a cross-validation curve
#'
#' The held-out loss against the candidate bandwidths, with the selected value
#' marked. Worth a glance every time: if the curve is still falling at an
#' endpoint of the grid, the "selected" bandwidth is where you stopped looking
#' rather than a minimum, and the grid should be widened.
#'
#' @param x A `cv.skmle` or `cv.kee_async` object.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly. Called for the plot.
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 150)
#' cv <- d$y |>
#'   kee_async_cv(d$x, y ~ x,
#'     id = id, time = time,
#'     h_grid = c(0.1, 0.2, 0.3, 0.5), K = 3, seed = 1, quiet = TRUE
#'   )
#' plot(cv)
#' @importFrom graphics abline points
#' @name plot.cv.skmle
#' @export
plot.cv.kee_async <- function(x, ...) .plot_cv(x$cv_results, x$h_cv, ...)

#' @rdname plot.cv.skmle
#' @export
plot.cv.skmle <- function(x, ...) .plot_cv(x$cv_results, x$h_cv, ...)

.plot_cv <- function(res, h_cv, ...) {
  ok <- is.finite(res$cvloss)
  if (!any(ok)) stop("no candidate bandwidth produced a loss to plot", call. = FALSE)
  plot(res$h[ok], res$cvloss[ok],
    type = "b", pch = 19, log = "x",
    xlab = "bandwidth h (log scale)", ylab = "held-out loss",
    main = "Cross-validation", ...
  )
  abline(v = h_cv, col = "#1f5aa6", lty = 2)
  points(h_cv, res$cvloss[which.min(abs(res$h - h_cv))],
    pch = 19, cex = 1.4, col = "#1f5aa6"
  )
  invisible(res)
}

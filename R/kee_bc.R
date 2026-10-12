#' Bias-Corrected Kernel Estimating Equations by Extrapolation in the Bandwidth
#'
#' @description
#' Applies the bandwidth selection and bias correction of [skmle_bc()] to the
#' kernel estimating-equation estimators of Cao and coauthors: [kee_cox()] for
#' the proportional hazards model with sparse longitudinal covariates, and
#' [kee_async()] for regression of an asynchronous longitudinal response.
#' Each function fits its estimator over a grid of bandwidths, chooses a
#' bandwidth by a plug-in estimate of the mean squared error of the
#' coefficients, and extrapolates the coefficients to bandwidth zero. The
#' standard error accounts for the correlation between fits at different
#' bandwidths and for the data-driven choice of bandwidth. See [skmle_bc()]
#' for details of the method.
#'
#' @param formula,data,id,obs_times As in [kee_cox()].
#' @param data_y,data_x,time As in [kee_async()].
#' @param h_grid Optional vector of candidate bandwidths. If `NULL`, `n_h`
#'   log-spaced values are used. `kee_cox_bc()` uses the default grid of
#'   [skmle_bc()]. `kee_async_bc()` uses values from
#'   \eqn{2\,\mathrm{IQR}\,n^{-0.6}} to the range of the pooled observation
#'   times, where IQR is the interquartile range of those times.
#' @param n_h Number of candidate bandwidths when `h_grid` is `NULL`.
#' @param degree Degree in \eqn{h} of the bias model and of the extrapolation,
#'   as in [skmle_bc()]. The default `1` removes the first-order bias of a
#'   one-sided kernel, and of a two-sided kernel when the covariate path is not
#'   smooth. When the bias of a two-sided kernel is of order \eqn{h^2}, linear
#'   extrapolation leaves a remainder of order \eqn{\hat h^2}, which is small
#'   relative to the standard error at the selected bandwidth. `2` removes
#'   this remainder but gives a larger variance.
#' @param n_boot,seed,quiet As in [skmle_bc()].
#' @param ... Further arguments passed to [kee_cox()] or [kee_async()], such as
#'   `one_sided`, `weight` or `link`.
#'
#' @return An object of class `c("kee_bc", "bc")` with the components
#'   described in [skmle_bc()]; `fit` is the [kee_cox()] or [kee_async()] fit
#'   at the selected bandwidth.
#'
#' @seealso [skmle_bc()], [kee_cox()], [kee_async()], [kee_async_cv()]
#' @name kee_bc
NULL

#' @rdname kee_bc
#' @export
kee_cox_bc <- function(formula, data, id, obs_times, h_grid = NULL, n_h = 12,
                       degree = 1, n_boot = 2000, seed = NULL, quiet = FALSE,
                       ...) {
  if (missing(formula) || missing(data) || missing(id) || missing(obs_times)) {
    stop("formula, data, id and obs_times must all be supplied")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  bc_check_n_boot(n_boot)
  degree <- bc_check_degree(degree)
  call <- match.call()
  env <- parent.frame()

  m <- match(c("formula", "data", "id", "obs_times"), names(call), 0L)
  mf_call <- call[c(1L, m)]
  mf_call[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_call, env)
  Y <- stats::model.response(mf)
  if (!inherits(Y, "Surv")) {
    stop("Response must be a survival object created with Surv()")
  }
  if (is.null(h_grid)) {
    h_grid <- bc_grid_surv(Y[, 1], stats::model.extract(mf, "obs_times"),
                           stats::model.extract(mf, "id"), n_h)
  }

  fit_call <- call
  fit_call[[1]] <- quote(skmle::kee_cox)
  fit_call$data <- data
  fit_call[c("h_grid", "n_h", "degree", "n_boot", "seed", "quiet")] <- NULL

  res <- bc_fit_grid(fit_call, env, h_grid, n_boot, seed, quiet, degree)
  res$fit$call <- fit_call
  res$fit$call$data <- call$data
  res$fit$call$h <- res$est$h

  out <- c(res$est, list(fit = res$fit, call = call))
  class(out) <- c("kee_bc", "bc")
  out
}

#' @rdname kee_bc
#' @export
kee_async_bc <- function(data_y, data_x, formula, id, time, h_grid = NULL,
                         n_h = 12, degree = 1, n_boot = 2000, seed = NULL,
                         quiet = FALSE, ...) {
  bc_check_n_boot(n_boot)
  degree <- bc_check_degree(degree)
  call <- match.call()
  env <- parent.frame()

  if (is.null(h_grid)) {
    d <- .async_parse(
      formula, data_y, data_x,
      .async_name(rlang::enquo(id)), .async_name(rlang::enquo(time)), TRUE
    )
    pooled <- c(d$yt, d$xt)
    iqr <- diff(stats::quantile(pooled, c(0.25, 0.75), names = FALSE))
    span <- diff(range(pooled))
    if (!is.finite(iqr) || iqr <= 0) iqr <- span
    if (!is.finite(span) || span <= 0) {
      stop("observation times are constant; supply 'h_grid'", call. = FALSE)
    }
    h_grid <- exp(seq(log(2 * iqr * d$n^(-0.6)), log(span),
                      length.out = n_h))
  }

  fit_call <- call
  fit_call[[1]] <- quote(skmle::kee_async)
  fit_call$data_y <- data_y
  fit_call$data_x <- data_x
  fit_call[c("h_grid", "n_h", "degree", "n_boot", "seed", "quiet")] <- NULL

  res <- bc_fit_grid(fit_call, env, h_grid, n_boot, seed, quiet, degree)
  res$fit$call <- fit_call
  res$fit$call$data_y <- call$data_y
  res$fit$call$data_x <- call$data_x
  res$fit$call$h <- res$est$h

  out <- c(res$est, list(fit = res$fit, call = call))
  class(out) <- c("kee_bc", "bc")
  out
}

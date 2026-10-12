# Why this exists ----------------------------------------------------------
#
# The SMKLE at bandwidth h pairs the hazard at t with a covariate observed at
# lag u in (0, h], so it estimates a parameter attenuated by that lag. For a
# covariate that jumps, the bias is linear in h whatever the kernel's shape.
# skmle_cv() cannot correct this: its held-out score carries covariates
# forward, so it rewards fits attenuated by the visit process's lag, and its
# selection converges to a constant rather than to zero (dev/sim_cv/theory).
#
# skmle_bc() works from the curve h -> beta_h instead. Its mean is
# beta0 + c h + e / h, so a weighted regression of beta_h on (1, h, 1/h) across
# a grid estimates the bias without beta0, and the intercept of a regression on
# (1, h) near the selected bandwidth extrapolates the estimate to h = 0. Both
# are linear combinations of the grid estimates, so their variance needs the
# covariance of estimates at two bandwidths, which the stacked estimating
# equations give from the per-subject scores every fit already computes.

#' Bias-Corrected SMKLE by Extrapolation in the Bandwidth
#'
#' @description
#' Fits [skmle()] over a grid of bandwidths, chooses a bandwidth by a plug-in
#' estimate of the mean squared error of the coefficients, and reports the
#' coefficients extrapolated to bandwidth zero. The standard error accounts
#' for the correlation between fits at different bandwidths.
#'
#' @inheritParams skmle
#' @param h_grid Optional increasing vector of candidate bandwidths. If `NULL`,
#'   `n_h` values log-spaced from
#'   \eqn{\max\{\min(X_i - T_{ij} : X_i > T_{ij}),\ \tau n^{-0.6}\}} to
#'   \eqn{\tau = \max_i X_i}.
#' @param n_h Number of candidate bandwidths when `h_grid` is `NULL`.
#' @param degree Degree in \eqn{h} of the bias model and of the
#'   extrapolation: `1` (the default) removes the first-order bias of a
#'   one-sided kernel, and `2` also removes a term in \eqn{h^2}.
#' @param n_boot Number of simulated grids used to include in `var` the
#'   variability from the data-driven choice of bandwidth. With `0`, `var` is
#'   the variance that treats the selected bandwidth and window as fixed.
#' @param seed Optional integer seed for those simulated grids. If `NULL`, the
#'   current RNG state is used.
#' @param quiet Logical; if `TRUE`, suppress progress output.
#' @param ... Further arguments passed to [skmle()], such as `nknots`,
#'   `maxeval` or `optimizer`.
#'
#' @details
#' # The bias
#'
#' At bandwidth \eqn{h} the estimator evaluates the hazard at each time with
#' covariates observed up to \eqn{h} earlier. The lag acts as measurement
#' error and attenuates the coefficients. For a covariate that changes by
#' jumps, the resulting bias is linear in \eqn{h}, and a small-sample term of
#' order \eqn{1/(nh)} is added to it. The mean of \eqn{\hat\beta_h} is
#' therefore modelled as \eqn{\beta_0 + c h + e/h}.
#'
#' # The selected bandwidth
#'
#' For each coefficient, \eqn{c} and \eqn{e} are estimated by a regression of
#' \eqn{\hat\beta_h} on \eqn{(1, h, 1/h)} across the grid, weighted by the
#' inverse sandwich variances. The selected bandwidth `h` is the grid point
#' that minimises the estimated mean squared error summed over coefficients,
#' \eqn{\sum_k \{(\hat c_k h + \hat e_k / h)^2 + \widehat{\mathrm{var}}(\hat\beta_{k,h})\}}.
#' `fit` is the [skmle()] fit at that bandwidth.
#'
#' # The extrapolated estimate
#'
#' `coefficients` is the intercept of a weighted regression of
#' \eqn{\hat\beta_h} on \eqn{(1, h)} over the grid points in
#' \eqn{[h/2, 2h]}. This removes the linear term of the bias. The intercept is
#' a linear combination \eqn{\sum_j a_j \hat\beta_{h_j}} of the grid
#' estimates, so with the bandwidth and window held fixed its covariance
#' `var_fixed` follows from the joint sandwich
#' \deqn{\mathrm{Cov}(\hat\beta_{h_j}, \hat\beta_{h_m}) =
#'   A_{h_j}^{-1} \Bigl(\tfrac{1}{n}\textstyle\sum_i S_{i,h_j} S_{i,h_m}^\top\Bigr)
#'   A_{h_m}^{-1} / n,}
#' where \eqn{S_{i,h}} is subject \eqn{i}'s score contribution at bandwidth
#' \eqn{h}.
#'
#' Because the bandwidth and window are chosen from the same estimates, the
#' extrapolated estimate has additional variance that `var_fixed` omits.
#' `var` includes it. `n_boot` grids are drawn from a normal distribution
#' whose mean is the grid estimates and whose covariance is the joint
#' sandwich; the selection and extrapolation are rerun on each simulated grid,
#' and `var` is the covariance of the resulting estimates. No model is
#' refitted, so this step takes a few seconds. In simulation at \eqn{s = 1},
#' the data-driven selection added 6--9\% to the standard deviation of the
#' estimate.
#'
#' The extrapolated estimate has smaller bias and larger variance than the
#' coefficients of `fit`. Use it with `var` for confidence intervals.
#'
#' A selected bandwidth at an end of the grid raises a warning. Grid points
#' whose fit has no usable variance are dropped, also with a warning.
#'
#' @return
#' An object of class `c("skmle_bc", "bc")` with components:
#'
#' * `coefficients`, `var`: the extrapolated estimate and its covariance,
#'   including the variance from choosing the bandwidth,
#' * `var_fixed`: the covariance treating the selected bandwidth and window as
#'   fixed,
#' * `h`: the selected bandwidth,
#' * `fit`: the [skmle()] fit at `h`,
#' * `window`: the bandwidths used in the extrapolation,
#' * `weights`: the matrix of extrapolation weights \eqn{a_j}, one row per grid
#'   point and one column per coefficient,
#' * `grid`: a data frame with one row per bandwidth and coefficient giving the
#'   estimate, its standard error, the estimated bias and the estimated mean
#'   squared error,
#' * `cov_all`: the joint covariance of the estimates over the grid, ordered by
#'   bandwidth and then coefficient,
#' * `n_boot`: the number of simulated grids behind `var`,
#' * `degree`: the degree of the bias model and extrapolation,
#' * `call`: the matched call.
#'
#' @examples
#' \donttest{
#' library(survival)
#' set.seed(123)
#' dat <- sim_skmle_data(
#'   n = 150,
#'   mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'   mu_bar = 8,
#'   alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'   beta = c(1, -0.5),
#'   s = 0,
#'   cen = 0.7
#' )
#' bc <- skmle_bc(Surv(X, delta) ~ covariates, data = dat, id = id,
#'                obs_times = obs_times, s = 0, n_h = 8, quiet = TRUE)
#' bc
#' summary(bc)
#' }
#'
#' @seealso [skmle()], [skmle_cv()]
#' @export
skmle_bc <- function(formula, data, id, obs_times, s = 0, h_grid = NULL,
                     n_h = 12, degree = 1, n_boot = 2000, seed = NULL,
                     quiet = FALSE, ...) {
  if (missing(formula) || missing(data) || missing(id) || missing(obs_times)) {
    stop("formula, data, id and obs_times must all be supplied")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  bc_check_n_boot(n_boot)
  degree <- bc_check_degree(degree)
  call <- match.call()
  env <- parent.frame()

  # Times and subjects, for the default grid.
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

  # Freeze `data` into the call so every fit sees the same data frame, as
  # skmle_cv() does for its refit.
  fit_call <- call
  fit_call[[1]] <- quote(skmle::skmle)
  fit_call$data <- data
  fit_call[c("h_grid", "n_h", "degree", "n_boot", "seed", "quiet")] <- NULL

  res <- bc_fit_grid(fit_call, env, h_grid, n_boot, seed, quiet, degree)
  res$fit$call <- fit_call
  res$fit$call$data <- call$data
  res$fit$call$h <- res$est$h

  out <- c(res$est, list(fit = res$fit, call = call))
  class(out) <- c("skmle_bc", "bc")
  out
}

# Default grid for the survival estimators: log-spaced from the larger of the
# smallest positive lag and tau n^-0.6 up to tau = max(X).
bc_grid_surv <- function(X_time, obs_times, id, n_h) {
  obs_times <- as.numeric(obs_times)
  n <- length(unique(id))
  tau <- max(X_time)
  lags <- X_time - obs_times
  lags <- lags[lags > 0]
  if (length(lags) == 0) {
    stop("no observation time precedes its event or censoring time; ",
         "cannot form a default 'h_grid'", call. = FALSE)
  }
  hmin <- max(min(lags), tau * n^(-0.6))
  exp(seq(log(hmin), log(tau), length.out = n_h))
}

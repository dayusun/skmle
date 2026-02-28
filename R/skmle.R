#' Fit skmle model
#'
#' @description Fits a transformed hazards model for survival data with sparsely
#' and intermittently observed longitudinal covariates using the sieve maximum
#' kernel-weighted log-likelihood estimator (SMKLE).
#'
#' @param formula A formula object, with the response on the left of a `~` operator, and the terms on the right. The response must be a survival object as returned by the `Surv` function.
#' @param data A data.frame in which to interpret the variables named in the formula, or in the `id` argument.
#' @param id A vector identifying subjects in the data.
#' @param obs_times The observation times for the subjects.
#' @param s The s parameter for the Box-Cox transformation function. `s = 0` corresponds to the proportional hazards model, and `s = 1` corresponds to the additive hazards model.
#' @param h The bandwidth parameter for kernel smoothing.
#' @param nknots The number of knots for the B-splines used to approximate the baseline hazard function. Default is 3.
#' @param norder The order of the B-splines. Default is 3 (cubic splines).
#' @param lq_nodes The number of Legendre-Gauss quadrature nodes for the integral of the baseline hazard. Default is 64.
#' @param maxeval Maximum number of evaluations for nlopt optimization. Default is 10000.
#' @param xtol_rel Relative tolerance for nlopt optimization. Default is 1e-6.
#'
#' @details
#' The `skmle` function implements the Sieve Maximum Kernel-weighted Log-likelihood
#' Estimator proposed by Sun et al. (2025). The method combines kernel weighting
#' on the log-likelihood (to handle sparse longitudinal covariates evaluated at event
#' or censoring times) with a B-spline sieve approximation for the unknown
#' nonparametric baseline hazard function. The estimation is optimized using an
#' efficient C++ backend via `Rcpp` and `nloptr`.
#'
#' @references
#' Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao.
#' "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
#' *Journal of the American Statistical Association* (2025): 1-12.
#'
#' @return A list of class `skmle` containing the estimation results, including
#' estimated coefficients, variance-covariance matrix, log-likelihood, and optimization details.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' # Suppose dat is a dataframe with sparse longitudinal data
#' # id: subject ID
#' # obs_times: times at which covariates were measured
#' # time: survival/censoring time
#' # status: event indicator (1 = event, 0 = censored)
#' # x1, x2: time-dependent covariates measured at obs_times
#'
#' # Fit the SMKLE model with Box-Cox parameter s = 0 (proportional hazards)
#' fit <- skmle(Surv(time, status) ~ x1 + x2,
#'   data = dat, id = id, obs_times = obs_times,
#'   s = 0, h = 0.5, nknots = 3
#' )
#' summary(fit)
#' plot(fit)
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom splines ns
#' @importFrom gaussquad legendre.quadrature.rules
#' @export
skmle <- function(formula, data, id, obs_times, s, h, nknots = 3, norder = 3, lq_nodes = 64, maxeval = 10000, xtol_rel = 1e-6) {
  mf <- model.frame(formula, data)
  Y <- model.response(mf)

  if (!inherits(Y, "Surv")) {
    stop("Response must be a survival object created with Surv()")
  }

  X_time <- Y[, 1]
  delta <- Y[, 2]

  Z <- model.matrix(formula, data)

  # Remove intercept if present
  if (colnames(Z)[1] == "(Intercept)") {
    Z <- Z[, -1, drop = FALSE]
  }

  id_vec <- as.numeric(data[[deparse(substitute(id))]])
  if (is.null(id_vec)) {
    id_vec <- as.numeric(eval(substitute(id), data, parent.frame()))
  }

  obs_times_vec <- as.numeric(data[[deparse(substitute(obs_times))]])
  if (is.null(obs_times_vec)) {
    obs_times_vec <- as.numeric(eval(substitute(obs_times), data, parent.frame()))
  }

  kerfun <- function(xx) {
    pmax((1 - xx^2) * 0.75, 0)
  }

  n <- length(unique(id_vec))

  kerval <- kerfun((X_time - obs_times_vec) / h) / h * as.numeric(X_time > obs_times_vec)

  knots <- (1:nknots) / (nknots + 1)

  bsmat <- splines::ns(X_time, knots = knots, intercept = TRUE, Boundary.knots = c(0, 1))

  lqrule <- gaussquad::legendre.quadrature.rules(lq_nodes)[[lq_nodes]]

  lq_x <- lqrule$x
  lq_w <- lqrule$w

  # Precompute for quadrature
  n_quad <- length(lq_x)
  bsmat_tt_all <- array(0, dim = c(n_quad, ncol(bsmat), 1))
  kerval_tt_all <- matrix(0, nrow = n_quad, ncol = length(X_time))

  for (q in 1:n_quad) {
    tt <- 0.5 * lq_x[q] + 0.5

    dist_tt <- tt - obs_times_vec
    kerval_tt_all[q, ] <- kerfun(dist_tt / h) / h * as.numeric(dist_tt > 0)

    bsmat_tt <- splines::ns(rep(tt, 2), knots = knots, intercept = TRUE, Boundary.knots = c(0, 1))[1, , drop = FALSE]
    bsmat_tt_all[q, , 1] <- bsmat_tt
  }

  # Reorganize array for C++
  bsmat_tt_mat <- matrix(0, nrow = n_quad, ncol = ncol(bsmat))
  for (q in 1:n_quad) {
    bsmat_tt_mat[q, ] <- bsmat_tt_all[q, , 1]
  }


  ineqmat <- matrix(0, nrow = length(X_time), ncol = ncol(Z) + ncol(bsmat))
  if (s != 0) {
    ineqmat <- cbind(Z, as.matrix(bsmat))
    ineqmat <- ineqmat[(X_time > obs_times_vec) & (abs(X_time - obs_times_vec) <= h), ]
  }

  res <- skmle_cpp_fit(
    n = n,
    p = ncol(Z),
    gammap = ncol(bsmat),
    s = s,
    h = h,
    covariates = as.matrix(Z),
    bsmat = as.matrix(bsmat),
    X = as.numeric(X_time),
    obs_times = as.numeric(obs_times_vec),
    delta = as.numeric(delta),
    kerval = as.numeric(kerval),
    lq_x = as.numeric(lq_x),
    lq_w = as.numeric(lq_w),
    bsmat_tt_all = bsmat_tt_mat,
    kerval_tt_all = as.matrix(kerval_tt_all),
    ineqmat = as.matrix(ineqmat),
    maxeval = maxeval,
    xtol_rel = xtol_rel
  )

  beta_est <- res$solution[1:ncol(Z)]
  gamma_est <- res$solution[(ncol(Z) + 1):(ncol(Z) + ncol(bsmat))]

  A_est <- calc_A(
    beta = beta_est,
    gamma = gamma_est,
    s = s,
    h = h,
    covariates = as.matrix(Z),
    bsmat = as.matrix(bsmat),
    X = as.numeric(X_time),
    obs_times = as.numeric(obs_times_vec),
    delta = as.numeric(delta),
    kerval = as.numeric(kerval),
    bsmat_XX = as.matrix(bsmat),
    n_subj = n
  )

  B_est <- calc_B(
    beta = beta_est,
    gamma = gamma_est,
    s = s,
    h = h,
    covariates = as.matrix(Z),
    bsmat = as.matrix(bsmat),
    X = as.numeric(X_time),
    obs_times = as.numeric(obs_times_vec),
    delta = as.numeric(delta),
    kerval = as.numeric(kerval),
    id = as.numeric(id_vec),
    bsmat_XX = as.matrix(bsmat),
    lq_x = as.numeric(lq_x),
    lq_w = as.numeric(lq_w),
    bsmat_tt_all = bsmat_tt_mat,
    kerval_tt_all = as.matrix(kerval_tt_all),
    n_subj = n
  )

  var_est <- solve(A_est) %*% B_est %*% solve(A_est) / n

  names(beta_est) <- colnames(Z)
  dimnames(var_est) <- list(colnames(Z), colnames(Z))

  out <- list(
    coefficients = beta_est,
    var = var_est,
    gamma = gamma_est,
    loglik = -res$minimum,
    A = A_est,
    B = B_est,
    convergence = res$status,
    n = n,
    s = s,
    h = h,
    nknots = nknots,
    call = match.call()
  )

  class(out) <- "skmle"
  return(out)
}

#' Print skmle object
#'
#' @param x An object of class `skmle`.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.skmle <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' Summary for skmle object
#'
#' @param object An object of class `skmle`.
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats pnorm
#' @export
summary.skmle <- function(object, ...) {
  coef <- object$coefficients
  se <- sqrt(diag(object$var))
  z_val <- coef / se
  p_val <- 2 * pnorm(-abs(z_val))

  coef_table <- cbind(
    Estimate = coef,
    `Std. Error` = se,
    `z value` = z_val,
    `Pr(>|z|)` = p_val
  )

  res <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = object$loglik,
    convergence = object$convergence,
    n = object$n
  )
  class(res) <- "summary.skmle"
  return(res)
}

#' Print summary of skmle object
#'
#' @param x An object of class `summary.skmle`.
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats printCoefmat
#' @export
print.summary.skmle <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat(sprintf("  n= %d\n\n", x$n))

  if (nrow(x$coefficients) > 0) {
    stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  }

  cat("\nLog-likelihood:", format(x$loglik, digits = 4), "\n")
  if (x$convergence < 0) {
    cat("Optimization may not have converged (status:", x$convergence, ")\n")
  }
  invisible(x)
}

#' Plot the estimated baseline function for skmle model
#'
#' @param x An object of class `skmle`.
#' @param t_seq A numeric vector of time points to evaluate the baseline function. Default is `seq(0, 1, length.out = 100)`.
#' @param ... Further arguments passed to or from other methods.
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#' @importFrom splines ns
#' @export
plot.skmle <- function(x, t_seq = seq(0, 1, length.out = 100), ...) {
  if (!inherits(x, "skmle")) {
    stop("Object must be of class 'skmle'")
  }

  if (!is.null(x$nknots)) {
    nknots <- x$nknots
  } else {
    nknots <- length(x$gamma) - 2
  }

  knots <- (1:nknots) / (nknots + 1)

  bsmat_plot <- splines::ns(t_seq, knots = knots, intercept = TRUE, Boundary.knots = c(0, 1))

  baseline_est <- as.vector(bsmat_plot %*% x$gamma)

  plot_data <- data.frame(
    Time = t_seq,
    Baseline = baseline_est
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Baseline)) +
    ggplot2::geom_line(color = "blue", linewidth = 1) +
    ggplot2::labs(
      title = "Estimated Baseline Function",
      x = "Time",
      y = "Baseline Value"
    ) +
    ggplot2::theme_minimal()

  return(p)
}

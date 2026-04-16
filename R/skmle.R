#' Fit a Transformed Hazards Model by SMKLE
#'
#' @description
#' Fit a transformed hazards model for survival data with sparsely and intermittently
#' observed longitudinal covariates using the sieve maximum kernel-weighted
#' log-likelihood estimator (SMKLE).
#'
#' @param formula A model formula. The left-hand side must be a `survival::Surv()`
#'   response and the right-hand side must contain at least one covariate.
#' @param data Data frame containing the variables used in `formula`, `id`, and
#'   `obs_times`.
#' @param id Subject identifier. Non-numeric identifiers are allowed and are internally
#'   converted to integer subject codes.
#' @param obs_times Longitudinal observation times aligned row-wise with `data`.
#' @param s Box-Cox transformation parameter. `s = 0` gives the proportional hazards
#'   model and `s = 1` gives the additive hazards model.
#' @param h Positive kernel bandwidth.
#' @param nknots Number of interior knots used in the sieve approximation of the
#'   baseline component.
#' @param norder Order parameter supplied to the high-level interface for the spline
#'   approximation.
#' @param lq_nodes Number of Legendre-Gauss quadrature nodes used in numerical
#'   integration.
#' @param maxeval Maximum number of optimizer evaluations.
#' @param xtol_rel Relative convergence tolerance passed to the optimizer.
#'
#' @details
#' `skmle()` is the main model-fitting function in the package. It combines:
#'
#' * kernel weighting to handle intermittently observed longitudinal covariates,
#' * a sieve approximation for the unknown baseline component, and
#' * a C++-backed numerical optimizer for the joint estimation problem.
#'
#' The returned object follows the usual R model pattern: print the fitted coefficients
#' with `print()`, obtain inferential output with `summary()`, and visualize the
#' estimated baseline component with `plot()`.
#'
#' @references
#' Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao.
#' "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
#' *Journal of the American Statistical Association* (2025): 1-12.
#'
#' @return
#' An object of class `skmle` containing:
#'
#' * `coefficients`: regression coefficient estimates,
#' * `var`: estimated variance-covariance matrix for the regression coefficients,
#' * `gamma`: estimated spline coefficients for the baseline component,
#' * `loglik`: maximized log-likelihood,
#' * `convergence`: optimizer status code,
#' * model metadata such as `n`, `s`, `h`, and `call`.
#'
#' @examples
#' library(survival)
#'
#' set.seed(123)
#' dat <- sim_skmle_data(
#'   n = 80,
#'   mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'   mu_bar = 8,
#'   alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'   beta = c(1, -0.5),
#'   s = 0,
#'   cen = 0.7
#' )
#'
#' fit <- skmle(
#'   Surv(X, delta) ~ covariates,
#'   data = dat,
#'   id = id,
#'   obs_times = obs_times,
#'   s = 0,
#'   h = 0.5,
#'   nknots = 3
#' )
#'
#' fit
#' summary(fit)
#'
#' # If you prefer explicit covariate names, split the matrix column first.
#' dat$Z1 <- dat$covariates[, 1]
#' dat$Z2 <- dat$covariates[, 2]
#' fit_named <- skmle(
#'   Surv(X, delta) ~ Z1 + Z2,
#'   data = dat,
#'   id = id,
#'   obs_times = obs_times,
#'   s = 0,
#'   h = 0.5,
#'   nknots = 3
#' )
#'
#' summary(fit_named)
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom splines ns
#' @importFrom gaussquad legendre.quadrature.rules
#' @export
skmle <- function(formula, data, id, obs_times, s, h, nknots = 3, norder = 3, lq_nodes = 64, maxeval = 10000, xtol_rel = 1e-6) {
  # validate inputs --------------------------------------------------------
  if (missing(formula) || missing(data) || missing(id) ||
    missing(obs_times) || missing(s) || missing(h)) {
    stop("formula, data, id, obs_times, s and h must all be supplied")
  }
  if (!is.numeric(h) || length(h) != 1 || h <= 0) stop("'h' must be a positive number")
  if (!is.numeric(s) || length(s) != 1) stop("'s' must be a numeric scalar")
  if (!is.numeric(nknots) || length(nknots) != 1 || nknots < 1) stop("'nknots' must be >= 1")
  if (!is.numeric(norder) || length(norder) != 1 || norder < 1) stop("'norder' must be >= 1")
  if (!is.numeric(lq_nodes) || length(lq_nodes) != 1 || lq_nodes < 1) stop("'lq_nodes' must be >= 1")
  if (!is.data.frame(data)) stop("'data' must be a data.frame")

  # robust parsing of formula, id, and obs_times to synchronize NA dropping
  call <- match.call()
  m <- match(c("formula", "data", "id", "obs_times"), names(call), 0L)
  mf_call <- call[c(1L, m)]
  mf_call[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_call, parent.frame())

  Y <- model.response(mf)

  if (!inherits(Y, "Surv")) {
    stop("Response must be a survival object created with Surv()")
  }

  X_time <- Y[, 1]
  delta <- Y[, 2]
  if (anyNA(X_time) || anyNA(delta)) {
    stop("missing values in survival response not permitted")
  }

  Z <- model.matrix(formula, data = mf)
  if (ncol(Z) > 0 && colnames(Z)[1] == "(Intercept)") {
    Z <- Z[, -1, drop = FALSE]
  }
  if (ncol(Z) == 0) stop("model must contain at least one covariate")

  id_raw <- model.extract(mf, "id")
  obs_times_vec <- as.numeric(model.extract(mf, "obs_times"))

  if (length(id_raw) != length(X_time) || length(obs_times_vec) != length(X_time)) {
    stop("Length of 'id' and 'obs_times' must match number of rows in data/formula")
  }
  id_vec <- as.integer(factor(id_raw))

  kerfun <- function(xx) {
    res <- (1 - xx^2) * 0.75
    res[res < 0] <- 0
    res
  }

  n <- length(unique(id_vec))
  kerval <- kerfun((X_time - obs_times_vec) / h) / h * as.numeric(X_time > obs_times_vec)

  knots <- (1:nknots) / (nknots + 1)
  bsmat <- splines::ns(X_time, knots = knots, intercept = TRUE, Boundary.knots = c(0, 1))

  lqrule <- gaussquad::legendre.quadrature.rules(lq_nodes)[[lq_nodes]]
  lq_x <- lqrule$x
  lq_w <- lqrule$w

  # Precompute for quadrature points (matrix form only)
  n_quad <- length(lq_x)
  tts <- 0.5 * lq_x + 0.5

  dist_tt_mat <- outer(tts, obs_times_vec, "-")
  kerval_tt_all <- kerfun(dist_tt_mat / h) / h * as.numeric(dist_tt_mat > 0)
  bsmat_tt_mat <- as.matrix(splines::ns(tts, knots = knots, intercept = TRUE, Boundary.knots = c(0, 1)))

  # inequality constraints matrix, may be empty if no rows satisfy the filter
  ineqmat <- matrix(numeric(0), nrow = 0, ncol = ncol(Z) + ncol(bsmat))
  if (s != 0) {
    tmp <- cbind(Z, as.matrix(bsmat))
    rows <- (X_time > obs_times_vec) & (abs(X_time - obs_times_vec) <= h)
    if (any(rows)) ineqmat <- tmp[rows, , drop = FALSE]
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

  var_est <- tryCatch(
    solve(A_est) %*% B_est %*% solve(A_est) / n,
    error = function(e) {
      warning("variance estimation failed, A_est may be singular: ", e$message)
      matrix(NA_real_, ncol(A_est), ncol(A_est))
    }
  )

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
#' @return `x`, invisibly. Called for its side effect of printing a brief
#'   summary of the call and estimated coefficients.
#' @examples
#' \dontrun{
#' fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
#'              obs_times = obs_times, s = 0, h = 0.5)
#' print(fit)
#' }
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
#' @return An object of class `summary.skmle` containing the call,
#'   a coefficient table with estimates, standard errors, z-statistics
#'   and p-values, the log-likelihood, convergence status, and the sample
#'   size `n`.
#' @examples
#' \dontrun{
#' fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
#'              obs_times = obs_times, s = 0, h = 0.5)
#' summary(fit)
#' }
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
#' @return `x`, invisibly. Called for its side effect of printing the
#'   summary table.
#' @examples
#' \dontrun{
#' fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
#'              obs_times = obs_times, s = 0, h = 0.5)
#' print(summary(fit))
#' }
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
#' @return A `ggplot` object showing the estimated nonparametric baseline
#'   function evaluated on `t_seq`.
#' @examples
#' \dontrun{
#' fit <- skmle(Surv(X, delta) ~ covariates, data = dat, id = id,
#'              obs_times = obs_times, s = 0, h = 0.5)
#' plot(fit)
#' }
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#' @importFrom splines ns
#' @export
plot.skmle <- function(x, t_seq = seq(0, 1, length.out = 100), ...) {
  if (!inherits(x, "skmle")) {
    stop("Object must be of class 'skmle'")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting; please install it")
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

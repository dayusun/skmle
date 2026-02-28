#' Fit the Cox model using estimating equations
#'
#' @description Fits the proportional hazards (Cox) model for survival data with
#' sparse longitudinal covariates using a kernel-weighted estimating equations approach.
#'
#' @param formula A formula object, with the response on the left of a `~` operator, and the terms on the right. The response must be a survival object as returned by the `Surv` function.
#' @param data A data.frame in which to interpret the variables named in the formula, or in the `id` argument.
#' @param id A vector identifying subjects in the data.
#' @param obs_times The observation times for the subjects.
#' @param h The bandwidth parameter for kernel smoothing.
#'
#' @details
#' `kee_cox` uses a kernel-smoothed partial likelihood estimating equation to estimate
#' regression coefficients when covariates are measured intermittently and sparsely.
#' It circumvents the estimation of the unknown nonparametric baseline hazard function.
#'
#' @references
#' Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao.
#' "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
#' *Journal of the American Statistical Association* (2025): 1-12.
#'
#' Cao, Hongyuan, et al. "Inference for Cox models with sparse longitudinal covariates."
#' *Biometrika* (2015).
#'
#' @return A list of class `kee` containing the estimation results, including
#' estimated coefficients, variance-covariance matrix, and optimization details.
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' # Simulate data for 200 subjects under the proportional hazards model (s = 0)
#' set.seed(123)
#' dat <- sim_skmle_data(
#'     n = 200,
#'     mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'     mu_bar = 8,
#'     alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'     beta = c(1, -0.5), # True coefficients
#'     s = 0, # proportional hazards
#'     cen = 0.7 # censoring parameter
#' )
#'
#' # Fit the Cox model using kernel estimating equations
#' fit_cox <- kee_cox(Surv(X, delta) ~ covariates,
#'     data = dat, id = id, obs_times = obs_times,
#'     h = 0.5
#' )
#' summary(fit_cox)
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom nleqslv nleqslv
#' @export
kee_cox <- function(formula, data, id, obs_times, h) {
    mf <- model.frame(formula, data)
    Y <- model.response(mf)

    if (!inherits(Y, "Surv")) {
        stop("Response must be a survival object created with Surv()")
    }

    X_time <- Y[, 1]
    delta <- Y[, 2]

    Z <- model.matrix(formula, data)
    if (colnames(Z)[1] == "(Intercept)") Z <- Z[, -1, drop = FALSE]

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

    estequ <- function(beta) {
        kee_cox_estequ(
            beta = beta,
            covariates = as.matrix(Z),
            X = as.numeric(X_time),
            obs_times = as.numeric(obs_times_vec),
            delta = as.numeric(delta),
            kerval = as.numeric(kerval),
            h = h,
            n_subj = n
        )
    }

    init_beta <- rep(0, ncol(Z))
    estres <- nleqslv::nleqslv(init_beta, fn = estequ)
    beta_est <- estres$x

    var_res <- kee_cox_var(
        beta = beta_est,
        covariates = as.matrix(Z),
        X = as.numeric(X_time),
        obs_times = as.numeric(obs_times_vec),
        delta = as.numeric(delta),
        kerval = as.numeric(kerval),
        h = h,
        id = as.numeric(id_vec),
        n_subj = n
    )

    W <- var_res$W
    Sigma <- var_res$Sigma
    var_est <- solve(W) %*% Sigma %*% solve(W)

    names(beta_est) <- colnames(Z)
    dimnames(var_est) <- list(colnames(Z), colnames(Z))

    out <- list(
        coefficients = beta_est,
        var = var_est,
        W = W,
        Sigma = Sigma,
        convergence = estres$termcd,
        n = n,
        h = h,
        call = match.call()
    )
    class(out) <- "kee"
    return(out)
}

#' Fit the Additive model using estimating equations
#'
#' @description Fits the additive hazards model for survival data with
#' sparse longitudinal covariates using a kernel-weighted estimating equations approach.
#'
#' @param formula A formula object, with the response on the left of a `~` operator, and the terms on the right. The response must be a survival object as returned by the `Surv` function.
#' @param data A data.frame in which to interpret the variables named in the formula, or in the `id` argument.
#' @param id A vector identifying subjects in the data.
#' @param obs_times The observation times for the subjects.
#' @param h The bandwidth parameter for kernel smoothing.
#' @param lq_nodes The number of Legendre-Gauss quadrature nodes for the integral. Default is 64.
#'
#' @details
#' `kee_additive` uses a kernel-smoothed martingale-based estimating equation to estimate
#' regression coefficients in the additive hazards model when covariates are measured intermittently and sparsely.
#'
#' @references
#' Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao.
#' "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
#' *Journal of the American Statistical Association* (2025): 1-12.
#'
#' Sun, Dayu, Hongyuan Cao, and Yining Chen. "Additive hazards models with sparse
#' longitudinal covariates." *Lifetime Data Analysis* (2022).
#'
#' @return A list of class `kee` containing the estimation results, including
#' estimated coefficients, and variance-covariance matrix.
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' # Simulate data for 200 subjects under the additive hazards model (s = 1)
#' set.seed(123)
#' dat <- sim_skmle_data(
#'     n = 200,
#'     mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'     mu_bar = 8,
#'     alpha = function(tt) 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'     beta = c(1, -0.5), # True coefficients
#'     s = 1, # additive hazards
#'     cen = 0.7 # censoring parameter
#' )
#'
#' # Fit the additive hazards model using kernel estimating equations
#' fit_add <- kee_additive(Surv(X, delta) ~ covariates,
#'     data = dat, id = id, obs_times = obs_times,
#'     h = 0.5
#' )
#' summary(fit_add)
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom gaussquad legendre.quadrature.rules
#' @export
kee_additive <- function(formula, data, id, obs_times, h, lq_nodes = 64) {
    mf <- model.frame(formula, data)
    Y <- model.response(mf)

    if (!inherits(Y, "Surv")) {
        stop("Response must be a survival object created with Surv()")
    }

    X_time <- Y[, 1]
    delta <- Y[, 2]

    Z <- model.matrix(formula, data)
    if (colnames(Z)[1] == "(Intercept)") Z <- Z[, -1, drop = FALSE]

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

    lqrule <- gaussquad::legendre.quadrature.rules(lq_nodes)[[lq_nodes]]
    lq_x <- lqrule$x
    lq_w <- lqrule$w

    res <- kee_additive_est(
        covariates = as.matrix(Z),
        X = as.numeric(X_time),
        obs_times = as.numeric(obs_times_vec),
        delta = as.numeric(delta),
        kerval = as.numeric(kerval),
        h = h,
        id = as.numeric(id_vec),
        lq_x = as.numeric(lq_x),
        lq_w = as.numeric(lq_w),
        n_subj = n
    )

    beta_est <- as.vector(res$est)
    A <- res$A_est
    B <- res$B_est
    Sigma <- res$Sigma_est
    var_est <- solve(A) %*% Sigma %*% solve(A)

    names(beta_est) <- colnames(Z)
    dimnames(var_est) <- list(colnames(Z), colnames(Z))

    out <- list(
        coefficients = beta_est,
        var = var_est,
        A = A,
        B = B,
        Sigma = Sigma,
        n = n,
        h = h,
        call = match.call()
    )
    class(out) <- "kee"
    return(out)
}

#' Print kee object
#'
#' @param x An object of class `kee`.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.kee <- function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    invisible(x)
}

#' Summary for kee object
#'
#' @param object An object of class `kee`.
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats pnorm
#' @export
summary.kee <- function(object, ...) {
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
        convergence = if (is.null(object$convergence)) NA else object$convergence,
        n = object$n
    )
    class(res) <- "summary.kee"
    return(res)
}

#' Print summary of kee object
#'
#' @param x An object of class `summary.kee`.
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats printCoefmat
#' @export
print.summary.kee <- function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat(sprintf("  n= %d\n\n", x$n))

    if (nrow(x$coefficients) > 0) {
        stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    }

    if (!is.na(x$convergence) && x$convergence != 1) {
        cat("\nOptimization status:", x$convergence, "\n")
    }
    invisible(x)
}

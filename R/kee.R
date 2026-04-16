#' Fit a Cox-Type KEE Model
#'
#' @description
#' Fit the proportional hazards model for sparse longitudinal covariate data using
#' a kernel estimating-equation approach.
#'
#' @param formula A model formula with a `survival::Surv()` response.
#' @param data Data frame containing all variables used in the fit.
#' @param id Subject identifier aligned row-wise with `data`.
#' @param obs_times Longitudinal observation times aligned row-wise with `data`.
#' @param h Positive kernel bandwidth.
#'
#' @details
#' `kee_cox()` targets the proportional hazards case without estimating a nonparametric
#' baseline component. It is therefore a useful specialized alternative to `skmle()`
#' when the scientific model is Cox-type and the main interest is in the regression
#' coefficients.
#'
#' @references
#' Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao.
#' "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
#' *Journal of the American Statistical Association* (2025): 1-12.
#'
#' Cao, Hongyuan, et al. "Inference for Cox models with sparse longitudinal covariates."
#' *Biometrika* (2015).
#'
#' @return
#' An object of class `kee` containing coefficient estimates, the estimated
#' variance-covariance matrix, the estimating-equation matrices, convergence status,
#' and the original function call.
#'
#' @examples
#' \donttest{
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
#' fit_cox <- kee_cox(
#'   Surv(X, delta) ~ covariates,
#'   data = dat,
#'   id = id,
#'   obs_times = obs_times,
#'   h = 0.5
#' )
#'
#' summary(fit_cox)
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom nleqslv nleqslv
#' @export
kee_cox <- function(formula, data, id, obs_times, h) {
    # basic input validation
    if (missing(formula) || missing(data) || missing(id) || missing(obs_times) || missing(h)) {
        stop("formula, data, id, obs_times and h must all be supplied")
    }
    if (!is.numeric(h) || length(h) != 1 || h <= 0) {
        stop("bandwidth 'h' must be a positive number")
    }
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame")
    }

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
    # drop intercept if present
    if (ncol(Z) > 0 && colnames(Z)[1] == "(Intercept)") {
        Z <- Z[, -1, drop = FALSE]
    }
    if (ncol(Z) == 0) {
        stop("model must contain at least one covariate")
    }

    id_raw <- model.extract(mf, "id")
    obs_times_vec <- as.numeric(model.extract(mf, "obs_times"))

    if (length(id_raw) != length(X_time) || length(obs_times_vec) != length(X_time)) {
        stop("Length of 'id' and 'obs_times' must match number of rows in data/formula")
    }
    # convert identifier to integer codes for safe use in C++
    id_vec <- as.integer(factor(id_raw))

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
    var_est <- safe_sandwich(W, Sigma, scale = 1, what = "kee_cox variance")

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

#' Fit an Additive Hazards KEE Model
#'
#' @description
#' Fit the additive hazards model for sparse longitudinal covariate data using
#' a kernel estimating-equation approach.
#'
#' @param formula A model formula with a `survival::Surv()` response.
#' @param data Data frame containing all variables used in the fit.
#' @param id Subject identifier aligned row-wise with `data`.
#' @param obs_times Longitudinal observation times aligned row-wise with `data`.
#' @param h Positive kernel bandwidth.
#' @param lq_nodes Number of quadrature nodes used in the numerical integration step.
#'
#' @details
#' `kee_additive()` is the specialized additive-hazards counterpart to `kee_cox()`.
#' It uses a kernel-smoothed martingale estimating equation and typically runs faster
#' than the general `skmle(s = 1)` fit because it solves a more specialized problem.
#'
#' @references
#' Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao.
#' "Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
#' *Journal of the American Statistical Association* (2025): 1-12.
#'
#' Sun, Dayu, Hongyuan Cao, and Yining Chen. "Additive hazards models with sparse
#' longitudinal covariates." *Lifetime Data Analysis* (2022).
#'
#' @return
#' An object of class `kee` containing coefficient estimates, an estimated
#' variance-covariance matrix, intermediate matrices used for sandwich variance
#' estimation, and model metadata.
#'
#' @examples
#' \donttest{
#' library(survival)
#'
#' set.seed(123)
#' dat <- sim_skmle_data(
#'   n = 80,
#'   mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'   mu_bar = 8,
#'   alpha = function(tt) 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'   beta = c(1, -0.5),
#'   s = 1,
#'   cen = 0.7
#' )
#'
#' fit_add <- kee_additive(
#'   Surv(X, delta) ~ covariates,
#'   data = dat,
#'   id = id,
#'   obs_times = obs_times,
#'   h = 0.5
#' )
#'
#' summary(fit_add)
#' }
#'
#' @importFrom survival Surv
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom gaussquad legendre.quadrature.rules
#' @export
kee_additive <- function(formula, data, id, obs_times, h, lq_nodes = 64) {
    # basic validation
    if (missing(formula) || missing(data) || missing(id) || missing(obs_times) || missing(h)) {
        stop("formula, data, id, obs_times and h must all be supplied")
    }
    if (!is.numeric(h) || length(h) != 1 || h <= 0) {
        stop("bandwidth 'h' must be a positive number")
    }
    if (!is.numeric(lq_nodes) || length(lq_nodes) != 1 || lq_nodes <= 0) {
        stop("'lq_nodes' must be a positive integer")
    }
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame")
    }

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
    if (ncol(Z) > 0 && colnames(Z)[1] == "(Intercept)") Z <- Z[, -1, drop = FALSE]
    if (ncol(Z) == 0) stop("model must contain at least one covariate")

    id_raw <- model.extract(mf, "id")
    obs_times_vec <- as.numeric(model.extract(mf, "obs_times"))

    if (length(id_raw) != length(X_time) || length(obs_times_vec) != length(X_time)) {
        stop("Length of 'id' and 'obs_times' must match number of rows in data/formula")
    }
    id_vec <- as.integer(factor(id_raw))

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
    var_est <- safe_sandwich(A, Sigma, scale = 1, what = "kee_additive variance")

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
#' @return `x`, invisibly. Called for its side effect of printing a brief
#'   summary of the call and estimated coefficients.
#' @examples
#' \dontrun{
#' fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
#'                obs_times = obs_times, h = 0.5)
#' print(fit)
#' }
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
#' @return An object of class `summary.kee` containing the call,
#'   a coefficient table with estimates, standard errors, z-statistics
#'   and p-values, the convergence code, and the sample size `n`.
#' @examples
#' \dontrun{
#' fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
#'                obs_times = obs_times, h = 0.5)
#' summary(fit)
#' }
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
#' @return `x`, invisibly. Called for its side effect of printing the
#'   summary table.
#' @examples
#' \dontrun{
#' fit <- kee_cox(Surv(X, delta) ~ covariates, data = dat, id = id,
#'                obs_times = obs_times, h = 0.5)
#' print(summary(fit))
#' }
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

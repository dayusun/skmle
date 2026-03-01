#' Generate Non-Homogeneous Poisson Process
#'
#' @description Generates event times from a Non-Homogeneous Poisson Process (NHPP)
#' with a specified intensity function. Used internally for simulating
#' observation times.
#'
#' @param lambdat A function specifying the intensity rate at time `t`.
#' @param censor The maximum follow-up or censoring time. Event times are generated up to this point.
#' @param lambda.bar An upper bound for the intensity function `lambdat` over the interval `[0, censor]`.
#' Needs to be strictly greater than or equal to `max(lambdat(t))` for valid generation.
#'
#' @return A numeric vector of sorted generated event times, or an empty vector if no events occurred.
#'
#' @importFrom stats rpois runif
#' @noRd
NonHPP_gen <- function(lambdat, censor, lambda.bar) {
    nn <- 0
    while (nn <= 0) {
        nn <- rpois(1, censor * lambda.bar)
        if (nn > 0) {
            tt <- sort(runif(nn, min = 0, max = censor))

            ind <- (sapply(tt, lambdat) / lambda.bar) > runif(nn)
            tt <- tt[ind]
        } else {
            tt <- numeric(0)
        }
        nn <- length(tt)
    }

    tt
}

#' Simulate data for Transformed Hazards Models
#'
#' @description Simulates right-censored survival data with sparsely and
#' intermittently observed longitudinal covariates under a transformed hazards model.
#'
#' @param n The number of subjects to simulate data for.
#' @param mu A function specifying the intensity rate for the Non-Homogeneous Poisson Process generating the observation times.
#' @param mu_bar An upper bound for the intensity function `mu`.
#' @param alpha A function specifying the baseline hazard component.
#' @param beta A numeric vector of true regression coefficients for the generated covariates (length 2 in this implementation).
#'   The simulation currently constructs precisely two covariates;
#'   supplying a vector of other length will trigger an error.
#' @param s The Box-Cox transformation parameter (0 for proportional hazards, 1 for additive hazards).
#' @param cen The maximum follow-up or censoring time parameter.
#' @param nstep Number of steps for generating the piecewise constant longitudinal covariate process. Default is 20.
#'
#' @details
#' `sim_skmle_data` generates survival data and sparsely observed longitudinal covariates.
#' The covariate process is generated as a piece-wise constant step function. Survival times
#' are simulated utilizing the specified transformed hazards model and the `gaussquad` numerical
#' integration technique to invert the cumulative hazard. Observation times for the longitudinal
#' covariates are generated using a Non-Homogeneous Poisson Process.
#'
#' @return A tibble (data frame) in long format containing the following columns:
#'   \item{id}{Subject identifier.}
#'   \item{X}{The observed survival or censoring time.}
#'   \item{delta}{Event indicator (TRUE if event occurred, FALSE if censored).}
#'   \item{covariates}{A two‑column matrix (time‑varying covariate in column 1,
#'     constant indicator in column 2) recording the observed covariates.}
#'   \item{obs_times}{The exact times at which the covariates were observed.}
#'   \item{censoring}{The randomly generated right-censoring time for the subject.}
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' # Simulate a dataset of 200 subjects under the proportional hazards model (s = 0)
#' set.seed(123)
#' sim_data <- sim_skmle_data(
#'     n = 200,
#'     mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'     mu_bar = 8,
#'     alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'     beta = c(1, -0.5),
#'     s = 0,
#'     cen = 0.7
#' )
#' head(sim_data)
#' }
#'
#' @importFrom stats pnorm runif stepfun
#' @importFrom MASS mvrnorm
#' @importFrom nleqslv nleqslv
#' @importFrom gaussquad legendre.quadrature legendre.quadrature.rules
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @export
sim_skmle_data <- function(n,
                           mu,
                           mu_bar,
                           alpha,
                           beta,
                           s,
                           cen,
                           nstep = 20) {
    ## input validation -------------------------------------------------------
    if (!is.numeric(n) || length(n) != 1 || n <= 0) {
        stop("'n' must be a positive integer")
    }
    if (!is.function(mu) || !is.function(alpha)) {
        stop("'mu' and 'alpha' must be functions of time")
    }
    if (!is.numeric(mu_bar) || length(mu_bar) != 1 || mu_bar <= 0) {
        stop("'mu_bar' must be a positive numeric value")
    }
    if (!is.numeric(beta) || length(beta) < 1) {
        stop("'beta' must be a numeric vector of length at least 1")
    }
    if (length(beta) != 2L) {
        stop("This simulation currently only supports two covariates (length(beta) == 2)")
    }
    if (!is.numeric(s) || length(s) != 1) {
        stop("'s' must be a numeric scalar")
    }
    if (!is.numeric(cen) || length(cen) != 1 || cen <= 0) {
        stop("'cen' must be a positive numeric value")
    }
    if (!is.numeric(nstep) || length(nstep) != 1 || nstep <= 0) {
        stop("'nstep' must be a positive integer")
    }

    # internal helper --------------------------------------------------------
    sim_single_subject <- function() {
        # 1. Generate censoring times (uniform over [cen,1] but not exceeding 1)
        cen_i <- runif(1, cen, 1.5)
        cen_i <- min(cen_i, 1)

        # 2. Generate Z(t) and h(t) as a step function
        Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)

        # Generate the first covariate (time-varying step function)
        z <- 2 * (pnorm(c(MASS::mvrnorm(1, rep(0, nstep), Sigmat_z))) - 0.5)
        left_time_points <- (0:(nstep - 1)) / nstep
        z_fun <- stepfun(left_time_points, c(0, z))

        # Generate the second covariate as a binary indicator
        z_cons_raw <- as.numeric((mean(z) + rnorm(1, 0, 1)) > 0)

        # The aggregated linear predictor part (used for survival time)
        h_fun <- function(x) {
            beta[1] * z_fun(x) + beta[2] * z_cons_raw
        }

        # 3. Generate failure time based on transformation parameter 's'
        if (s == 0) {
            lam_fun <- function(tt) exp(alpha(tt) + h_fun(tt))
        } else {
            lam_fun <- function(tt) (s * (alpha(tt) + h_fun(tt)) + 1)^(1 / s)
        }

        u <- runif(1)
        lqrule64 <- gaussquad::legendre.quadrature.rules(64)[[64]]

        # Solve for failure time using nleqslv and Legendre-Gauss Quadrature
        fail_time <- nleqslv::nleqslv(cen_i / 2, function(ttt) {
            gaussquad::legendre.quadrature(lam_fun, lower = 0, upper = ttt, lqrule64) + log(u)
        })$x

        X <- min(fail_time, cen_i)

        # 4. Generate R_ik (Observation times) via NHPP
        obs_times <- NonHPP_gen(mu, cen_i, mu_bar)
        if (length(obs_times) == 0) obs_times <- cen_i

        # Combine covariates evaluated at observation times
        covariates_obscov <- cbind(
            z_fun(obs_times),
            rep(z_cons_raw, length(obs_times))
        )
        colnames(covariates_obscov) <- paste0("cov", seq_len(ncol(covariates_obscov)))

        # Return results for the subject
        tibble::tibble(
            X = X,
            delta = fail_time < cen_i,
            covariates = covariates_obscov,
            obs_times = obs_times,
            censoring = cen_i
        )
    }

    # Simulate for 'n' subjects and bind rows
    res_list <- replicate(n, sim_single_subject(), simplify = FALSE)
    sim_data <- dplyr::bind_rows(res_list, .id = "id")

    return(sim_data)
}

# broom generics.
#
# Re-exported from 'generics' rather than taken as a dependency on 'broom'
# itself, which is the usual arrangement: a modelling package supplies the
# methods, broom supplies nothing but the names.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment


#' Summarise a fit as a tibble
#'
#' One row per coefficient for `skmle` and `kee` fits, and one row per
#' (target time, coefficient) for `kee_td`, so a coefficient curve goes
#' straight into `ggplot2` without reshaping.
#'
#' @param x A fitted `skmle`, `kee` or `kee_td` object.
#' @param conf.int Logical, add `conf.low` and `conf.high` columns.
#' @param conf.level Confidence level for those columns.
#' @param ... Unused.
#'
#' @return A tibble with columns `term`, `estimate`, `std.error`, `statistic`
#'   and `p.value`, preceded by `time` for a `kee_td` fit, and followed by
#'   `conf.low` and `conf.high` when `conf.int = TRUE`.
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 150)
#' fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
#' tidy(fit, conf.int = TRUE)
#' @importFrom stats qnorm pnorm
#' @name tidy.skmle
#' @export
tidy.skmle <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
    .tidy_flat(stats::coef(x), sqrt(diag(x$var)), conf.int, conf.level)
}

#' @rdname tidy.skmle
#' @export
tidy.kee <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
    .tidy_flat(stats::coef(x), sqrt(diag(x$var)), conf.int, conf.level)
}

#' @rdname tidy.skmle
#' @export
tidy.kee_td <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
    nms <- colnames(x$coefficients)
    out <- tibble::tibble(
        time = rep(x$times, times = length(nms)),
        term = rep(nms, each = length(x$times)),
        estimate = as.numeric(x$coefficients),
        std.error = as.numeric(x$se),
        statistic = as.numeric(x$z),
        p.value = as.numeric(x$p.value)
    )
    if (conf.int) {
        zq <- qnorm(1 - (1 - conf.level) / 2)
        out$conf.low <- out$estimate - zq * out$std.error
        out$conf.high <- out$estimate + zq * out$std.error
    }
    out
}

.tidy_flat <- function(est, se, conf.int, conf.level) {
    z <- est / se
    out <- tibble::tibble(
        term = names(est),
        estimate = as.numeric(est),
        std.error = as.numeric(se),
        statistic = as.numeric(z),
        p.value = 2 * pnorm(-abs(as.numeric(z)))
    )
    if (conf.int) {
        zq <- qnorm(1 - (1 - conf.level) / 2)
        out$conf.low <- out$estimate - zq * out$std.error
        out$conf.high <- out$estimate + zq * out$std.error
    }
    out
}


#' One-row summary of a fit
#'
#' @param x A fitted `skmle`, `kee` or `kee_td` object.
#' @param ... Unused.
#'
#' @return A one-row tibble. `nobs` is the number of **subjects**, which is the
#'   unit the asymptotics are in. No `AIC` or `logLik` column is reported: the
#'   criterion these estimators optimise is a kernel-weighted pseudo-likelihood,
#'   so an information criterion built from it would not mean what it usually
#'   means.
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 150)
#' fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
#' glance(fit)
#' @name glance.skmle
#' @export
glance.skmle <- function(x, ...) {
    tibble::tibble(
        nobs = as.integer(x$n),
        nterms = length(stats::coef(x)),
        s = x$s %||% NA_real_,
        h = x$h,
        convergence = as.integer(x$convergence %||% NA_integer_)
    )
}

#' @rdname glance.skmle
#' @export
glance.kee <- function(x, ...) {
    tibble::tibble(
        nobs = as.integer(x$n),
        nterms = length(stats::coef(x)),
        h = x$h %||% NA_real_,
        npair = as.integer(x$npair %||% NA_integer_),
        link = x$link %||% NA_character_,
        one_sided = x$one_sided %||% NA,
        convergence = as.integer(x$convergence %||% NA_integer_)
    )
}

#' @rdname glance.skmle
#' @export
glance.kee_td <- function(x, ...) {
    tibble::tibble(
        nobs = as.integer(x$n),
        nterms = ncol(x$coefficients),
        ntimes = length(x$times),
        h = x$h,
        h2 = x$h2,
        link = x$link,
        one_sided = x$one_sided,
        nconverged = sum(x$convergence == 0L)
    )
}


#' Add fitted means to the covariate table
#'
#' Fitted means belong to covariate occasions. The estimating equation
#' evaluates the link at each observed covariate vector, and the kernel weight
#' is what ties that vector to a response occasion. So `augment()` returns
#' `data_x` with a `.fitted` column, \eqn{g(X_i(S_{ik})^\top \hat\beta)}.
#'
#' There is deliberately no `.resid`. A residual needs a response value at
#' \eqn{S_{ik}}, and asynchronous data has none; any residual reported here
#' would have to be invented.
#'
#' The fit does not retain its data, so `data_x` must be supplied.
#'
#' @param x A `kee_async` fit.
#' @param data_x The covariate table the fit was made from.
#' @param ... Unused.
#'
#' @return `data_x` as a tibble with a `.fitted` column added.
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 150)
#' fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
#' augment(fit, d$x)
#' @export
augment.kee_async <- function(x, data_x, ...) {
    if (missing(data_x)) {
        stop(
            "'data_x' must be supplied: the fit does not retain the data it ",
            "was built from",
            call. = FALSE
        )
    }
    if (!is.data.frame(data_x)) stop("'data_x' must be a data.frame", call. = FALSE)

    b <- stats::coef(x)
    nms <- setdiff(names(b), "(Intercept)")
    missing_cols <- setdiff(nms, names(data_x))
    if (length(missing_cols)) {
        stop("'data_x' is missing the fitted covariates: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    eta <- if ("(Intercept)" %in% names(b)) b[["(Intercept)"]] else 0
    for (nm in nms) eta <- eta + b[[nm]] * data_x[[nm]]
    out <- tibble::as_tibble(data_x)
    out$.fitted <- .async_mu(eta, .async_link_code(x$link))
    out
}

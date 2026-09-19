#' Kernel weight functions
#'
#' @description
#' A weight function says how much a covariate observation at lag `u` from a
#' time contributes to the estimating equation at that time, where `u` is the
#' lag divided by the bandwidth. Every estimator in this package takes one
#' through its `weight` argument, and any R function of `u` will do: the
#' Epanechnikov kernels here are a starting point, not the only option.
#'
#' @details
#' # Writing your own
#'
#' A weight is a function of one argument, the standardised lag
#' \eqn{u = (t - r)/h}, returning a value of the same shape. It is called on
#' vectors and on matrices, and it must preserve the shape it is given: the
#' sieve quadrature evaluates it on a matrix of node-by-observation lags. Zero
#' outside its own support is the function's own responsibility.
#'
#' ```r
#' gauss_half <- function(u) exp(-u^2 / 2) * (u >= 0 & u <= 3)
#' kee_cox(Surv(X, delta) ~ z, data, id = id, obs_times = r,
#'         h = 0.3, weight = gauss_half)
#' ```
#'
#' Nothing about the function has to be polynomial, nonnegative or symmetric.
#' A weight that takes negative values is supported throughout: the estimating
#' equations admit a row on whether its weight is nonzero, not on whether it is
#' positive.
#'
#' `support` is the one attribute worth setting, as a length-2 numeric giving
#' the interval outside which the weight is zero. It is used to enumerate pairs
#' in the asynchronous estimators, where scanning every pair would otherwise be
#' quadratic. Without it a conservative default is assumed and the fit is
#' correct but slower.
#'
#' # Scaling
#'
#' The weight is used as \eqn{W(u)/h}. A weight that integrates to something
#' other than 1 is not wrong: the normalisation cancels wherever the weight
#' appears in a ratio, which is everywhere in the risk-set averages, and is
#' absorbed into `h` elsewhere. `w_epan_half()` integrates to `1/2` for that
#' reason and is the weight every published result in this package used.
#'
#' # The one-sided restriction is separate
#'
#' `one_sided` is a modelling choice, not a property of the weight. It is the
#' risk-set restriction of a hazard model: only covariate observations strictly
#' before a time may inform it. It is applied on the sign of the lag after the
#' weight has been evaluated, so a two-sided weight under `one_sided = TRUE` is
#' truncated at zero rather than rejected.
#'
#' @return A function of one argument, the standardised lag
#'   \eqn{u = (t - r)/h} as a vector or a matrix, returning the weight with the
#'   shape of its argument preserved and carrying a `support` attribute.
#'
#' @examples
#' W <- w_epan_half()
#' W(c(-0.5, 0, 0.5, 1, 1.5))
#' attr(W, "support")
#'
#' # Any function of the standardised lag is a weight.
#' tri <- function(u) pmax(1 - abs(u), 0)
#' tri(c(-0.5, 0, 0.5))
#' @name kernel-weights
NULL

#' @describeIn kernel-weights The Epanechnikov half kernel
#'   \eqn{0.75(1 - u^2)_+} on \eqn{[0, 1]}, the default for the survival
#'   estimators and the weight every published result in this package was
#'   computed with.
#' @export
w_epan_half <- function() {
    # Zeroed below 0 as well as above 1.  The survival code always multiplied
    # by (lag > 0) and so never saw the difference, but a weight that claims
    # support [0, 1] and then returns 0.53 at u = -0.5 would quietly become a
    # full kernel anywhere the one-sided restriction is off.
    W <- function(u) pmax((1 - u^2) * 0.75, 0) * (u >= 0 & u <= 1)
    attr(W, "support") <- c(0, 1)
    attr(W, "order") <- 1L
    attr(W, "coef") <- c(0.75, 0, -0.75)
    attr(W, "mirror") <- FALSE
    W
}

#' @describeIn kernel-weights The Epanechnikov full kernel
#'   \eqn{0.75(1 - u^2)_+} on \eqn{[-1, 1]}, the default for the asynchronous
#'   longitudinal estimators and the kernel of Cao, Zeng and Fine (2015).
#'   Symmetry kills the first moment, so it is order 2.
#' @export
w_epan_full <- function() {
    W <- function(u) pmax((1 - u^2) * 0.75, 0)
    attr(W, "support") <- c(-1, 1)
    attr(W, "order") <- 2L
    attr(W, "coef") <- c(0.75, 0, -0.75)
    attr(W, "mirror") <- FALSE
    W
}

#' Resolve the weight argument
#'
#' `NULL` means "the Epanechnikov implied by `one_sided`", which is what every
#' call made before `weight` existed was asking for. Resolving it here rather
#' than as a default argument keeps that backward compatibility exact: a
#' literal default of `w_epan_half()` would silently turn
#' `one_sided = FALSE` into a one-sided fit, because the support would still
#' be `[0, 1]`.
#'
#' @param weight A weight function, or `NULL`.
#' @param one_sided Logical, used only when `weight` is `NULL`.
#' @return A weight function.
#' @keywords internal
resolve_weight <- function(weight, one_sided) {
    if (is.null(weight)) {
        return(if (isTRUE(one_sided)) w_epan_half() else w_epan_full())
    }
    if (!is.function(weight)) {
        stop("'weight' must be a function of the standardised lag, or NULL",
            call. = FALSE
        )
    }
    probe <- tryCatch(weight(c(-0.5, 0, 0.5)), error = function(e) e)
    if (inherits(probe, "error")) {
        stop("'weight' could not be evaluated on a numeric vector: ",
            conditionMessage(probe),
            call. = FALSE
        )
    }
    if (length(probe) != 3L || !is.numeric(probe) || anyNA(probe)) {
        stop(
            "'weight' must return one finite numeric value per element of its ",
            "argument; got ", length(probe), " value(s) of type ",
            typeof(probe),
            call. = FALSE
        )
    }
    # The sieve quadrature calls the weight on a matrix of node-by-observation
    # lags and indexes the result as a matrix.  A weight built with an
    # apply() or a scalar if() returns a bare vector and would be silently
    # recycled into the wrong cells, so the shape is checked once, here.
    m <- matrix(c(-0.5, 0, 0.5, 1), nrow = 2L)
    shape <- tryCatch(weight(m), error = function(e) e)
    if (inherits(shape, "error") || !identical(dim(shape), dim(m))) {
        stop(
            "'weight' must preserve the shape of its argument: called on a ",
            "2x2 matrix it must return a 2x2 matrix. Build it from ",
            "vectorised arithmetic (pmax, ifelse, arithmetic operators) ",
            "rather than from apply() or a scalar if().",
            call. = FALSE
        )
    }
    weight
}

#' The support of a weight, with a conservative fallback
#'
#' @param weight A weight function.
#' @return Length-2 numeric.
#' @keywords internal
weight_support <- function(weight) {
    supp <- attr(weight, "support")
    if (is.null(supp) || length(supp) != 2L || anyNA(supp) || supp[1] >= supp[2]) {
        return(c(-1, 1))
    }
    as.numeric(supp)
}

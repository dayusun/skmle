# Shared front end for the two asynchronous longitudinal estimators.
#
# Both take the response and the covariate as separate tables, because the whole
# point is that they are not on a common grid and so cannot be a single data
# frame without inventing rows.
.async_prep <- function(y_id, y_time, y, x_id, x_time, X, intercept) {
    X <- as.matrix(X)
    if (!is.numeric(X)) stop("'X' must be numeric")
    if (length(x_id) != nrow(X) || length(x_time) != nrow(X)) {
        stop("'x_id' and 'x_time' must have one entry per row of 'X'")
    }
    if (length(y_id) != length(y_time) || length(y_id) != length(y)) {
        stop("'y_id', 'y_time' and 'y' must have equal length")
    }
    if (!nrow(X)) stop("'X' has no rows")
    if (!length(y)) stop("'y' is empty")

    nm <- colnames(X)
    if (is.null(nm)) nm <- rep("", ncol(X))
    nm[!nzchar(nm)] <- paste0("X", which(!nzchar(nm)))
    colnames(X) <- nm
    if (intercept) X <- cbind(`(Intercept)` = 1, X)
    p <- ncol(X)

    # Common subject coding across the two tables, so a subject present in only
    # one of them still occupies a row of the score matrix and contributes zero.
    lev <- sort(unique(c(as.character(y_id), as.character(x_id))))
    yi <- match(as.character(y_id), lev)
    xi <- match(as.character(x_id), lev)

    # Sort covariate rows by subject, then by time: the C++ pair enumeration
    # walks a contiguous block per subject.
    o <- order(xi, x_time)
    list(
        yi = yi, yt = as.numeric(y_time), y = as.numeric(y),
        xi = xi[o], xt = as.numeric(x_time)[o], X = X[o, , drop = FALSE],
        n = length(lev), p = p, lev = lev
    )
}

.async_link_code <- function(link) {
    switch(link,
        identity = 0L,
        log = 1L,
        logistic = 2L,
        stop("unsupported link")
    )
}

# A^-1 B A^-1 with the same NA-on-singular behaviour the rest of the package uses.
.async_sandwich <- function(bread, U, nms) {
    p <- ncol(bread)
    Ai <- tryCatch(solve(bread), error = function(e) matrix(NA_real_, p, p))
    v <- Ai %*% crossprod(U) %*% t(Ai)
    dimnames(v) <- list(nms, nms)
    v
}

.async_check_h <- function(h, what = "h") {
    if (!is.numeric(h) || length(h) != 1L || is.na(h) || h <= 0) {
        stop(sprintf("'%s' must be a single positive number", what))
    }
    invisible(TRUE)
}

.async_check_one_sided <- function(one_sided) {
    if (!is.logical(one_sided) || length(one_sided) != 1L || is.na(one_sided)) {
        stop("'one_sided' must be TRUE or FALSE")
    }
    invisible(TRUE)
}


#' Asynchronous longitudinal regression, time-invariant coefficients
#'
#' @description
#' Fit a regression of a sparsely observed longitudinal response on a sparsely
#' observed longitudinal covariate when the two are measured on **different**
#' time grids. Responses at \eqn{T_{ij}} are paired with the same subject's
#' covariate observations at \eqn{S_{ik}}, and each pair is weighted by its time
#' separation, so no pair has to be discarded and no value has to be carried
#' forward.
#'
#' This is the time-invariant coefficient estimator of Cao, Zeng and Fine
#' (2015), which solves
#' \deqn{U_n(\beta) = n^{-1} \sum_i \sum_j \sum_k W_h(T_{ij} - S_{ik})
#'       X_i(S_{ik}) [ Y_i(T_{ij}) - g\{X_i(S_{ik})^\top \beta\} ] = 0,}
#' with \eqn{W_h(t) = W(t/h)/h}. See [kee_async_td()] for the companion
#' estimator with time-dependent coefficients \eqn{\beta(t)}.
#'
#' @details
#' The link is evaluated at each **observed** covariate vector and the weight
#' multiplies the resulting contribution. The covariate is never smoothed and
#' then substituted, which would be regression calibration and is not consistent
#' under sparse sampling. The estimator converges at the univariate smoothing
#' rate \eqn{(nh)^{1/2}}, not at \eqn{\sqrt n}.
#'
#' **Half kernel or full kernel.** `one_sided = FALSE`, the default, is the full
#' two-sided Epanechnikov kernel of the paper: a response is paired with
#' covariate observations on either side of it. `one_sided = TRUE` admits only
#' covariate observations that strictly precede the response, which is what a
#' causal reading of the covariate path requires and what the survival
#' estimators in this package do by default.
#'
#' Cao, Zeng and Fine assume the covariance of the covariate process is twice
#' differentiable and obtain a bias of order \eqn{h^2}. Their Section 6 notes
#' that relaxing this to one-sided differentiability, which admits processes
#' with independent increments such as the Ornstein-Uhlenbeck process, leaves a
#' bias of order \eqn{h} instead. Watch for an estimate that drifts as `h`
#' grows.
#'
#' Pairs are enumerated in C++ and then collapsed onto covariate rows before the
#' coefficients are touched, so the full pair design is never materialised and
#' the Newton iterations cost \eqn{O(n_x p^2)} rather than \eqn{O(n_{pair} p^2)}.
#'
#' @param y_id,y_time,y Subject identifier, observation time and value for each
#'   response occasion. Vectors of equal length.
#' @param x_id,x_time Subject identifier and observation time for each covariate
#'   occasion. Vectors of equal length, matching the rows of `X`.
#' @param X Covariate matrix, one row per covariate occasion.
#' @param h Positive bandwidth.
#' @param one_sided Logical. `FALSE` (the default) is the full two-sided
#'   Epanechnikov kernel of the paper. `TRUE` is the half kernel: only covariate
#'   observations strictly before the response contribute.
#' @param link One of `"identity"` (closed form), `"log"` or `"logistic"`
#'   (Newton-Raphson).
#' @param intercept Logical, prepend a column of ones to `X`.
#' @param maxit,tol Newton-Raphson controls, ignored when `link = "identity"`.
#'
#' @return An object of class `kee` with `coefficients`, the sandwich `var`,
#'   the bread `A`, the meat `Sigma`, `n` (subjects), `npair` (weighted pairs
#'   contributing) and the call.
#'
#' @references
#' Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
#' asynchronous longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B} 77, 755-776.
#'
#' @seealso [kee_async_td()] for time-dependent coefficients.
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 100, beta = c(0.5, 1.5))
#' fit <- kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x, h = 0.25)
#' summary(fit)
#'
#' # Half kernel: only covariate observations that precede the response.
#' kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
#'   h = 0.25, one_sided = TRUE
#' )
#' @export
kee_async <- function(y_id, y_time, y, x_id, x_time, X, h, one_sided = FALSE,
                      link = c("identity", "log", "logistic"),
                      intercept = TRUE, maxit = 50L, tol = 1e-8) {
    link <- match.arg(link)
    .async_check_h(h)
    .async_check_one_sided(one_sided)
    supp <- if (one_sided) c(0, 1) else c(-1, 1)
    d <- .async_prep(y_id, y_time, y, x_id, x_time, X, intercept)

    # Half-open [lo, hi) block of sorted covariate rows owned by each subject,
    # 0-based for C++.
    lo <- findInterval(seq_len(d$n) - 1L, d$xi)
    hi <- findInterval(seq_len(d$n), d$xi)

    pr <- async_pairs_cpp(
        as.integer(d$yi - 1L), d$yt, d$xt,
        as.integer(lo), as.integer(hi), h, supp[1], supp[2]
    )
    if (!length(pr$j)) {
        stop("no response occasion has a covariate observation inside the window")
    }

    w <- epan_kernel(pr$z) / h
    # Strictly before, matching the survival side: a covariate observed at
    # exactly the response time is not "in the past".
    if (one_sided) w[pr$z <= 0] <- 0
    ac <- async_pair_accum_cpp(pr$j, pr$k, w, d$y, nrow(d$X))
    fit <- async_solve_cpp(
        d$X, as.integer(d$xi - 1L), ac$omega, ac$c, d$n,
        .async_link_code(link), as.integer(maxit), tol, numeric(d$p)
    )

    if (fit$convergence == 1L) {
        stop("the weighted bread matrix is singular; widen 'h' or drop collinear covariates")
    }
    if (fit$convergence == 2L) {
        warning("Newton-Raphson reached 'maxit' without meeting 'tol'")
    }

    beta <- as.numeric(fit$beta)
    names(beta) <- colnames(d$X)

    out <- list(
        coefficients = beta,
        var = .async_sandwich(fit$bread, fit$U, colnames(d$X)),
        W = fit$bread,
        Sigma = crossprod(fit$U),
        A = fit$bread,
        convergence = fit$convergence,
        n = d$n,
        npair = if (one_sided) sum(pr$z > 0) else length(pr$j),
        h = h,
        one_sided = one_sided,
        link = link,
        call = match.call()
    )
    class(out) <- "kee"
    out
}


#' Asynchronous longitudinal regression, time-dependent coefficients
#'
#' @description
#' Fit \eqn{E\{Y(t) \mid X(t)\} = g\{X(t)^\top \beta(t)\}} from asynchronous
#' sparse longitudinal data, estimating the coefficient curve \eqn{\beta(t)}
#' pointwise. This is the time-dependent coefficient estimator of Cao, Zeng and
#' Fine (2015), the companion to [kee_async()].
#'
#' At a target time \eqn{t} the estimating equation weights a response occasion
#' and a covariate occasion by their separate distances to \eqn{t},
#' \deqn{U_n\{\beta(t)\} = n^{-1} \sum_i \sum_j \sum_k W_{h_1}(T_{ij} - t)
#'       W_{h_2}(S_{ik} - t) X_i(S_{ik})
#'       [ Y_i(T_{ij}) - g\{X_i(S_{ik})^\top \beta(t)\} ] = 0.}
#'
#' @details
#' Because both time arguments are smoothed, the estimator converges at the
#' **bivariate** smoothing rate \eqn{(n h_1 h_2)^{1/2}}, slower than the
#' \eqn{(nh)^{1/2}} of the time-invariant fit and much slower than the usual
#' varying-coefficient rate available under synchronous sampling. Expect wide
#' bands and choose bandwidths accordingly.
#'
#' The product weight factorises over the two occasion indices, so no pair is
#' ever enumerated: the response side collapses to two per-subject scalars and
#' the covariate side to one weight per row. Each target time then costs
#' \eqn{O((n_y + n_x p^2))}. Fits at successive target times are warm-started
#' from the previous solution when the link is nonlinear.
#'
#' Standard errors are the pointwise sandwich at each target time; they are not
#' simultaneous bands, and no undersmoothing correction is applied, so the
#' intervals are centred on the smoothed curve rather than on \eqn{\beta(t)}.
#'
#' @inheritParams kee_async
#' @param times Numeric vector of target times at which to estimate
#'   \eqn{\beta(t)}.
#' @param h Bandwidth for the response side, \eqn{h_1}.
#' @param h2 Bandwidth for the covariate side, \eqn{h_2}. Defaults to `h`.
#' @param one_sided Logical. `FALSE` (the default) is the full two-sided kernel.
#'   `TRUE` restricts the **covariate** side to observations strictly before the
#'   target time, so \eqn{\beta(t)} is estimated from covariate values that were
#'   already available at `t`. The response side is always two-sided: it is a
#'   local average around `t`, not a filtering step. Both lags are measured as
#'   "how far in the past", matching [kee_async()] and the survival
#'   estimators.
#'
#' @return An object of class `kee_td`:
#' \describe{
#'   \item{coefficients}{`length(times)` by `p` matrix of \eqn{\hat\beta(t)}.}
#'   \item{se, z, p.value}{Pointwise standard errors and Wald statistics, same
#'     shape as `coefficients`.}
#'   \item{var}{List of \eqn{p \times p} sandwich matrices, one per target time.}
#'   \item{times, h, h2, link, n, call}{Fit metadata.}
#'   \item{convergence}{Integer per target time: 0 converged, 1 singular,
#'     2 hit `maxit`, 3 no data in the window.}
#' }
#'
#' @references
#' Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
#' asynchronous longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B} 77, 755-776.
#'
#' @seealso [kee_async()] for time-invariant coefficients, [plot.kee_td()].
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
#' fit <- kee_async_td(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
#'   times = seq(0.2, 0.8, by = 0.1), h = 0.3
#' )
#' coef(fit)
#' @export
kee_async_td <- function(y_id, y_time, y, x_id, x_time, X, times, h, h2 = h,
                         one_sided = FALSE,
                         link = c("identity", "log", "logistic"),
                         intercept = TRUE, maxit = 50L, tol = 1e-8) {
    link <- match.arg(link)
    .async_check_one_sided(one_sided)
    .async_check_h(h)
    .async_check_h(h2, "h2")
    if (!is.numeric(times) || !length(times) || anyNA(times)) {
        stop("'times' must be a non-empty numeric vector without NA")
    }
    d <- .async_prep(y_id, y_time, y, x_id, x_time, X, intercept)

    lk <- .async_link_code(link)
    nt <- length(times)
    nms <- colnames(d$X)

    beta <- se <- matrix(NA_real_, nt, d$p, dimnames = list(NULL, nms))
    vars <- vector("list", nt)
    conv <- integer(nt)
    npair <- integer(nt)
    start <- numeric(d$p)

    for (m in seq_len(nt)) {
        # Lags are measured as "how far in the past", the same convention the
        # survival code and kee_async() use, so a half-support weight means the
        # same thing in all three places.
        wy <- epan_kernel((times[m] - d$yt) / h) / h
        lag_x <- times[m] - d$xt
        vx <- epan_kernel(lag_x / h2) / h2
        if (one_sided) vx[lag_x <= 0] <- 0

        # A subject contributes only if it has weight on both sides.
        ac <- async_td_accum_cpp(
            as.integer(d$yi - 1L), wy, d$y, as.integer(d$xi - 1L), vx, d$n
        )
        npair[m] <- sum(ac$omega != 0 | ac$c != 0)
        if (!npair[m]) {
            conv[m] <- 3L
            next
        }

        fit <- async_solve_cpp(
            d$X, as.integer(d$xi - 1L), ac$omega, ac$c, d$n,
            lk, as.integer(maxit), tol, start
        )
        conv[m] <- fit$convergence
        if (fit$convergence == 1L) next

        b <- as.numeric(fit$beta)
        beta[m, ] <- b
        if (lk != 0L && all(is.finite(b))) start <- b
        v <- .async_sandwich(fit$bread, fit$U, nms)
        vars[[m]] <- v
        dv <- diag(v)
        se[m, ] <- ifelse(is.na(dv) | dv < 0, NA_real_, sqrt(pmax(dv, 0)))
    }

    if (all(conv == 3L)) {
        stop("no target time has data in its window; widen 'h'/'h2' or move 'times'")
    }
    if (any(conv == 2L)) {
        warning("Newton-Raphson reached 'maxit' at ", sum(conv == 2L), " target time(s)")
    }

    z <- beta / se
    out <- list(
        coefficients = beta,
        se = se,
        z = z,
        p.value = 2 * pnorm(-abs(z)),
        var = vars,
        times = times,
        h = h,
        h2 = h2,
        one_sided = one_sided,
        link = link,
        n = d$n,
        nobs = c(response = length(d$y), covariate = nrow(d$X)),
        nactive = npair,
        convergence = conv,
        call = match.call()
    )
    class(out) <- "kee_td"
    out
}


#' @param x A `kee_td` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @rdname kee_async_td
#' @export
print.kee_td <- function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nTime-dependent coefficients, ", x$link, " link\n", sep = "")
    cat("subjects: ", x$n, "   bandwidths: h = ", format(x$h),
        ", h2 = ", format(x$h2), "\n\n",
        sep = ""
    )
    tab <- cbind(t = x$times, x$coefficients)
    print(tab, digits = 4)
    bad <- x$convergence != 0L
    if (any(bad)) {
        cat("\n", sum(bad), " of ", length(bad),
            " target times did not produce an estimate\n",
            sep = ""
        )
    }
    invisible(x)
}


#' Plot estimated coefficient curves
#'
#' Draws \eqn{\hat\beta_j(t)} against `t` with pointwise 95% Wald bands, one
#' panel per coefficient. The bands are pointwise and are not corrected for
#' smoothing bias, so they cover \eqn{E\hat\beta(t)}, not \eqn{\beta(t)}.
#'
#' @param x A `kee_td` object from [kee_async_td()].
#' @param which Integer or character vector selecting coefficients to draw.
#'   Defaults to all.
#' @param level Confidence level for the pointwise bands.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
#' fit <- kee_async_td(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
#'   times = seq(0.2, 0.8, by = 0.05), h = 0.3
#' )
#' plot(fit)
#' @importFrom graphics abline lines par polygon
#' @importFrom stats qnorm pnorm
#' @export
plot.kee_td <- function(x, which = seq_len(ncol(x$coefficients)),
                        level = 0.95, ...) {
    if (!inherits(x, "kee_td")) stop("object must be of class 'kee_td'")
    if (is.character(which)) which <- match(which, colnames(x$coefficients))
    which <- which[!is.na(which)]
    if (!length(which)) stop("'which' selects no coefficient")
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
        stop("'level' must be a single number in (0, 1)")
    }

    zq <- qnorm(1 - (1 - level) / 2)
    nc <- ceiling(sqrt(length(which)))
    op <- par(mfrow = c(ceiling(length(which) / nc), nc))
    on.exit(par(op), add = TRUE)

    for (j in which) {
        b <- x$coefficients[, j]
        s <- x$se[, j]
        ok <- is.finite(b)
        lo <- b - zq * s
        hi <- b + zq * s
        rng <- range(c(b[ok], lo[is.finite(lo)], hi[is.finite(hi)]), finite = TRUE)
        plot(x$times, b,
            type = "n", ylim = rng,
            xlab = "t", ylab = colnames(x$coefficients)[j],
            main = colnames(x$coefficients)[j], ...
        )
        band <- ok & is.finite(lo) & is.finite(hi)
        if (any(band)) {
            polygon(c(x$times[band], rev(x$times[band])),
                c(lo[band], rev(hi[band])),
                col = "grey85", border = NA
            )
        }
        abline(h = 0, col = "grey50", lty = 3)
        lines(x$times[ok], b[ok], lwd = 2, col = "#1f5aa6")
    }
    invisible(x)
}

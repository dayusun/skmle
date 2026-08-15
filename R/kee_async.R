# Front end for the asynchronous longitudinal estimators.
#
# The response and the covariate arrive as two separate tables, because the
# whole point of the setting is that they are not on a common grid and so
# cannot be one data frame without inventing rows.  The formula spans both:
# its left-hand side is resolved in `data_y`, its right-hand side in `data_x`.

# `id` and `time` may be given unquoted (id = subject), as a string
# (id = "subject"), or embraced from an enclosing function (id = {{ col }}).
# Capturing with enquo() rather than substitute() is what makes the third work:
# substitute() would see the literal `{{ col }}` and deparse it.
.async_name <- function(q) {
    if (rlang::quo_is_symbol(q)) {
        return(rlang::as_name(q))
    }
    v <- try(rlang::eval_tidy(q), silent = TRUE)
    if (is.character(v) && length(v) == 1L && !is.na(v)) {
        return(v)
    }
    stop(
        "column arguments must be a column name, given unquoted, as a string, ",
        "or embraced with {{ }}",
        call. = FALSE
    )
}

.async_col <- function(data, nm, what, side) {
    if (!nm %in% names(data)) {
        stop(sprintf(
            "'%s' column \"%s\" not found in %s (columns: %s)",
            what, nm, side, paste(names(data), collapse = ", ")
        ), call. = FALSE)
    }
    data[[nm]]
}

# Parse the two-table formula into the vectors the C++ backend consumes.
.async_parse <- function(formula, data_y, data_x, id_nm, time_nm, intercept) {
    if (!inherits(formula, "formula") || length(formula) != 3L) {
        stop("'formula' must be two-sided, e.g. y ~ x1 + x2", call. = FALSE)
    }
    if (!is.data.frame(data_y)) stop("'data_y' must be a data.frame", call. = FALSE)
    if (!is.data.frame(data_x)) stop("'data_x' must be a data.frame", call. = FALSE)

    y_id <- .async_col(data_y, id_nm, "id", "data_y")
    y_time <- .async_col(data_y, time_nm, "time", "data_y")
    x_id <- .async_col(data_x, id_nm, "id", "data_x")
    x_time <- .async_col(data_x, time_nm, "time", "data_x")

    y <- eval(formula[[2L]], data_y, environment(formula))
    if (!is.numeric(y)) stop("the response must be numeric", call. = FALSE)
    if (length(y) != nrow(data_y)) {
        stop("the response must have one value per row of 'data_y'", call. = FALSE)
    }

    # Covariates come from data_x only.  na.pass keeps the rows aligned with
    # x_id and x_time so the three can be filtered together; letting
    # model.matrix drop rows silently would misalign them.
    rhs <- stats::reformulate(attr(stats::terms(formula, data = data_x), "term.labels"))
    mfx <- stats::model.frame(rhs, data = data_x, na.action = stats::na.pass)
    X <- stats::model.matrix(rhs, data = mfx)
    if (ncol(X) && colnames(X)[1L] == "(Intercept)") X <- X[, -1L, drop = FALSE]
    if (!ncol(X)) stop("the model must contain at least one covariate", call. = FALSE)

    ok_y <- !is.na(y) & !is.na(y_id) & !is.na(y_time)
    ok_x <- stats::complete.cases(X) & !is.na(x_id) & !is.na(x_time)
    if (!any(ok_y)) stop("every row of 'data_y' has a missing value", call. = FALSE)
    if (!any(ok_x)) stop("every row of 'data_x' has a missing value", call. = FALSE)

    .async_prep(
        y_id[ok_y], y_time[ok_y], y[ok_y],
        x_id[ok_x], x_time[ok_x], X[ok_x, , drop = FALSE], intercept
    )
}

.async_prep <- function(y_id, y_time, y, x_id, x_time, X, intercept) {
    X <- as.matrix(X)
    if (!is.numeric(X)) stop("covariates must be numeric", call. = FALSE)
    if (!nrow(X)) stop("no covariate rows remain", call. = FALSE)
    if (!length(y)) stop("no response rows remain", call. = FALSE)

    nm <- colnames(X)
    if (is.null(nm)) nm <- rep("", ncol(X))
    nm[!nzchar(nm)] <- paste0("X", which(!nzchar(nm)))
    colnames(X) <- nm
    if (intercept) X <- cbind(`(Intercept)` = 1, X)

    # Common subject coding across the two tables, so a subject present in only
    # one of them still occupies a row of the score matrix and contributes zero.
    lev <- sort(unique(c(as.character(y_id), as.character(x_id))))
    yi <- match(as.character(y_id), lev)
    xi <- match(as.character(x_id), lev)

    # Covariate rows sorted by subject then time: the C++ pair enumeration
    # walks a contiguous block per subject.
    o <- order(xi, x_time)
    list(
        yi = yi, yt = as.numeric(y_time), y = as.numeric(y),
        xi = xi[o], xt = as.numeric(x_time)[o], X = X[o, , drop = FALSE],
        n = length(lev), p = ncol(X), lev = lev
    )
}

# Restrict a prepared object to a subset of subjects, renumbering the codes.
.async_subset <- function(d, keep) {
    map <- integer(d$n)
    map[keep] <- seq_along(keep)
    ky <- map[d$yi] > 0L
    kx <- map[d$xi] > 0L
    list(
        yi = map[d$yi][ky], yt = d$yt[ky], y = d$y[ky],
        xi = map[d$xi][kx], xt = d$xt[kx], X = d$X[kx, , drop = FALSE],
        n = length(keep), p = d$p, lev = d$lev[keep]
    )
}

.async_link_code <- function(link) {
    switch(link,
        identity = 0L,
        log = 1L,
        logistic = 2L,
        stop("unsupported link", call. = FALSE)
    )
}

.async_mu <- function(eta, lk) {
    if (lk == 1L) {
        exp(eta)
    } else if (lk == 2L) {
        stats::plogis(eta)
    } else {
        eta
    }
}

# A^-1 B A^-1, with the same NA-on-singular behaviour as the rest of the package.
.async_sandwich <- function(bread, U, nms) {
    p <- ncol(bread)
    Ai <- tryCatch(solve(bread), error = function(e) matrix(NA_real_, p, p))
    v <- Ai %*% crossprod(U) %*% t(Ai)
    dimnames(v) <- list(nms, nms)
    v
}

.async_check_h <- function(h, what = "h") {
    if (!is.numeric(h) || length(h) != 1L || is.na(h) || h <= 0) {
        stop(sprintf("'%s' must be a single positive number", what), call. = FALSE)
    }
    invisible(TRUE)
}

.async_check_one_sided <- function(one_sided) {
    if (!is.logical(one_sided) || length(one_sided) != 1L || is.na(one_sided)) {
        stop("'one_sided' must be TRUE or FALSE", call. = FALSE)
    }
    invisible(TRUE)
}

# Collapse the in-window pairs onto covariate rows.  Everything downstream --
# the fit, the sandwich, the cross-validation loss -- is a function of these
# quantities alone, so the pair design is formed once and discarded.
.async_collapse <- function(d, h, one_sided) {
    lo <- findInterval(seq_len(d$n) - 1L, d$xi)
    hi <- findInterval(seq_len(d$n), d$xi)
    supp <- if (one_sided) c(0, 1) else c(-1, 1)

    pr <- async_pairs_cpp(
        as.integer(d$yi - 1L), d$yt, d$xt,
        as.integer(lo), as.integer(hi), h, supp[1], supp[2]
    )
    if (!length(pr$j)) {
        return(NULL)
    }
    w <- epan_kernel(pr$z) / h
    # Strictly before, matching the survival side: a covariate observed at
    # exactly the response time is not "in the past".
    if (one_sided) w[pr$z <= 0] <- 0
    ac <- async_pair_accum_cpp(pr$j, pr$k, w, d$y, nrow(d$X))
    ac$npair <- if (one_sided) sum(pr$z > 0) else length(pr$j)
    ac$nresp <- length(unique(pr$j))
    ac
}

# Weighted squared error of the fitted mean, evaluated from the collapsed
# quantities:  sum_m w_m (Y_m - mu_k)^2 = sy2 - 2 sum_k c_k mu_k + sum_k w_k mu_k^2.
.async_loss <- function(ac, X, beta, lk) {
    mu <- .async_mu(as.numeric(X %*% beta), lk)
    num <- ac$sy2 - 2 * sum(ac$c * mu) + sum(ac$omega * mu^2)
    den <- sum(ac$omega)
    if (!is.finite(num) || !is.finite(den) || den <= 0) {
        return(NA_real_)
    }
    num / den
}

.async_fit_ti <- function(d, h, one_sided, lk, maxit, tol, start = NULL) {
    ac <- .async_collapse(d, h, one_sided)
    if (is.null(ac)) {
        return(NULL)
    }
    fit <- async_solve_cpp(
        d$X, as.integer(d$xi - 1L), ac$omega, ac$c, d$n,
        lk, as.integer(maxit), tol,
        if (is.null(start)) numeric(d$p) else start
    )
    fit$npair <- ac$npair
    fit$nresp <- ac$nresp
    fit
}

# A bandwidth that admits almost nothing is nearly always a units mistake --
# times in days with h chosen as if they were on the unit interval.
.async_warn_sparse <- function(nresp, ny, h) {
    if (nresp < 0.05 * ny) {
        warning(sprintf(
            paste0(
                "only %d of %d response occasions have a covariate observation ",
                "within h = %s. Check that 'h' is on the same scale as the ",
                "observation times."
            ),
            nresp, ny, format(h)
        ), call. = FALSE)
    }
    invisible(NULL)
}


#' Asynchronous longitudinal regression with time-invariant coefficients
#'
#' @description
#' Regress a sparsely observed longitudinal response on a sparsely observed
#' longitudinal covariate when the two are measured on **different** time grids
#' and are never seen together.
#'
#' Each response occasion \eqn{T_{ij}} is paired with every covariate occasion
#' \eqn{S_{ik}} of the same subject, and each pair contributes in proportion to
#' its time separation, so no pair is discarded and no value is carried forward.
#' The estimator solves the kernel-weighted estimating equation of Cao, Zeng and
#' Fine (2015),
#' \deqn{U_n(\beta) = n^{-1} \sum_i \sum_j \sum_k W_h(T_{ij} - S_{ik})
#'       X_i(S_{ik}) [ Y_i(T_{ij}) - g\{X_i(S_{ik})^\top \beta\} ] = 0,}
#' with \eqn{W_h(u) = W(u/h)/h} and \eqn{W} the Epanechnikov kernel.
#'
#' @details
#' # Why not something simpler
#'
#' Two obvious shortcuts are inconsistent here, which is the reason this
#' function exists.
#'
#' *Last value carried forward* replaces \eqn{X_i(T_{ij})} by the most recent
#' observed value. Under sparse sampling the gap back to that value does not
#' shrink as \eqn{n} grows, so the substituted covariate keeps a non-vanishing
#' error and the coefficient is attenuated toward zero.
#'
#' *Regression calibration* smooths the covariate path first and substitutes the
#' smoothed value. Under sparse sampling the smoothing window never empties, so
#' the same attenuation appears; for a nonlinear link there is a Jensen bias as
#' well, because the estimating equation needs the average of \eqn{g(\cdot)}
#' rather than \eqn{g} of the average.
#'
#' The kernel-weighted equation instead evaluates the link at each **observed**
#' covariate vector, and lets the weight express how much that observation says
#' about the response occasion.
#'
#' # Rate and bandwidth
#'
#' The estimator is consistent and asymptotically normal at the smoothing rate
#' \eqn{(nh)^{1/2}}, not at \eqn{\sqrt n}. The bandwidth therefore matters: too
#' small and few pairs contribute, too large and the bias grows. Use
#' [kee_async_cv()] rather than guessing.
#'
#' `h` is on the **same scale as the observation times**. Unlike the survival
#' estimators there is no requirement that times lie on \eqn{[0, 1]}, but `h`
#' must match whatever scale is used. A warning is issued when fewer than 5% of
#' response occasions have any covariate observation in their window, which is
#' the signature of a units mistake.
#'
#' # Half kernel or full kernel
#'
#' `one_sided = FALSE`, the default, is the full two-sided kernel of the paper.
#' `one_sided = TRUE` admits only covariate observations strictly **before** the
#' response, which is what a causal reading of the covariate path requires and
#' what the survival estimators in this package do by default. It roughly halves
#' the number of contributing pairs, so it needs a larger `h` to reach the same
#' effective sample size.
#'
#' # Assumptions worth checking
#'
#' Cao, Zeng and Fine assume the covariance function of the covariate process is
#' twice differentiable, which gives a bias of order \eqn{h^2}. Their Section 6
#' notes that relaxing this to one-sided differentiability, which admits
#' processes with independent increments such as the Ornstein-Uhlenbeck process,
#' leaves a bias of order \eqn{h} instead. The practical symptom is an estimate
#' that drifts steadily as `h` grows, so fit at several bandwidths and look
#' before trusting one number.
#'
#' Observation times are assumed independent of the response and covariate
#' processes. Informative observation times need a different estimator.
#'
#' @param data_y Data frame of response occasions: one row per time the
#'   response was recorded. First argument, so the function pipes.
#' @param data_x Data frame of covariate occasions: one row per time the
#'   covariates were recorded. Its time grid need not overlap `data_y`'s at all.
#' @param formula A two-sided formula, `y ~ x1 + x2`. **The left-hand side is
#'   looked up in `data_y` and the right-hand side in `data_x`**, which is the
#'   one unusual thing about this interface and follows from the data being in
#'   two tables: there is no single frame holding both.
#' @param id Subject identifier naming a column present in **both** tables.
#'   Give it unquoted (`id = subject`), as a string (`id = "subject"`), or
#'   embraced (`id = {{ col }}`) when calling from inside another function.
#'   Non-numeric identifiers are allowed.
#' @param time Observation time, naming a column present in **both** tables.
#'   Accepts the same three forms as `id`.
#' @param h Positive bandwidth, on the same scale as `time`. If omitted, a
#'   rule-of-thumb value is read off the observation times and reported in a
#'   message. That is a starting point, not a tuned choice: use
#'   [kee_async_cv()] to select it from the data.
#' @param one_sided Logical. `FALSE` (default) is the full two-sided kernel of
#'   the paper; `TRUE` is the half kernel, admitting only covariate observations
#'   strictly before the response.
#' @param link One of `"identity"` (closed form), `"log"` or `"logistic"`
#'   (Newton-Raphson).
#' @param intercept Logical, include an intercept.
#' @param maxit,tol Newton-Raphson controls, ignored when `link = "identity"`.
#'
#' @return An object of class `kee` with components
#' \describe{
#'   \item{coefficients}{Named coefficient vector.}
#'   \item{var}{Sandwich covariance matrix \eqn{A^{-1} B A^{-1}}.}
#'   \item{A, Sigma}{The bread and the meat of the sandwich.}
#'   \item{n}{Number of subjects, the unit the asymptotics are in.}
#'   \item{npair}{Number of weighted (response, covariate) pairs contributing.}
#'   \item{h, one_sided, link}{Fit settings.}
#'   \item{convergence}{`0` converged, `1` singular, `2` hit `maxit`.}
#'   \item{call}{The matched call.}
#' }
#' `coef()`, `vcov()`, `confint()`, `nobs()`, `print()` and `summary()` methods
#' are available.
#'
#' @references
#' Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
#' asynchronous longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B} 77, 755-776.
#'
#' @seealso [kee_async_cv()] to choose `h`, [kee_async_td()] for a coefficient
#'   curve \eqn{\beta(t)}, [sim_async_data()] to simulate.
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 150, beta = c(0.5, 1.5))
#'
#' # The two tables are on different time grids and share only `id`.
#' head(d$y, 3)
#' head(d$x, 3)
#'
#' fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.25)
#' summary(fit)
#' confint(fit)
#' tidy(fit, conf.int = TRUE)
#'
#' # Data comes first, so it pipes.
#' d$y |>
#'   kee_async(d$x, y ~ x, id = id, time = time, h = 0.25) |>
#'   glance()
#'
#' # Half kernel: only covariate observations preceding the response.
#' kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.25, one_sided = TRUE)
#' @export
kee_async <- function(data_y, data_x, formula, id, time, h = NULL, one_sided = FALSE,
                      link = c("identity", "log", "logistic"),
                      intercept = TRUE, maxit = 50L, tol = 1e-8) {
    link <- match.arg(link)
    .async_check_one_sided(one_sided)
    d <- .async_parse(
        formula, data_y, data_x,
        .async_name(rlang::enquo(id)), .async_name(rlang::enquo(time)), intercept
    )
    if (is.null(h)) {
        h <- default_bandwidth_async(c(d$yt, d$xt), d$n)
        announce_bandwidth(h, "kee_async_cv")
    }
    .async_check_h(h)

    fit <- .async_fit_ti(d, h, one_sided, .async_link_code(link), maxit, tol)
    if (is.null(fit)) {
        stop(
            "no response occasion has a covariate observation inside the ",
            "window; widen 'h' or check that it is on the same scale as 'time'",
            call. = FALSE
        )
    }
    if (fit$convergence == 1L) {
        stop(
            "the weighted bread matrix is singular; widen 'h' or drop ",
            "collinear covariates",
            call. = FALSE
        )
    }
    if (fit$convergence == 2L) {
        warning("Newton-Raphson reached 'maxit' without meeting 'tol'", call. = FALSE)
    }
    .async_warn_sparse(fit$nresp, length(d$y), h)

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
        npair = fit$npair,
        h = h,
        one_sided = one_sided,
        link = link,
        model = sprintf(
            "Asynchronous longitudinal regression (%s link, %s kernel)",
            link, if (one_sided) "half" else "full"
        ),
        call = match.call()
    )
    # Tagged ahead of "kee" so augment() can dispatch on it: a fitted mean per
    # covariate occasion is well defined here and is not for a hazard fit.
    class(out) <- c("kee_async", "kee")
    out
}


#' Choose the bandwidth for an asynchronous longitudinal fit
#'
#' @description
#' K-fold cross-validation over **subjects** for [kee_async()]. For each
#' candidate bandwidth the estimator is fitted on the training subjects and
#' scored on the held-out subjects by the kernel-weighted squared error of the
#' fitted mean,
#' \deqn{ \frac{\sum_{i \in \mathrm{test}} \sum_j \sum_k W_h(T_{ij} - S_{ik})
#'        [Y_i(T_{ij}) - g\{X_i(S_{ik})^\top \hat\beta\}]^2}
#'       {\sum_{i \in \mathrm{test}} \sum_j \sum_k W_h(T_{ij} - S_{ik})}.}
#'
#' @details
#' # How the folds are formed
#'
#' Folds are formed over subjects, never over rows. Rows within a subject come
#' from one trajectory, so splitting by row would put a subject on both sides of
#' the split and leak the answer into the held-out score.
#'
#' # Why the loss is normalised
#'
#' Dividing by the total weight is what makes bandwidths comparable. A larger
#' `h` admits more pairs and so accumulates a larger unnormalised error whatever
#' the quality of the fit; without the denominator the criterion would select
#' the smallest candidate every time.
#'
#' # Why squared error and not a likelihood
#'
#' The estimator solves an estimating equation, not a likelihood, so there is no
#' held-out log-likelihood to evaluate. The same weighted squared error on the
#' response scale is used for all three links, applied after `link`.
#'
#' # The default grid
#'
#' When `h_grid` is `NULL` the grid is log-spaced over
#' \eqn{[2 (Q_3 - Q_1) n^{-0.7},\, 2 (Q_3 - Q_1) n^{-0.3}]}, where
#' \eqn{Q_1, Q_3} are quartiles of the pooled observation times of both tables
#' and \eqn{n} is the number of subjects. This is the rule used by the reference
#' implementation of Cao, Zeng and Fine's method, and it adapts to the scale of
#' `time` on its own.
#'
#' Always look at `cv_results`. If the minimum sits at an endpoint of the grid a
#' warning is issued, and the grid should be widened by hand: the selected value
#' is then a boundary artefact rather than a minimum.
#'
#' @inheritParams kee_async
#' @param h_grid Optional numeric vector of candidate bandwidths. When `NULL` a
#'   grid is generated from the observation times; see Details.
#' @param n_h Number of candidates to generate when `h_grid` is `NULL`.
#' @param K Number of folds. Reduced to the number of subjects if larger.
#' @param seed Optional integer seed for the subject-to-fold assignment. Pass it
#'   (or call [set.seed()] first) to make the selection reproducible.
#' @param quiet Logical; suppress progress output.
#'
#' @return An object of class `cv.kee_async`:
#' \describe{
#'   \item{h_cv}{The selected bandwidth.}
#'   \item{fit}{The [kee_async()] fit refitted on all subjects at `h_cv`.}
#'   \item{cv_results}{Data frame of `h`, the mean held-out loss `cvloss`, and
#'     `nfold_used`, the number of folds that produced a usable fit.}
#'   \item{h_grid, fold_id, seed, call}{Selection metadata.}
#' }
#'
#' @references
#' Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
#' asynchronous longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B} 77, 755-776.
#'
#' @seealso [kee_async()]
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 150, beta = c(0.5, 1.5))
#' cv <- d$y |>
#'   kee_async_cv(d$x, y ~ x,
#'     id = id, time = time,
#'     h_grid = c(0.15, 0.25, 0.40), K = 3, seed = 1, quiet = TRUE
#'   )
#' cv
#' cv$cv_results
#' @export
kee_async_cv <- function(data_y, data_x, formula, id, time, h_grid = NULL,
                         n_h = 10L, K = 5L, one_sided = FALSE,
                         link = c("identity", "log", "logistic"),
                         intercept = TRUE, maxit = 50L, tol = 1e-8,
                         seed = NULL, quiet = FALSE) {
    link <- match.arg(link)
    .async_check_one_sided(one_sided)
    call <- match.call()
    d <- .async_parse(
        formula, data_y, data_x,
        .async_name(rlang::enquo(id)), .async_name(rlang::enquo(time)), intercept
    )
    lk <- .async_link_code(link)

    if (is.null(h_grid)) {
        pooled <- c(d$yt, d$xt)
        iqr <- diff(stats::quantile(pooled, c(0.25, 0.75), names = FALSE))
        if (!is.finite(iqr) || iqr <= 0) iqr <- diff(range(pooled))
        if (!is.finite(iqr) || iqr <= 0) {
            stop("observation times are constant; supply 'h_grid'", call. = FALSE)
        }
        h_grid <- exp(seq(
            log(2 * iqr * d$n^(-0.7)), log(2 * iqr * d$n^(-0.3)),
            length.out = max(2L, as.integer(n_h))
        ))
    }
    h_grid <- sort(unique(as.numeric(h_grid)))
    if (!length(h_grid) || anyNA(h_grid) || any(h_grid <= 0)) {
        stop("'h_grid' must be positive and free of NA", call. = FALSE)
    }

    K <- as.integer(K)
    if (K > d$n) {
        if (!quiet) {
            message(
                "Requested K = ", K, " exceeds n = ", d$n,
                " subjects; using K = ", d$n, " (leave-one-out)."
            )
        }
        K <- d$n
    }
    if (K < 2L) stop("'K' must be at least 2", call. = FALSE)

    if (!is.null(seed)) set.seed(seed)
    fold_id <- sample(rep_len(seq_len(K), d$n))

    losses <- matrix(NA_real_, K, length(h_grid))
    for (m in seq_along(h_grid)) {
        if (!quiet) cat(sprintf("Evaluating bandwidth h = %g\n", h_grid[m]))
        for (f in seq_len(K)) {
            tr <- which(fold_id != f)
            te <- which(fold_id == f)
            if (!length(tr) || !length(te)) next
            fit <- .async_fit_ti(
                .async_subset(d, tr), h_grid[m], one_sided, lk, maxit, tol
            )
            if (is.null(fit) || fit$convergence == 1L) next
            if (!all(is.finite(fit$beta))) next
            dte <- .async_subset(d, te)
            ac <- .async_collapse(dte, h_grid[m], one_sided)
            if (is.null(ac)) next
            losses[f, m] <- .async_loss(ac, dte$X, as.numeric(fit$beta), lk)
        }
    }

    cvloss <- colMeans(losses, na.rm = TRUE)
    nfold_used <- colSums(!is.na(losses))
    cvloss[nfold_used == 0L] <- NA_real_
    if (all(is.na(cvloss))) {
        stop(
            "no candidate bandwidth produced a usable fit on any fold; widen ",
            "'h_grid' or check that it is on the same scale as 'time'",
            call. = FALSE
        )
    }
    best_h <- h_grid[which.min(cvloss)]
    if (!quiet) cat(sprintf("Selected h = %g\n", best_h))
    if (length(h_grid) > 1L && best_h %in% range(h_grid)) {
        warning(
            "the selected bandwidth is at an endpoint of 'h_grid'; widen the ",
            "grid to check that the minimum is interior",
            call. = FALSE
        )
    }

    fit_call <- call
    fit_call[[1L]] <- quote(skmle::kee_async)
    fit_call$h <- best_h
    fit_call$h_grid <- NULL
    fit_call$n_h <- NULL
    fit_call$K <- NULL
    fit_call$seed <- NULL
    fit_call$quiet <- NULL
    fit <- eval(fit_call, parent.frame())

    out <- list(
        h_cv = best_h,
        fit = fit,
        cv_results = tibble::tibble(
            h = h_grid, cvloss = cvloss, nfold_used = nfold_used
        ),
        h_grid = h_grid,
        fold_id = fold_id,
        seed = seed,
        call = call
    )
    class(out) <- "cv.kee_async"
    out
}

#' @param x A `cv.kee_async` object.
#' @param ... Ignored.
#' @rdname kee_async_cv
#' @export
print.cv.kee_async <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n", length(unique(x$fold_id)),
    "-fold subject-level cross-validation\n\n",
    sep = ""
  )
  print(x$cv_results, row.names = FALSE, digits = 5)
  cat("\nSelected h = ", format(x$h_cv), "\n", sep = "")
  cat("\nCoefficients at the refit:\n")
  print(x$fit$coefficients)
  invisible(x)
}


#' Asynchronous longitudinal regression with time-dependent coefficients
#'
#' @description
#' Fit \eqn{E\{Y(t) \mid X(t)\} = g\{X(t)^\top \beta(t)\}} from asynchronous
#' sparse longitudinal data, estimating the coefficient curve \eqn{\beta(t)}
#' pointwise. This is the time-dependent coefficient estimator of Cao, Zeng and
#' Fine (2015), the companion to [kee_async()].
#'
#' At a target time \eqn{t} a response occasion and a covariate occasion are
#' weighted by their **separate** distances to \eqn{t},
#' \deqn{U_n\{\beta(t)\} = n^{-1} \sum_i \sum_j \sum_k W_{h_1}(t - T_{ij})
#'       W_{h_2}(t - S_{ik}) X_i(S_{ik})
#'       [ Y_i(T_{ij}) - g\{X_i(S_{ik})^\top \beta(t)\} ] = 0.}
#' The lag between the two occasions never enters, only their distances to
#' \eqn{t}. That is the difference from [kee_async()], where the lag is
#' everything.
#'
#' @details
#' # A slower rate than you may expect
#'
#' Because both time arguments are smoothed, the estimator converges at the
#' **bivariate** rate \eqn{(n h_1 h_2)^{1/2}}. That is slower than the
#' \eqn{(nh)^{1/2}} of the time-invariant fit, and much slower than the usual
#' varying-coefficient rate available when response and covariate are observed
#' together. Expect wide bands. Do not read a wiggle in the curve as structure
#' without checking it against them.
#'
#' A practical consequence: a sample that gives a comfortable time-invariant fit
#' can be far too small for a credible curve. If the bands cover a horizontal
#' line across the whole range, the honest reading is that the data do not
#' resolve time variation, not that \eqn{\beta(t)} is flat.
#'
#' # Reading the output
#'
#' `se` holds pointwise sandwich standard errors at each target time. They are
#' pointwise, not simultaneous, so a curve leaving the band at one or two target
#' times is not evidence of anything on its own.
#'
#' No undersmoothing correction is applied, so the intervals are centred on the
#' smoothed curve \eqn{E\hat\beta(t)} rather than on \eqn{\beta(t)}. With a large
#' bandwidth the curve is flattened toward a constant and the bands do not widen
#' to say so. Refit at a smaller `h` to see how much of the shape is smoothing.
#'
#' Target times within a bandwidth of the edge of the observed time range draw on
#' a one-sided window and are the least reliable part of the curve. Keep `times`
#' inside the data.
#'
#' # Cost
#'
#' The product weight factorises over the two occasion indices, so no pair is
#' ever enumerated: the response side collapses to two per-subject scalars and
#' the covariate side to one weight per row. Each target time costs
#' \eqn{O(n_y + n_x p^2)}. Fits at successive target times are warm-started from
#' the previous solution when the link is nonlinear.
#'
#' @inheritParams kee_async
#' @param times Numeric vector of target times at which to estimate
#'   \eqn{\beta(t)}. Defaults to 25 points spanning the 10th to 90th percentile
#'   of the observed response times, which keeps them inside the data: target
#'   times near the ends of the range are fitted from a one-sided window and are
#'   the least reliable part of the curve.
#' @param h Bandwidth for the response side, \eqn{h_1}. If omitted, a
#'   rule-of-thumb value is used and reported in a message; see [kee_async()].
#' @param h2 Bandwidth for the covariate side, \eqn{h_2}. Defaults to `h`.
#' @param one_sided Logical. `FALSE` (default) is the full two-sided kernel.
#'   `TRUE` restricts the **covariate** side to observations strictly before the
#'   target time, so \eqn{\beta(t)} uses only covariate values already available
#'   at `t`. The response side stays two-sided: it is a local average around
#'   `t`, not a filtering step. Both lags are measured as "how far in the past",
#'   matching [kee_async()] and the survival estimators.
#'
#' @return An object of class `kee_td`:
#' \describe{
#'   \item{coefficients}{`length(times)` by `p` matrix of \eqn{\hat\beta(t)}.}
#'   \item{se, z, p.value}{Pointwise standard errors and Wald statistics, the
#'     same shape as `coefficients`.}
#'   \item{var}{List of \eqn{p \times p} sandwich matrices, one per target time.}
#'   \item{times, h, h2, one_sided, link, n}{Fit settings.}
#'   \item{nactive}{Covariate rows carrying nonzero weight at each target time.}
#'   \item{convergence}{Per target time: `0` converged, `1` singular, `2` hit
#'     `maxit`, `3` no data in the window.}
#' }
#' `coef()`, `vcov()`, `confint()`, `nobs()`, `print()` and `plot()` methods are
#' available.
#'
#' @references
#' Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
#' asynchronous longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B} 77, 755-776.
#'
#' @seealso [kee_async()], [plot.kee_td()]
#'
#' @examples
#' set.seed(1)
#' d <- sim_async_data(n = 200, beta = function(tt) cbind(0.5, 1 + tt))
#' fit <- d$y |>
#'   kee_async_td(d$x, y ~ x,
#'     id = id, time = time,
#'     times = seq(0.2, 0.8, by = 0.1), h = 0.3
#'   )
#' coef(fit)
#' tidy(fit)
#' @export
kee_async_td <- function(data_y, data_x, formula, id, time, times = NULL,
                         h = NULL, h2 = h, one_sided = FALSE,
                         link = c("identity", "log", "logistic"),
                         intercept = TRUE, maxit = 50L, tol = 1e-8) {
    link <- match.arg(link)
    .async_check_one_sided(one_sided)
    d <- .async_parse(
        formula, data_y, data_x,
        .async_name(rlang::enquo(id)), .async_name(rlang::enquo(time)), intercept
    )
    if (is.null(h)) {
        h <- default_bandwidth_async(c(d$yt, d$xt), d$n)
        announce_bandwidth(h, "kee_async_cv")
        if (missing(h2)) h2 <- h
    }
    .async_check_h(h)
    .async_check_h(h2, "h2")
    if (is.null(times)) {
        # Inside the data: the ends of the range are fitted from a one-sided
        # window and are the least reliable part of any curve.
        times <- seq(
            stats::quantile(d$yt, 0.1, names = FALSE),
            stats::quantile(d$yt, 0.9, names = FALSE),
            length.out = 25L
        )
    }
    if (!is.numeric(times) || !length(times) || anyNA(times)) {
        stop("'times' must be a non-empty numeric vector without NA", call. = FALSE)
    }

    lk <- .async_link_code(link)
    nt <- length(times)
    nms <- colnames(d$X)

    beta <- se <- matrix(NA_real_, nt, d$p, dimnames = list(NULL, nms))
    vars <- vector("list", nt)
    conv <- integer(nt)
    nactive <- integer(nt)
    start <- numeric(d$p)

    for (m in seq_len(nt)) {
        # Both lags are "how far in the past", the package-wide convention.
        wy <- epan_kernel((times[m] - d$yt) / h) / h
        lag_x <- times[m] - d$xt
        vx <- epan_kernel(lag_x / h2) / h2
        if (one_sided) vx[lag_x <= 0] <- 0

        ac <- async_td_accum_cpp(
            as.integer(d$yi - 1L), wy, d$y, as.integer(d$xi - 1L), vx, d$n
        )
        nactive[m] <- sum(ac$omega != 0 | ac$c != 0)
        if (!nactive[m]) {
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
        stop(
            "no target time has data in its window; widen 'h'/'h2', move ",
            "'times' inside the observed range, or check the time scale",
            call. = FALSE
        )
    }
    if (any(conv == 2L)) {
        warning("Newton-Raphson reached 'maxit' at ", sum(conv == 2L),
            " target time(s)",
            call. = FALSE
        )
    }
    if (any(conv == 3L)) {
        warning(sum(conv == 3L), " target time(s) had no data in their window",
            call. = FALSE
        )
    }

    z <- beta / se
    out <- list(
        coefficients = beta,
        se = se,
        z = z,
        p.value = 2 * stats::pnorm(-abs(z)),
        var = vars,
        times = times,
        h = h,
        h2 = h2,
        one_sided = one_sided,
        link = link,
        n = d$n,
        nobs = c(response = length(d$y), covariate = nrow(d$X)),
        nactive = nactive,
        convergence = conv,
        call = match.call()
    )
    class(out) <- "kee_td"
    out
}


#' @param x A `kee_td` object.
#' @param ... Ignored.
#' @rdname kee_async_td
#' @export
print.kee_td <- function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nTime-dependent coefficients, ", x$link, " link\n", sep = "")
    cat("subjects: ", x$n, "   bandwidths: h = ", format(x$h),
        ", h2 = ", format(x$h2),
        "   kernel: ", if (x$one_sided) "half" else "full", "\n\n",
        sep = ""
    )
    print(cbind(t = x$times, x$coefficients), digits = 4)
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
#' Draws \eqn{\hat\beta_j(t)} against `t` with pointwise Wald bands, one panel
#' per coefficient.
#'
#' The bands are pointwise, not simultaneous: reading them as a confidence
#' region for the whole curve overstates the evidence. They are also not
#' corrected for smoothing bias, so they cover \eqn{E\hat\beta(t)} rather than
#' \eqn{\beta(t)}; a large bandwidth flattens the curve without widening the
#' band to say so.
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
#' fit <- d$y |>
#'   kee_async_td(d$x, y ~ x,
#'     id = id, time = time,
#'     times = seq(0.2, 0.8, by = 0.05), h = 0.3
#'   )
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

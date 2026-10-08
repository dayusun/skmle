# Asynchrony-corrected kernel estimating equation and the synchronous-loss
# cross-validation that selects its bandwidth.
#
# The kernel equation evaluates the link at the covariate observed at S, while
# the response depends on the covariate at T. For a stationary Gaussian
# covariate with variogram gamma(u) = E{x(t) - x(t - u)}^2, mean m and variance
# s2, x(T) | x(S) is normal with mean m + rho (x(S) - m) and variance
# s2 (1 - rho^2), rho = 1 - gamma(u) / (2 s2). Replacing g(eta_S) by
# E[g(eta_T) | x(S)] removes the asynchrony bias at every lag, so the corrected
# fit is unbiased at any fixed h.

# Every within-subject (response, covariate) pair, as a matrix with the lag
# u = T - S, the response, the covariate, and the subject code. With
# one_sided = TRUE only pairs with u > 0 are kept, since no other pair is ever
# weighted.
.async_all_pairs <- function(d, xcol, one_sided) {
    ys <- split(seq_along(d$yi), factor(d$yi, levels = seq_len(d$n)))
    xs <- split(seq_along(d$xi), factor(d$xi, levels = seq_len(d$n)))
    xv <- d$X[, xcol]
    P <- do.call(rbind, lapply(seq_len(d$n), function(i) {
        yi <- ys[[i]]
        xi <- xs[[i]]
        if (!length(yi) || !length(xi)) return(NULL)
        cbind(
            u = as.vector(outer(d$yt[yi], d$xt[xi], "-")),
            y = rep(d$y[yi], times = length(xi)),
            x = rep(xv[xi], each = length(yi)),
            sid = i
        )
    }))
    if (is.null(P)) {
        stop("no subject has both a response and a covariate observation",
             call. = FALSE)
    }
    if (one_sided) P <- P[P[, "u"] > 0, , drop = FALSE]
    P
}

# Isotonic (nondecreasing) regression of {x(S_k) - x(S_k')}^2 on |S_k - S_k'|
# over all within-subject covariate pairs. It uses no outcome and no
# bandwidth. Returns the variogram as a step function (0 below the smallest
# observed lag), together with the covariate mean and variance.
.async_variogram <- function(d, xcol) {
    xv <- d$X[, xcol]
    cp <- do.call(rbind, lapply(split(seq_along(d$xi), d$xi), function(ix) {
        if (length(ix) < 2L) return(NULL)
        cb <- utils::combn(ix, 2L)
        cbind(abs(d$xt[cb[1L, ]] - d$xt[cb[2L, ]]), (xv[cb[1L, ]] - xv[cb[2L, ]])^2)
    }))
    if (is.null(cp)) {
        stop("the asynchrony correction needs at least one subject with two ",
             "covariate observations", call. = FALSE)
    }
    o <- order(cp[, 1L])
    lag <- cp[o, 1L]
    fit <- stats::isoreg(lag, cp[o, 2L])$yf
    m <- mean(xv)
    list(
        gamma = function(u) {
            i <- findInterval(u, lag)
            out <- fit[pmax(i, 1L)]
            out[i == 0L] <- 0
            out
        },
        m = m,
        s2 = mean((xv - m)^2)
    )
}

# Gauss-Hermite nodes and weights for E f(Z), Z ~ N(0, 1) (Golub-Welsch).
.gh_normal <- function(k = 30L) {
    b <- sqrt(seq_len(k - 1L) / 2)
    J <- matrix(0, k, k)
    J[cbind(seq_len(k - 1L), 2:k)] <- b
    J[cbind(2:k, seq_len(k - 1L))] <- b
    e <- eigen(J, symmetric = TRUE)
    list(z = e$values * sqrt(2), w = e$vectors[1L, ]^2)
}

# Solve sum_p W_h(u_p) X_p {y_p - E[g(eta_T) | x_S]} = 0 by Newton-Raphson for
# a model with one covariate (column xcol of the design) and an optional
# intercept. P is a pair matrix from .async_all_pairs(). Returns the estimate,
# the bread, the per-subject scores (one row per subject code 1..n) and a
# convergence code as in async_solve_cpp: 0 converged, 1 singular, 2 maxit.
.async_fit_corr <- function(P, h, n, p, xcol, lk, vg, start, maxit, tol,
                            one_sided) {
    w <- epan_kernel(P[, "u"] / h) / h
    if (one_sided) w[P[, "u"] <= 0] <- 0
    k <- w > 0
    if (!any(k)) return(NULL)
    P <- P[k, , drop = FALSE]
    w <- w[k]
    X <- matrix(0, nrow(P), p)
    if (p == 2L) X[, 1L] <- 1
    X[, xcol] <- P[, "x"]
    y <- P[, "y"]
    rho <- pmin(pmax(1 - vg$gamma(abs(P[, "u"])) / (2 * vg$s2), -1), 1)
    cm <- vg$m + rho * (P[, "x"] - vg$m)
    v <- vg$s2 * (1 - rho^2)
    # E[g(eta_T) | x_S] and its derivative in beta; X_T replaces x by cm + sd Z.
    Zt <- X
    Zt[, xcol] <- cm
    gh <- if (lk == 2L) .gh_normal() else NULL
    ef <- function(b) {
        bx <- b[xcol]
        if (lk == 0L) {
            mu <- as.vector(Zt %*% b)
            D <- Zt
        } else if (lk == 1L) {
            mu <- exp(as.vector(Zt %*% b) + bx^2 * v / 2)
            D <- Zt * mu
            D[, xcol] <- mu * (cm + bx * v)
        } else {
            mu <- 0
            D <- 0 * X
            sd <- sqrt(v)
            for (q in seq_along(gh$z)) {
                Zq <- Zt
                Zq[, xcol] <- cm + sd * gh$z[q]
                pq <- stats::plogis(as.vector(Zq %*% b))
                mu <- mu + gh$w[q] * pq
                D <- D + gh$w[q] * pq * (1 - pq) * Zq
            }
        }
        list(U = w * X * (y - mu), J = crossprod(X * w, D))
    }
    b <- if (is.null(start)) numeric(p) else start
    conv <- 2L
    for (it in seq_len(maxit)) {
        e <- ef(b)
        st <- tryCatch(solve(e$J, colSums(e$U)), error = function(err) NULL)
        if (is.null(st)) return(list(beta = b, convergence = 1L))
        b <- b + st
        if (!all(is.finite(b))) return(list(beta = b, convergence = 1L))
        if (max(abs(st)) < tol) {
            conv <- 0L
            break
        }
    }
    e <- ef(b)
    U <- rowsum(e$U, factor(P[, "sid"], levels = seq_len(n)), reorder = TRUE)
    U[is.na(U)] <- 0
    list(beta = b, bread = e$J, U = unname(U), convergence = conv,
         npair = nrow(P), nresp = NA_integer_)
}

# Held-out loss whose minimiser over h is the fitted mean's distance from the
# truth: squared error (identity), Poisson deviance up to terms free of beta
# (log), Bernoulli log loss (logistic).
.async_pair_loss <- function(y, eta, lk) {
    if (lk == 0L) {
        (y - eta)^2
    } else if (lk == 1L) {
        exp(eta) - y * eta
    } else {
        log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta
    }
}

# The smallest h at which, in every training fold of every split, at least
# ceiling(sqrt(n_train)) subjects have a pair inside the kernel. Below it a
# fold fit rests on a handful of pairs, its error does not vanish, and the
# held-out intercept is meaningless.
.async_grid_floor <- function(P, n, folds) {
    mi <- tapply(abs(P[, "u"]), factor(P[, "sid"], levels = seq_len(n)), min)
    mi[is.na(mi)] <- Inf
    if (is.null(folds)) {
        # Whole sample: ceiling(sqrt(n)) subjects.
        m <- sort(mi)
        return(m[min(ceiling(sqrt(n)), length(m))] * (1 + 1e-8))
    }
    fl <- 0
    for (r in seq_len(ncol(folds))) {
        for (f in sort(unique(folds[, r]))) {
            tr <- folds[, r] != f
            m <- sort(mi[tr])
            k <- min(ceiling(sqrt(sum(tr))), length(m))
            if (is.finite(m[k])) fl <- max(fl, m[k])
        }
    }
    fl * (1 + 1e-8)
}

# Saturated loss min_eta loss(y, eta): 0 for squared error and log loss,
# y - y log y for the Poisson loss. The held-out intercept estimates an
# expected loss, which cannot be below the mean saturated loss; a candidate
# whose intercept falls below it on some fold comes from a degenerate fold fit
# (a handful of pairs), and is dropped.
.async_sat_loss <- function(y, lk) {
    if (lk == 1L) ifelse(y > 0, y - y * log(y), 0) else numeric(length(y))
}

# Synchronous-loss cross-validation. For each split, fold and candidate h:
# fit on the training subjects (`fitter`), compute the unweighted loss of every
# held-out pair at its covariate value, regress the losses on
# (1, gdes(|u|)) -- the estimated variogram (matrix) at the pair's lag -- and
# keep the intercept, the held-out loss at zero asynchrony. Intercepts are
# averaged over folds and then over splits. Fits run from the largest h down,
# each starting from the previous estimate.
.async_cv_sync <- function(P, design_of, gdes, fitter, h_grid, folds, lk,
                           start, quiet) {
    G <- length(h_grid)
    R <- ncol(folds)
    nfold_used <- integer(G)
    Lr <- matrix(NA_real_, R, G)
    for (r in seq_len(R)) {
        Ks <- sort(unique(folds[, r]))
        losses <- matrix(NA_real_, length(Ks), G)
        for (fi in seq_along(Ks)) {
            te <- which(folds[, r] == Ks[fi])
            H <- P[P[, "sid"] %in% te, , drop = FALSE]
            if (!nrow(H) || length(te) == nrow(folds)) next
            Hx <- design_of(H)
            QR <- qr(cbind(1, gdes(abs(H[, "u"]))))
            sat <- mean(.async_sat_loss(H[, "y"], lk))
            st <- start
            for (m in rev(seq_len(G))) {
                if (!quiet && r == 1L && fi == 1L) {
                    cat(sprintf("Evaluating bandwidth h = %g\n", h_grid[m]))
                }
                fit <- fitter(te, h_grid[m], st)
                if (is.null(fit) || fit$convergence != 0L ||
                    !all(is.finite(fit$beta))) next
                st <- fit$beta
                ell <- .async_pair_loss(H[, "y"], as.vector(Hx %*% fit$beta), lk)
                ic <- qr.coef(QR, ell)[1L]
                if (is.finite(ic) && ic >= sat) losses[fi, m] <- ic
            }
        }
        nfold_used <- nfold_used + colSums(!is.na(losses))
        # A candidate counts only when every fold produced a loss: the mean of
        # a subset of fold intercepts is not comparable across candidates.
        Lr[r, ] <- colMeans(losses)
    }
    cvloss <- colMeans(Lr)
    if (all(is.na(cvloss))) {
        stop(
            "no candidate bandwidth produced a usable fit on every fold; widen ",
            "'h_grid' or check that it is on the same scale as 'time'",
            call. = FALSE
        )
    }
    list(cvloss = cvloss, nfold_used = nfold_used,
         index = which.min(replace(cvloss, is.na(cvloss), Inf)),
         folds = folds)
}

# Minimum-variance bandwidth. Fit the corrected equation on all subjects at
# every candidate (largest h first, each fit starting from the previous
# estimate) and return tr{Sx V(h)}, V the sandwich. The corrected fit is
# unbiased at every h, so this is the MSE criterion; V-hat(h) converges at
# every grid point, so the minimiser converges to the minimiser of V.
.async_select_var <- function(fitter, Sx, h_grid, nms, start, quiet) {
    G <- length(h_grid)
    crit <- rep(NA_real_, G)
    st <- start
    for (m in rev(seq_len(G))) {
        if (!quiet) cat(sprintf("Evaluating bandwidth h = %g\n", h_grid[m]))
        fit <- fitter(h_grid[m], st)
        if (is.null(fit) || fit$convergence != 0L || !all(is.finite(fit$beta))) next
        st <- fit$beta
        V <- .async_sandwich(fit$bread, fit$U, nms)
        if (all(is.finite(V))) crit[m] <- sum(Sx * V)
    }
    if (all(is.na(crit))) {
        stop(
            "no candidate bandwidth produced a usable fit; widen 'h_grid' or ",
            "check that it is on the same scale as 'time'",
            call. = FALSE
        )
    }
    list(crit = crit, index = which.min(replace(crit, is.na(crit), Inf)))
}

# The single covariate the Gaussian correction is defined for; errors otherwise.
.async_corr_col <- function(d) {
    xcol <- which(colnames(d$X) != "(Intercept)")
    if (length(xcol) != 1L) {
        stop(
            "the Gaussian correction is implemented for a single covariate; ",
            "use correct = \"pairs\"",
            call. = FALSE
        )
    }
    xcol
}

## ----------------------------------------------------------------------------
## Pair-moment estimator
##
## Each subject's own covariate pairs stand in for the unobserved x(T). For an
## ordered pair of covariate occasions (S_k, S_k'), (x(S_k), x(S_k')) has the
## law of (x(S), x(T)) at the same two times, so
##   U_h(b) = sum_{(T,S)} W_h(T - S) X(S) Y(T) / sum W_h(T - S)
##          - sum_{k != k'} W_h(S_k' - S_k) r(S_k') X(S_k) g{X(S_k')'b}
##                                      / sum W_h(S_k' - S_k) r(S_k')
## has mean zero at the true b for every h, with r = lambda_T / lambda_S the
## ratio of the response and covariate visit intensities. No model for the
## covariate process is involved.
## ----------------------------------------------------------------------------

# Response-covariate pairs with the covariate row index k (all covariates).
.async_ts_pairs <- function(d, one_sided) {
    ys <- split(seq_along(d$yi), factor(d$yi, levels = seq_len(d$n)))
    xs <- split(seq_along(d$xi), factor(d$xi, levels = seq_len(d$n)))
    P <- do.call(rbind, lapply(seq_len(d$n), function(i) {
        a <- ys[[i]]
        b <- xs[[i]]
        if (!length(a) || !length(b)) return(NULL)
        cbind(u = as.vector(outer(d$yt[a], d$xt[b], "-")),
              y = rep(d$y[a], times = length(b)),
              sid = i, k = rep(b, each = length(a)))
    }))
    if (is.null(P)) {
        stop("no subject has both a response and a covariate observation",
             call. = FALSE)
    }
    if (one_sided) P <- P[P[, "u"] > 0, , drop = FALSE]
    P
}

# Ordered within-subject covariate pairs: lag v = S_b - S_a, subject, rows a
# (the instrument, playing S) and b (the stand-in for T).
.async_cov_pairs <- function(d, one_sided) {
    Q <- do.call(rbind, lapply(split(seq_along(d$xi), d$xi), function(ix) {
        if (length(ix) < 2L) return(NULL)
        g <- expand.grid(a = ix, b = ix)
        g <- g[g$a != g$b, ]
        cbind(v = d$xt[g$b] - d$xt[g$a], sid = d$xi[g$a], a = g$a, b = g$b)
    }))
    if (is.null(Q)) {
        stop("correct = \"pairs\" needs at least one subject with two ",
             "covariate observations", call. = FALSE)
    }
    if (one_sided) Q <- Q[Q[, "v"] > 0, , drop = FALSE]
    Q
}

# lambda_T / lambda_S at each covariate occasion, from the observation times
# alone: logistic regression of "is a response occasion" on a natural spline of
# time with `df` degrees of freedom.
.async_intensity_ratio <- function(d, df = 4L) {
    tt <- c(d$yt, d$xt)
    lab <- rep(1:0, c(length(d$yt), length(d$xt)))
    if (length(unique(tt)) <= df + 1L) return(rep(length(d$yt) / length(d$xt), length(d$xt)))
    B <- splines::ns(tt, df = df)
    f <- stats::glm.fit(cbind(1, B), lab, family = stats::binomial())
    # When the two schedules barely overlap in time (a covariate collected in
    # one study phase only), the fit separates: an aliased spline column comes
    # back NA and the fitted probability reaches 0 or 1, so the ratio is NaN or
    # Inf and the pair fit fails on an NA condition. N occasions cannot resolve
    # a probability finer than 1/N, so bound it there and say so.
    cf <- f$coefficients
    cf[is.na(cf)] <- 0
    pr <- stats::plogis(as.vector(cbind(1, stats::predict(B, d$xt)) %*% cf))
    eps <- 0.5 / length(tt)
    nb <- sum(pr < eps | pr > 1 - eps)
    if (nb) {
        warning(nb, " covariate occasion(s) had an estimated response/covariate ",
                "intensity ratio at the edge of what the observation times can ",
                "resolve; bounded at ", signif(eps, 3), call. = FALSE)
        pr <- pmin(pmax(pr, eps), 1 - eps)
    }
    pr / (1 - pr)
}

# Columns that vary within some subject (excluding the intercept).
.async_time_varying <- function(d) {
    cols <- which(colnames(d$X) != "(Intercept)")
    tv <- cols[vapply(cols, function(j) {
        any(tapply(d$X[, j], d$xi, function(v) diff(range(v)) > 0))
    }, TRUE)]
    if (!length(tv)) {
        stop("no covariate varies within a subject, so there is no asynchrony ",
             "to correct", call. = FALSE)
    }
    tv
}

# vech of the isotonic variogram matrix Gamma(u) = E[dx dx'] of the
# time-varying covariates, by polarisation: isotonic variograms of x_j and of
# x_j + x_k, Gamma_jk = (gamma_{j+k} - gamma_jj - gamma_kk) / 2. Returns a
# function of the lag giving one column per entry. Used only by the CV.
.async_gamma_design <- function(d, tv) {
    cp <- do.call(rbind, lapply(split(seq_along(d$xi), d$xi), function(ix) {
        if (length(ix) < 2L) return(NULL)
        cb <- utils::combn(ix, 2L)
        cbind(abs(d$xt[cb[1L, ]] - d$xt[cb[2L, ]]),
              d$X[cb[1L, ], tv, drop = FALSE] - d$X[cb[2L, ], tv, drop = FALSE])
    }))
    if (is.null(cp)) {
        stop("the synchronous-loss criterion needs at least one subject with ",
             "two covariate observations", call. = FALSE)
    }
    lag <- cp[, 1L]
    D <- cp[, -1L, drop = FALSE]
    o <- order(lag)
    ls <- lag[o]
    iso <- function(z) stats::isoreg(ls, z[o])$yf
    q <- length(tv)
    dg <- lapply(seq_len(q), function(j) iso(D[, j]^2))
    cols <- list()
    for (j in seq_len(q)) {
        for (k in j:q) {
            cols[[length(cols) + 1L]] <- if (j == k) dg[[j]] else
                (iso((D[, j] + D[, k])^2) - dg[[j]] - dg[[k]]) / 2
        }
    }
    Y <- do.call(cbind, cols)
    function(u) {
        i <- findInterval(u, ls)
        out <- Y[pmax(i, 1L), , drop = FALSE]
        out[i == 0L, ] <- 0
        out
    }
}

# Solve the pair-moment equation by Newton-Raphson. P: response-covariate
# pairs (.async_ts_pairs), Q: covariate pairs (.async_cov_pairs), r: intensity
# ratio per covariate row. Returns the estimate, the bread, per-subject
# influence rows (one per subject code 1..n) and a convergence code.
.async_fit_pairs <- function(P, Q, h, d, lk, r, start, maxit, tol, one_sided) {
    w <- epan_kernel(P[, "u"] / h) / h
    if (one_sided) w[P[, "u"] <= 0] <- 0
    v <- epan_kernel(Q[, "v"] / h) / h * r[Q[, "b"]]
    if (one_sided) v[Q[, "v"] <= 0] <- 0
    kp <- w > 0
    kq <- v > 0
    if (!any(kp) || !any(kq)) return(NULL)
    P <- P[kp, , drop = FALSE]
    w <- w[kp]
    Q <- Q[kq, , drop = FALSE]
    v <- v[kq]
    n <- d$n
    full <- function(M, g) {
        out <- matrix(0, n, ncol(M))
        rr <- rowsum(M, g)
        out[as.integer(rownames(rr)), ] <- rr
        out
    }
    XS <- d$X[P[, "k"], , drop = FALSE]
    Xa <- d$X[Q[, "a"], , drop = FALSE]
    Xb <- d$X[Q[, "b"], , drop = FALSE]
    WA <- sum(w)
    WB <- sum(v)
    Ai <- full(w * XS * P[, "y"], P[, "sid"])
    wa <- full(cbind(w), P[, "sid"])[, 1L]
    wb <- full(cbind(v), Q[, "sid"])[, 1L]
    A <- colSums(Ai) / WA
    dmu <- function(eta) {
        if (lk == 0L) rep(1, length(eta)) else if (lk == 1L) exp(eta) else {
            pp <- stats::plogis(eta)
            pp * (1 - pp)
        }
    }
    ef <- function(b) {
        eta <- as.vector(Xb %*% b)
        Bi <- full(v * Xa * .async_mu(eta, lk), Q[, "sid"])
        list(Bi = Bi, U = A - colSums(Bi) / WB,
             J = crossprod(Xa * v * dmu(eta), Xb) / WB)
    }
    b <- if (is.null(start)) numeric(d$p) else start
    conv <- 2L
    for (it in seq_len(maxit)) {
        e <- ef(b)
        st <- tryCatch(solve(e$J, e$U), error = function(err) NULL)
        if (is.null(st)) return(list(beta = b, convergence = 1L))
        b <- b + st
        if (!all(is.finite(b))) return(list(beta = b, convergence = 1L))
        if (max(abs(st)) < tol) {
            conv <- 0L
            break
        }
    }
    e <- ef(b)
    B <- colSums(e$Bi) / WB
    U <- (Ai - outer(wa, A)) / WA - (e$Bi - outer(wb, B)) / WB
    list(beta = b, bread = e$J, U = unname(U), convergence = conv,
         npair = nrow(P), nresp = NA_integer_)
}

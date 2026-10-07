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

# Synchronous-loss cross-validation. For each split, fold and candidate h:
# fit the corrected equation on the training subjects, compute the unweighted
# loss of every held-out pair at its covariate value, regress the losses on
# (1, gamma-hat(|u|)) and keep the intercept, the held-out loss at zero
# asynchrony. Intercepts are averaged over folds and then over splits. Fits
# run from the largest h down, each starting from the previous estimate.
.async_cv_sync <- function(d, P, vg, xcol, h_grid, folds, one_sided, lk,
                           start, maxit, tol, quiet) {
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
            Ptr <- P[!P[, "sid"] %in% te, , drop = FALSE]
            if (!nrow(H) || !nrow(Ptr)) next
            Hx <- matrix(0, nrow(H), d$p)
            if (d$p == 2L) Hx[, 1L] <- 1
            Hx[, xcol] <- H[, "x"]
            QR <- qr(cbind(1, vg$gamma(abs(H[, "u"]))))
            st <- start
            for (m in rev(seq_len(G))) {
                if (!quiet && r == 1L && fi == 1L) {
                    cat(sprintf("Evaluating bandwidth h = %g\n", h_grid[m]))
                }
                fit <- .async_fit_corr(Ptr, h_grid[m], d$n, d$p, xcol, lk, vg,
                                       st, maxit, tol, one_sided)
                if (is.null(fit) || fit$convergence != 0L ||
                    !all(is.finite(fit$beta))) next
                st <- fit$beta
                ell <- .async_pair_loss(H[, "y"], as.vector(Hx %*% fit$beta), lk)
                losses[fi, m] <- qr.coef(QR, ell)[1L]
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

# The single covariate the correction is defined for; errors otherwise.
.async_corr_col <- function(d) {
    xcol <- which(colnames(d$X) != "(Intercept)")
    if (length(xcol) != 1L) {
        stop(
            "the asynchrony correction is implemented for a single covariate; ",
            "use correct = FALSE (kee_async) or method = \"plain\" (kee_async_cv)",
            call. = FALSE
        )
    }
    xcol
}

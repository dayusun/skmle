## The asynchrony-corrected fits (kee_async(correct = "pairs" / "gaussian")) and
## the two rules that select their bandwidth: the minimum-variance default
## (kee_async_cv(method = "variance")) and synchronous-loss cross-validation
## (method = "sync").

# Every within-subject pair of response and covariate, built by hand.
hand_pairs <- function(d) {
    do.call(rbind, lapply(sort(unique(d$y$id)), function(i) {
        yy <- d$y[d$y$id == i, ]
        xx <- d$x[d$x$id == i, ]
        if (!nrow(xx)) return(NULL)
        g <- expand.grid(j = seq_len(nrow(yy)), k = seq_len(nrow(xx)))
        data.frame(id = i, u = yy$time[g$j] - xx$time[g$k], y = yy$y[g$j],
                   x = xx$x[g$k])
    }))
}

# The isotonic variogram, built by hand.
hand_vario <- function(d) {
    cp <- do.call(rbind, lapply(split(d$x, d$x$id), function(xx) {
        if (nrow(xx) < 2L) return(NULL)
        cb <- utils::combn(nrow(xx), 2L)
        cbind(abs(xx$time[cb[1, ]] - xx$time[cb[2, ]]), (xx$x[cb[1, ]] - xx$x[cb[2, ]])^2)
    }))
    o <- order(cp[, 1])
    iso <- stats::isoreg(cp[o, 1], cp[o, 2])
    list(gamma = function(u) {
        i <- findInterval(u, cp[o, 1])
        ifelse(i == 0L, 0, iso$yf[pmax(i, 1L)])
    }, m = mean(d$x$x), s2 = mean((d$x$x - mean(d$x$x))^2))
}

test_that("the identity-link corrected fit solves the corrected equation", {
    set.seed(11)
    d <- sim_async_data(n = 80)
    h <- 0.3
    fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h, correct = "gaussian")

    P <- hand_pairs(d)
    vg <- hand_vario(d)
    w <- pmax(0.75 * (1 - (P$u / h)^2), 0) / h
    rho <- pmin(pmax(1 - vg$gamma(abs(P$u)) / (2 * vg$s2), -1), 1)
    X <- cbind(1, P$x)
    Z <- cbind(1, vg$m + rho * (P$x - vg$m))
    b <- solve(crossprod(X * w, Z), crossprod(X * w, P$y))
    expect_equal(unname(coef(fit)), as.vector(b), tolerance = 1e-8)
    expect_identical(fit$correct, "gaussian")
    expect_match(fit$model, "Gaussian correction")
    expect_equal(fit$variogram$gamma(c(0.05, 0.2, 0.6)), vg$gamma(c(0.05, 0.2, 0.6)))

    # Sandwich from per-subject scores of the same equation.
    r <- as.vector(P$y - Z %*% b)
    U <- rowsum(X * w * r, P$id)
    A <- crossprod(X * w, Z)
    V <- solve(A) %*% crossprod(U) %*% t(solve(A))
    expect_equal(unname(vcov(fit)), V, tolerance = 1e-8)
})

test_that("the correction removes the attenuation at a wide bandwidth", {
    skip_on_cran()
    set.seed(12)
    d <- sim_async_data(n = 600, beta = c(0.5, 1.5))
    plain <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.6)
    corr <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.6, correct = "gaussian")
    # The plain fit is attenuated by about beta * gamma-bar / (2 var x).
    expect_lt(coef(plain)[["x"]], 1.35)
    expect_lt(abs(coef(corr)[["x"]] - 1.5), 3 * sqrt(vcov(corr)[2, 2]))
})

test_that("the log and logistic corrections run and stay near the plain fit at a small h", {
    set.seed(13)
    for (lk in c("log", "logistic")) {
        d <- sim_async_data(n = 150, beta = c(0.3, 0.5), link = lk)
        p <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.05, link = lk)
        c1 <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.05, link = lk,
                        correct = "gaussian")
        expect_equal(c1$convergence, 0L)
        expect_lt(max(abs(coef(c1) - coef(p))), 0.1)
    }
})

test_that("the Gaussian correction refuses more than one covariate", {
    set.seed(14)
    d <- sim_async_data(n = 40, beta = c(0.5, 1, -1))
    expect_error(
        kee_async(d$y, d$x, y ~ x1 + x2, id = id, time = time, h = 0.3, correct = "gaussian"),
        "single covariate"
    )
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x1 + x2, id = id, time = time,
                     h_grid = c(0.2, 0.3), K = 2, quiet = TRUE, correct = "gaussian"),
        "single covariate"
    )
})

test_that("the synchronous-loss criterion is the held-out intercept on gamma-hat", {
    skip_on_cran()
    set.seed(15)
    d <- sim_async_data(n = 60)
    hs <- c(0.15, 0.4)
    K <- 2
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        h_grid = hs, K = K, seed = 3, quiet = TRUE,
                                        correct = "gaussian", method = "sync"))
    expect_identical(cv$method, "sync")

    set.seed(3)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    fold <- sample(rep_len(seq_len(K), length(subj)))
    P <- hand_pairs(d)
    vg <- hand_vario(d)
    # Fold fits use the variogram of all subjects: it involves no outcome.
    corr_fit <- function(Q, h) {
        w <- pmax(0.75 * (1 - (Q$u / h)^2), 0) / h
        rho <- pmin(pmax(1 - vg$gamma(abs(Q$u)) / (2 * vg$s2), -1), 1)
        X <- cbind(1, Q$x)
        Z <- cbind(1, vg$m + rho * (Q$x - vg$m))
        as.vector(solve(crossprod(X * w, Z), crossprod(X * w, Q$y)))
    }
    hand <- vapply(hs, function(h) mean(vapply(seq_len(K), function(f) {
        tr <- as.character(P$id) %in% subj[fold != f]
        b <- corr_fit(P[tr, ], h)
        H <- P[!tr, ]
        ell <- (H$y - b[1] - b[2] * H$x)^2
        coef(lm(ell ~ vg$gamma(abs(H$u))))[[1]]
    }, 1)), 1)
    expect_equal(cv$cv_results$cvloss, hand, tolerance = 1e-8)
    expect_equal(cv$h_cv, hs[which.min(cv$cv_results$cvloss)])
    # The returned fit is the corrected fit at the selected bandwidth.
    ref <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = cv$h_cv, correct = "gaussian")
    expect_equal(coef(cv$fit), coef(ref))
    expect_identical(cv$fit$correct, "gaussian")
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "synchronous loss")
})

test_that("the Gaussian default grid starts where every training fold has sqrt(n_train) subjects in the kernel", {
    skip_on_cran()
    set.seed(16)
    d <- sim_async_data(n = 100)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        n_h = 8, K = 4, seed = 2, quiet = TRUE,
                                        correct = "gaussian", method = "sync"))
    P <- hand_pairs(d)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    mi <- tapply(abs(P$u), factor(as.character(P$id), levels = subj), min)
    lo <- max(vapply(1:4, function(f) {
        tr <- cv$fold_id != f
        sort(mi[tr])[ceiling(sqrt(sum(tr)))]
    }, 1))
    expect_equal(cv$h_grid[1], lo, tolerance = 1e-6)
    expect_equal(max(cv$h_grid), max(abs(P$u)))
    expect_length(cv$h_grid, 8L)
    # At the floor, every training fold has that many subjects inside the kernel.
    for (f in 1:4) {
        tr <- subj[cv$fold_id != f]
        inside <- unique(P$id[abs(P$u) < cv$h_grid[1] & as.character(P$id) %in% tr])
        expect_gte(length(inside), ceiling(sqrt(length(tr))))
    }
})

## ------------------------------------------------------------ pair-moment fit

# Ordered within-subject covariate pairs, by hand.
hand_cov_pairs <- function(d) {
    do.call(rbind, lapply(split(d$x, d$x$id), function(xx) {
        if (nrow(xx) < 2L) return(NULL)
        g <- expand.grid(a = seq_len(nrow(xx)), b = seq_len(nrow(xx)))
        g <- g[g$a != g$b, ]
        data.frame(id = xx$id[1], v = xx$time[g$b] - xx$time[g$a], a = I(as.matrix(xx[g$a, -(1:2)])),
                   b = I(as.matrix(xx[g$b, -(1:2)])), tb = xx$time[g$b])
    }))
}

# lambda_T / lambda_S at time t: logistic regression of the occasion type on ns(time, 4).
hand_ratio <- function(d, t) {
    tt <- c(d$y$time, d$x$time)
    B <- splines::ns(tt, df = 4)
    f <- glm(rep(1:0, c(nrow(d$y), nrow(d$x))) ~ B, family = binomial())
    p <- plogis(cbind(1, predict(B, t)) %*% coef(f))
    as.vector(p / (1 - p))
}

test_that("the identity-link pair-moment fit solves the pair-moment equation", {
    set.seed(21)
    d <- sim_async_data(n = 120, beta = c(0.5, 1, -0.7))
    h <- 0.3
    fit <- kee_async(d$y, d$x, y ~ x1 + x2, id = id, time = time, h = h, correct = "pairs")
    expect_identical(fit$correct, "pairs")
    expect_match(fit$model, "pair-moment")

    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    # Response-covariate pairs.
    P <- do.call(rbind, lapply(sort(unique(d$y$id)), function(i) {
        yy <- d$y[d$y$id == i, ]
        xx <- d$x[d$x$id == i, ]
        g <- expand.grid(j = seq_len(nrow(yy)), k = seq_len(nrow(xx)))
        data.frame(id = i, u = yy$time[g$j] - xx$time[g$k], y = yy$y[g$j],
                   x1 = xx$x1[g$k], x2 = xx$x2[g$k])
    }))
    Q <- hand_cov_pairs(d)
    w <- W(P$u / h) / h
    v <- W(Q$v / h) / h * hand_ratio(d, Q$tb)
    XS <- cbind(1, P$x1, P$x2)
    Xa <- cbind(1, Q$a)
    Xb <- cbind(1, Q$b)
    A <- colSums(w * XS * P$y) / sum(w)
    M <- crossprod(Xa * v, Xb) / sum(v)
    b <- solve(M, A)
    expect_equal(unname(coef(fit)), as.vector(b), tolerance = 1e-8)

    # Sandwich from per-subject influence of the two ratio estimators.
    subj <- sort(unique(c(d$y$id, d$x$id)))
    sumby <- function(Z, id) {
        out <- matrix(0, length(subj), ncol(Z))
        r <- rowsum(Z, id)
        out[match(as.numeric(rownames(r)), subj), ] <- r
        out
    }
    Bv <- as.vector(Xb %*% b)
    psi <- (sumby(w * XS * P$y, P$id) - outer(sumby(cbind(w), P$id)[, 1], A)) / sum(w) -
        (sumby(v * Xa * Bv, Q$id) - outer(sumby(cbind(v), Q$id)[, 1], as.vector(M %*% b))) / sum(v)
    V <- solve(M) %*% crossprod(psi) %*% t(solve(M))
    expect_equal(unname(vcov(fit)), unname(V), tolerance = 1e-8)
})

test_that("the pair-moment fit takes several covariates and time-invariant ones", {
    set.seed(22)
    d <- sim_async_data(n = 150, beta = c(0.3, 0.6, -0.4), link = "logistic")
    d$x$g <- as.numeric(d$x$id %% 2)
    fit <- kee_async(d$y, d$x, y ~ x1 + x2 + g, id = id, time = time, h = 0.3,
                     link = "logistic", correct = "pairs")
    expect_equal(fit$convergence, 0L)
    expect_named(coef(fit), c("(Intercept)", "x1", "x2", "g"))
    expect_true(all(is.finite(sqrt(diag(vcov(fit))))))

    # A model with no covariate that varies within a subject has no asynchrony.
    expect_error(
        kee_async(d$y, d$x, y ~ g, id = id, time = time, h = 0.3, correct = "pairs"),
        "varies within a subject"
    )
})

test_that("the pair-moment correction removes the attenuation at a wide bandwidth", {
    skip_on_cran()
    set.seed(23)
    d <- sim_async_data(n = 1500, beta = c(0.5, 1.5))
    plain <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.6)
    pm <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.6, correct = "pairs")
    expect_lt(coef(plain)[["x"]], 1.35)
    expect_lt(abs(coef(pm)[["x"]] - 1.5), 3 * sqrt(vcov(pm)[2, 2]))
})

test_that("the saturated loss bounds the held-out intercept from below", {
    sat <- skmle:::.async_sat_loss
    expect_equal(sat(c(0, 1, 2), 0L), c(0, 0, 0))
    expect_equal(sat(c(0, 1), 2L), c(0, 0))
    expect_equal(sat(c(0, 1, 3), 1L), c(0, 1, 3 - 3 * log(3)))
})

test_that("sync CV selects h for the pair-moment fit by the held-out intercept", {
    skip_on_cran()
    set.seed(24)
    d <- sim_async_data(n = 80)
    hs <- c(0.15, 0.4)
    K <- 2
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        h_grid = hs, K = K, seed = 4, quiet = TRUE,
                                        method = "sync"))
    expect_identical(cv$correct, "pairs")

    set.seed(4)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    fold <- sample(rep_len(seq_len(K), length(subj)))
    P <- hand_pairs(d)
    Q <- hand_cov_pairs(d)
    rq <- hand_ratio(d, Q$tb)
    vg <- hand_vario(d)
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    hand <- vapply(hs, function(h) mean(vapply(seq_len(K), function(f) {
        trP <- as.character(P$id) %in% subj[fold != f]
        trQ <- as.character(Q$id) %in% subj[fold != f]
        w <- W(P$u[trP] / h) / h
        v <- W(Q$v[trQ] / h) / h * rq[trQ]
        A <- colSums(w * cbind(1, P$x[trP]) * P$y[trP]) / sum(w)
        M <- crossprod(cbind(1, Q$a[trQ, 1]) * v, cbind(1, Q$b[trQ, 1])) / sum(v)
        b <- solve(M, A)
        H <- P[!trP, ]
        ell <- (H$y - b[1] - b[2] * H$x)^2
        coef(lm(ell ~ vg$gamma(abs(H$u))))[[1]]
    }, 1)), 1)
    expect_equal(cv$cv_results$cvloss, hand, tolerance = 1e-6)
    expect_equal(cv$h_cv, hs[which.min(cv$cv_results$cvloss)])
    ref <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = cv$h_cv, correct = "pairs")
    expect_equal(coef(cv$fit), coef(ref))
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "pairs-corrected")
})

test_that("the pair-moment default grid floor covers both kinds of pair", {
    skip_on_cran()
    set.seed(25)
    d <- sim_async_data(n = 100)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        n_h = 8, K = 4, seed = 2, quiet = TRUE,
                                        method = "sync"))
    P <- hand_pairs(d)
    Q <- hand_cov_pairs(d)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    fl <- function(lag, id) {
        mi <- tapply(abs(lag), factor(as.character(id), levels = subj), min)
        mi[is.na(mi)] <- Inf
        max(vapply(1:4, function(f) {
            tr <- cv$fold_id != f
            sort(mi[tr])[ceiling(sqrt(sum(tr)))]
        }, 1))
    }
    expect_equal(cv$h_grid[1], max(fl(P$u, P$id), fl(Q$v, Q$id)), tolerance = 1e-6)
    expect_length(cv$h_grid, 8L)
})

test_that("the default selects the h minimising tr(Sigma_X V) of the pair-moment fit", {
    skip_on_cran()
    set.seed(26)
    d <- sim_async_data(n = 90)
    hs <- c(0.1, 0.25, 0.6)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        h_grid = hs, quiet = TRUE))
    expect_identical(cv$method, "variance")
    expect_identical(cv$correct, "pairs")
    expect_null(cv$fold_id)

    # Second moment of X at the response times: covariate rows weighted by the
    # intensity ratio at their own time.
    rx <- hand_ratio(d, d$x$time)
    X <- cbind(1, d$x$x)
    Sx <- crossprod(X * sqrt(rx)) / sum(rx)
    hand <- vapply(hs, function(h) {
        V <- vcov(kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h, correct = "pairs"))
        sum(Sx * V)
    }, 1)
    expect_equal(cv$cv_results$trSV, hand, tolerance = 1e-6)
    expect_equal(cv$h_cv, hs[which.min(hand)])
    ref <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = cv$h_cv, correct = "pairs")
    expect_equal(coef(cv$fit), coef(ref))
    expect_equal(vcov(cv$fit), vcov(ref))
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "Minimum-variance")
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    expect_invisible(plot(cv))
})

test_that("the minimum-variance default grid starts where sqrt(n) subjects have both kinds of pair", {
    skip_on_cran()
    set.seed(27)
    d <- sim_async_data(n = 100)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        n_h = 8, quiet = TRUE))
    P <- hand_pairs(d)
    Q <- hand_cov_pairs(d)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    fl <- function(lag, id) {
        mi <- tapply(abs(lag), factor(as.character(id), levels = subj), min)
        mi[is.na(mi)] <- Inf
        sort(mi)[ceiling(sqrt(length(subj)))]
    }
    expect_equal(cv$h_grid[1], max(fl(P$u, P$id), fl(Q$v, Q$id)), tolerance = 1e-6)
    expect_length(cv$h_grid, 8L)
})

test_that("the minimum-variance rule takes the Gaussian correction", {
    skip_on_cran()
    set.seed(28)
    d <- sim_async_data(n = 90)
    hs <- c(0.1, 0.25, 0.6)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        h_grid = hs, quiet = TRUE, correct = "gaussian"))
    X <- cbind(1, d$x$x)
    Sx <- crossprod(X) / nrow(X)
    hand <- vapply(hs, function(h) {
        sum(Sx * vcov(kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h,
                                correct = "gaussian")))
    }, 1)
    expect_equal(cv$cv_results$trSV, hand, tolerance = 1e-6)
    expect_identical(cv$fit$correct, "gaussian")
})

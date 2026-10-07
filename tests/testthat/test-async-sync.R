## The asynchrony-corrected fit (kee_async(correct = TRUE)) and the
## synchronous-loss cross-validation that selects its bandwidth
## (kee_async_cv(method = "sync")).

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
    fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h, correct = TRUE)

    P <- hand_pairs(d)
    vg <- hand_vario(d)
    w <- pmax(0.75 * (1 - (P$u / h)^2), 0) / h
    rho <- pmin(pmax(1 - vg$gamma(abs(P$u)) / (2 * vg$s2), -1), 1)
    X <- cbind(1, P$x)
    Z <- cbind(1, vg$m + rho * (P$x - vg$m))
    b <- solve(crossprod(X * w, Z), crossprod(X * w, P$y))
    expect_equal(unname(coef(fit)), as.vector(b), tolerance = 1e-8)
    expect_true(fit$correct)
    expect_match(fit$model, "asynchrony-corrected")
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
    corr <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.6, correct = TRUE)
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
                        correct = TRUE)
        expect_equal(c1$convergence, 0L)
        expect_lt(max(abs(coef(c1) - coef(p))), 0.1)
    }
})

test_that("the correction refuses more than one covariate", {
    set.seed(14)
    d <- sim_async_data(n = 40, beta = c(0.5, 1, -1))
    expect_error(
        kee_async(d$y, d$x, y ~ x1 + x2, id = id, time = time, h = 0.3, correct = TRUE),
        "single covariate"
    )
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x1 + x2, id = id, time = time,
                     h_grid = c(0.2, 0.3), K = 2, quiet = TRUE),
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
                                        h_grid = hs, K = K, seed = 3, quiet = TRUE))
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
    ref <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = cv$h_cv, correct = TRUE)
    expect_equal(coef(cv$fit), coef(ref))
    expect_true(cv$fit$correct)
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "synchronous loss")
})

test_that("the default grid starts where every training fold has sqrt(n_train) subjects in the kernel", {
    skip_on_cran()
    set.seed(16)
    d <- sim_async_data(n = 100)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                        n_h = 8, K = 4, seed = 2, quiet = TRUE))
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

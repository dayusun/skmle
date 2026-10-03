## The nearest-prior ("locf") and estimating-equation ("eeq") criteria of
## kee_async_cv(), rebuilt by hand from the data frames.

cv_fits <- function(d, hs, folds) {
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    lapply(hs, function(h) {
        lapply(seq_len(ncol(folds)), function(r) {
            lapply(seq_len(max(folds[, r])), function(f) {
                tr <- subj[folds[, r] != f]
                list(te = subj[folds[, r] == f], b = coef(kee_async(
                    d$y[as.character(d$y$id) %in% tr, ],
                    d$x[as.character(d$x$id) %in% tr, ],
                    y ~ x, id = id, time = time, h = h
                )))
            })
        })
    })
}

folds_for <- function(seed, n_subj, K, R) {
    set.seed(seed)
    matrix(replicate(R, sample(rep_len(seq_len(K), n_subj))), nrow = n_subj)
}

test_that("locf pairs each response with the last covariate at or before it", {
    skip_on_cran()
    set.seed(81)
    d <- sim_async_data(n = 40)
    hs <- c(0.1, 0.3)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time, h_grid = hs, K = 2, seed = 4, quiet = TRUE,
        criterion = "locf", repeats = 2
    ))
    folds <- folds_for(4, 40, 2, 2)
    fits <- cv_fits(d, hs, folds)
    pair_x <- function(id, t) {
        xs <- d$x[d$x$id == id, ]
        if (!nrow(xs)) return(NA_real_)
        before <- xs[xs$time <= t, ]
        if (nrow(before)) before$x[which.max(before$time)] else xs$x[which.min(xs$time)]
    }
    hand <- vapply(seq_along(hs), function(m) {
        mean(vapply(seq_len(2), function(r) {
            num <- 0
            den <- 0
            for (fit in fits[[m]][[r]]) {
                yt <- d$y[as.character(d$y$id) %in% fit$te, ]
                for (j in seq_len(nrow(yt))) {
                    xv <- pair_x(yt$id[j], yt$time[j])
                    if (is.na(xv)) next
                    num <- num + (yt$y[j] - fit$b[1] - fit$b[2] * xv)^2
                    den <- den + 1
                }
            }
            num / den
        }, numeric(1)))
    }, numeric(1))
    expect_equal(cv$cv_results$cvloss, hand, tolerance = 1e-8)
    expect_null(cv$h_eval)
})

test_that("eeq scores the squared norm of the held-out estimating function", {
    skip_on_cran()
    set.seed(83)
    d <- sim_async_data(n = 40)
    hs <- c(0.1, 0.3)
    he <- 0.15
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time, h_grid = hs, K = 2, seed = 6, quiet = TRUE,
        criterion = "eeq", h_eval = he, repeats = 2
    ))
    expect_equal(cv$h_eval, he)
    folds <- folds_for(6, 40, 2, 2)
    fits <- cv_fits(d, hs, folds)
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    hand <- vapply(seq_along(hs), function(m) {
        mean(vapply(seq_len(2), function(r) {
            mean(vapply(fits[[m]][[r]], function(fit) {
                yt <- d$y[as.character(d$y$id) %in% fit$te, ]
                xt <- d$x[as.character(d$x$id) %in% fit$te, ]
                U <- c(0, 0)
                D <- 0
                for (j in seq_len(nrow(yt))) {
                    k <- which(xt$id == yt$id[j])
                    if (!length(k)) next
                    w <- W((yt$time[j] - xt$time[k]) / he) / he
                    res <- yt$y[j] - fit$b[1] - fit$b[2] * xt$x[k]
                    U <- U + c(sum(w * res), sum(w * res * xt$x[k]))
                    D <- D + sum(w)
                }
                sum((U / D)^2)
            }, numeric(1)))
        }, numeric(1)))
    }, numeric(1))
    expect_equal(cv$cv_results$cvloss, hand, tolerance = 1e-8)

    # The default evaluation bandwidth for "eeq" is the lower end of the
    # default grid.
    cv2 <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time, h_grid = hs, K = 2, seed = 6, quiet = TRUE,
        criterion = "eeq"
    ))
    pooled <- c(d$y$time, d$x$time)
    expect_equal(cv2$h_eval,
        2 * diff(quantile(pooled, c(.25, .75), names = FALSE)) * 40^(-0.7)
    )
})

test_that("locf pairs fall back to the first covariate when none is earlier", {
    p <- list(
        n = 1L, yi = c(1L, 1L), yt = c(0.1, 0.6), y = c(1, 2),
        xi = c(1L, 1L), xt = c(0.3, 0.5), X = cbind(1, c(10, 20))
    )
    ac <- skmle:::.async_locf_pairs(p)
    # The response at 0.1 precedes both covariates: first covariate (row 1).
    # The response at 0.6 follows both: last at or before (row 2).
    expect_equal(ac$omega, c(1, 1))
    expect_equal(ac$c, c(1, 2))
    expect_equal(ac$sy2, 5)
})

## Bandwidth selection for the asynchronous estimators.
##
## The fixed criterion scores every candidate bandwidth on the same held-out
## pairs, formed at an evaluation bandwidth, so the loss is a sum of
## per-subject terms. These tests rebuild it by hand and check the selection
## rules built on it.

# Held-out loss by hand: for every split, refit without each fold and score its
# subjects. `h_score` NULL means score at the candidate's own h (matched).
hand_cv <- function(d, hs, folds, h_score = NULL, pooled = TRUE) {
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    one_h <- function(h) {
        per_rep <- apply(folds, 2, function(fold) {
            nd <- vapply(seq_len(max(fold)), function(f) {
                te <- subj[fold == f]
                tr <- subj[fold != f]
                b <- coef(kee_async(
                    d$y[as.character(d$y$id) %in% tr, ],
                    d$x[as.character(d$x$id) %in% tr, ],
                    y ~ x,
                    id = id, time = time, h = h
                ))
                hs_ <- if (is.null(h_score)) h else h_score
                yt <- d$y[as.character(d$y$id) %in% te, ]
                xt <- d$x[as.character(d$x$id) %in% te, ]
                num <- 0
                den <- 0
                for (j in seq_len(nrow(yt))) {
                    k <- which(xt$id == yt$id[j])
                    if (!length(k)) next
                    w <- W((yt$time[j] - xt$time[k]) / hs_) / hs_
                    num <- num + sum(w * (yt$y[j] - b[1] - b[2] * xt$x[k])^2)
                    den <- den + sum(w)
                }
                c(num, den)
            }, numeric(2))
            if (pooled) sum(nd[1, ]) / sum(nd[2, ]) else mean(nd[1, ] / nd[2, ])
        })
        mean(per_rep)
    }
    vapply(hs, one_h, numeric(1))
}

folds_for <- function(seed, n_subj, K, R) {
    set.seed(seed)
    matrix(replicate(R, sample(rep_len(seq_len(K), n_subj))), nrow = n_subj)
}

test_that("the fixed criterion pools every held-out pair at h_eval", {
    skip_on_cran()
    set.seed(53)
    d <- sim_async_data(n = 60)
    hs <- c(0.15, 0.4)
    # Two candidates: an endpoint minimum is certain and the warning expected.
    cv <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = hs, K = 2, seed = 11, quiet = TRUE, h_eval = 0.1,
        criterion = "fixed"
    ))
    expect_identical(cv$criterion, "fixed")
    expect_equal(cv$h_eval, 0.1)
    folds <- folds_for(11, 60, 2, 1)
    expect_equal(cv$cv_results$cvloss, hand_cv(d, hs, folds, h_score = 0.1),
        tolerance = 1e-8
    )

    # "matched" keeps the earlier behaviour: each candidate at its own h,
    # averaged over fold ratios.
    cvm <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = hs, K = 2, seed = 11, quiet = TRUE, criterion = "matched"
    ))
    expect_equal(cvm$cv_results$cvloss, hand_cv(d, hs, folds, pooled = FALSE),
        tolerance = 1e-8
    )
    expect_null(cvm$h_eval)

    expect_null(cv$fit$call$h_eval)
    expect_null(cv$fit$call$criterion)
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "h_eval")
})

test_that("repeats average the loss over independent splits", {
    skip_on_cran()
    set.seed(67)
    d <- sim_async_data(n = 50)
    hs <- c(0.12, 0.3)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = hs, K = 2, seed = 5, quiet = TRUE, h_eval = 0.1,
        criterion = "fixed", repeats = 3
    ))
    expect_equal(dim(cv$fold_id), c(50L, 3L))
    folds <- folds_for(5, 50, 2, 3)
    expect_equal(cv$cv_results$cvloss, hand_cv(d, hs, folds, h_score = 0.1),
        tolerance = 1e-8
    )
    expect_null(cv$fit$call$repeats)
})

test_that("several evaluation bandwidths average their losses", {
    skip_on_cran()
    set.seed(71)
    d <- sim_async_data(n = 50)
    hs <- c(0.12, 0.3)
    args <- list(d$y, d$x, y ~ x,
        id = quote(id), time = quote(time),
        h_grid = hs, K = 2, seed = 3, quiet = TRUE, criterion = "fixed"
    )
    run <- function(he) {
        suppressWarnings(do.call(kee_async_cv, c(args, h_eval = list(he))))$cv_results$cvloss
    }
    expect_equal(run(c(0.08, 0.15)), (run(0.08) + run(0.15)) / 2, tolerance = 1e-10)

    lower <- suppressWarnings(do.call(kee_async_cv, c(args, h_eval = "lower")))
    pooled <- c(d$y$time, d$x$time)
    expect_equal(lower$h_eval,
        2 * diff(quantile(pooled, c(.25, .75), names = FALSE)) * 50^(-0.7)
    )
})

test_that("the fixed criterion uses the canonical deviance for nonlinear links", {
    set.seed(59)
    d <- sim_async_data(n = 80, beta = c(0.3, 1.0), link = "logistic")
    b <- c(0.2, 0.8)
    h <- 0.2
    p <- skmle:::.async_parse(y ~ x, d$y, d$x, "id", "time", TRUE)
    ac <- skmle:::.async_collapse(p, h, FALSE)

    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    num <- 0
    den <- 0
    for (j in seq_along(p$yt)) {
        k <- which(p$xi == p$yi[j])
        if (!length(k)) next
        w <- W((p$yt[j] - p$xt[k]) / h) / h
        eta <- as.numeric(p$X[k, , drop = FALSE] %*% b)
        num <- num + sum(w * (log1p(exp(eta)) - p$y[j] * eta))
        den <- den + sum(w)
    }
    expect_equal(skmle:::.async_deviance(ac, p$X, b, 2L), num / den,
        tolerance = 1e-10
    )
})

test_that("the selection rules pick what they say", {
    sel <- skmle:::.async_select
    h <- exp(seq(log(0.01), log(1), length.out = 9))
    # A flat bottom from index 3 to 7 with a noise dip at 3.
    L <- c(5, 3, 1.00, 1.01, 1.02, 1.01, 1.03, 3, 5)
    expect_equal(sel(L, NULL, "min", h)$index, 3L)

    # Paired SE of 0.05 everywhere: indices 3..7 tie, centre is 5.
    U <- matrix(0, 4, 9)
    U[, ] <- rep(L / 4, each = 4)
    U <- U + outer(c(-1, 1, -1, 1) * 0.025, rep(1, 9)) * outer(rep(1, 4), (1:9 != 3))
    s_flat <- sel(L, U, "flat", h)
    expect_equal(s_flat$index, 5L)
    expect_true(all(s_flat$diff_se[-3] > 0))

    # Valley fit: falls then rises; lowest plateau centred.
    s_uni <- sel(c(5, 3, 1, 1.2, 1, 1.1, 3, 4, 5), NULL, "unimodal", h)
    expect_true(s_uni$index %in% 3:6)
    f <- s_uni$fitted
    expect_true(all(diff(f[1:which.min(f)]) <= 1e-12))
    expect_true(all(diff(f[which.min(f):9]) >= -1e-12))

    s_sm <- sel(L, NULL, "smooth", h)
    expect_true(s_sm$index %in% 3:7)
})

test_that("flat selection needs the fixed criterion, and h_eval is checked", {
    set.seed(61)
    d <- sim_async_data(n = 30)
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x,
            id = id, time = time,
            h_grid = c(0.2, 0.3), K = 2, seed = 1, quiet = TRUE, select = "flat"
        ),
        "needs criterion"
    )
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x,
            id = id, time = time,
            h_grid = c(0.2, 0.3), K = 2, seed = 1, quiet = TRUE, h_eval = 1e-9,
            criterion = "fixed"
        ),
        "h_eval"
    )
})

test_that("flat selection runs end to end and reports paired SEs", {
    skip_on_cran()
    set.seed(73)
    d <- sim_async_data(n = 80)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time,
        h_grid = exp(seq(log(0.03), log(0.6), length.out = 8)), K = 4,
        seed = 2, quiet = TRUE, criterion = "fixed", repeats = 2, select = "flat"
    ))
    r <- cv$cv_results
    expect_true(all(is.finite(r$diff_se)))
    i <- which(r$h == cv$h_cv)
    expect_lte(r$cvloss[i] - min(r$cvloss), r$diff_se[i] + 1e-12)
})

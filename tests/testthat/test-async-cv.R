## Bandwidth selection for the asynchronous estimators: plain K-fold
## cross-validation over subjects on a grid built from the data
## (method = "plain"). The synchronous-loss method is in test-async-sync.R.

# Held-out loss by hand: for every split, refit without each fold and score its
# subjects at the candidate's own bandwidth; mean over folds, then splits.
hand_cv <- function(d, hs, folds) {
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    vapply(hs, function(h) {
        mean(apply(folds, 2, function(fold) {
            mean(vapply(seq_len(max(fold)), function(f) {
                te <- subj[fold == f]
                tr <- subj[fold != f]
                b <- coef(kee_async(
                    d$y[as.character(d$y$id) %in% tr, ],
                    d$x[as.character(d$x$id) %in% tr, ],
                    y ~ x,
                    id = id, time = time, h = h
                ))
                yt <- d$y[as.character(d$y$id) %in% te, ]
                xt <- d$x[as.character(d$x$id) %in% te, ]
                num <- 0
                den <- 0
                for (j in seq_len(nrow(yt))) {
                    k <- which(xt$id == yt$id[j])
                    if (!length(k)) next
                    w <- W((yt$time[j] - xt$time[k]) / h) / h
                    num <- num + sum(w * (yt$y[j] - b[1] - b[2] * xt$x[k])^2)
                    den <- den + sum(w)
                }
                num / den
            }, numeric(1)))
        }))
    }, numeric(1))
}

# Every within-subject response-minus-covariate lag.
all_lags <- function(d) {
    unlist(lapply(unique(d$y$id), function(i) {
        outer(d$y$time[d$y$id == i], d$x$time[d$x$id == i], "-")
    }), use.names = FALSE)
}

test_that("the default grid runs from the smallest to the largest within-subject lag", {
    skip_on_cran()
    set.seed(91)
    d <- sim_async_data(n = 40)
    lg <- all_lags(d)

    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time, n_h = 12, K = 2, seed = 1, quiet = TRUE,
        method = "plain"
    ))
    two <- abs(lg[lg != 0])
    expect_equal(cv$h_grid, exp(seq(log(min(two)), log(max(two)), length.out = 12)))

    cv1 <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time, n_h = 12, K = 2, seed = 1, quiet = TRUE,
        one_sided = TRUE, method = "plain"
    ))
    pos <- lg[lg > 0]
    expect_equal(cv1$h_grid, exp(seq(log(min(pos)), log(max(pos)), length.out = 12)))
})

test_that("repeats average the plain held-out loss over independent splits", {
    skip_on_cran()
    set.seed(67)
    d <- sim_async_data(n = 50)
    hs <- c(0.12, 0.3)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = hs, K = 2, seed = 5, quiet = TRUE, repeats = 3,
        method = "plain"
    ))
    expect_equal(dim(cv$fold_id), c(50L, 3L))
    set.seed(5)
    folds <- matrix(replicate(3, sample(rep_len(1:2, 50))), nrow = 50)
    expect_equal(cv$cv_results$cvloss, hand_cv(d, hs, folds), tolerance = 1e-8)
    expect_equal(cv$h_cv, hs[which.min(cv$cv_results$cvloss)])
    expect_null(cv$fit$call$repeats)
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "3 random splits")
})

test_that("the experimental criteria are gone", {
    set.seed(61)
    d <- sim_async_data(n = 30)
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x,
            id = id, time = time, h_grid = c(0.2, 0.3), K = 2, quiet = TRUE,
            criterion = "fixed"
        ),
        "unused argument"
    )
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x,
            id = id, time = time, h_grid = c(0.2, 0.3), K = 2, quiet = TRUE,
            h_eval = 0.1
        ),
        "unused argument"
    )
})

test_that("a singular fit at a tiny bandwidth fails cleanly instead of crashing", {
    # The default grid starts at the smallest lag, where one or two pairs make
    # the identity-link system singular. The solver used to leave beta empty and
    # crash in its final pass; it now reports non-convergence.
    set.seed(91)
    d <- sim_async_data(n = 40)
    p <- skmle:::.async_parse(y ~ x, d$y, d$x, "id", "time", TRUE)
    g <- skmle:::.async_lag_grid(p, 12, FALSE)
    fits <- lapply(g[1:4], function(h) skmle:::.async_fit_ti(p, h, FALSE, 0L, 50L, 1e-8))
    conv <- vapply(fits, function(f) if (is.null(f)) NA_integer_ else f$convergence, 1L)
    expect_true(any(conv == 1L, na.rm = TRUE))
    for (f in fits) if (!is.null(f)) expect_length(f$beta, 2L)
})

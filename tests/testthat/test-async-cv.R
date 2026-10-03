## Bandwidth selection for the asynchronous estimators.
##
## The default criterion scores every candidate bandwidth on the same held-out
## pairs, formed at one evaluation bandwidth, so the comparison across
## candidates changes only beta-hat. These tests rebuild the held-out loss by
## hand under both criteria.

test_that("the fixed criterion scores every candidate on the same pairs", {
    skip_on_cran()
    set.seed(53)
    d <- sim_async_data(n = 60)
    hs <- c(0.15, 0.4)
    he <- 0.1
    K <- 2
    # Two candidates: an endpoint minimum is certain and the warning expected.
    cv <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = hs, K = K, seed = 11, quiet = TRUE, h_eval = he
    ))
    expect_identical(cv$criterion, "fixed")
    expect_equal(cv$h_eval, he)

    set.seed(11)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    fold <- sample(rep_len(seq_len(K), length(subj)))
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    loss_f <- function(f, h, h_score) {
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
            w <- W((yt$time[j] - xt$time[k]) / h_score) / h_score
            num <- num + sum(w * (yt$y[j] - b[1] - b[2] * xt$x[k])^2)
            den <- den + sum(w)
        }
        num / den
    }
    score <- function(h, h_score) {
        mean(vapply(seq_len(K), loss_f, numeric(1), h = h, h_score = h_score))
    }
    # Every candidate is scored at h_eval, not at its own bandwidth.
    hand <- vapply(hs, score, numeric(1), h_score = he)
    expect_equal(cv$cv_results$cvloss, hand, tolerance = 1e-8)

    # "matched" keeps the earlier behaviour: each candidate at its own h.
    cvm <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = hs, K = K, seed = 11, quiet = TRUE, criterion = "matched"
    ))
    hand_m <- vapply(hs, function(h) score(h, h), numeric(1))
    expect_equal(cvm$cv_results$cvloss, hand_m, tolerance = 1e-8)
    expect_null(cvm$h_eval)

    # The refit call carries neither CV argument.
    expect_null(cv$fit$call$h_eval)
    expect_null(cv$fit$call$criterion)
    expect_match(paste(capture.output(print(cv)), collapse = "\n"), "h_eval")
})

test_that("the fixed criterion uses the canonical deviance for nonlinear links", {
    set.seed(59)
    d <- sim_async_data(n = 80, beta = c(0.3, 1.0), link = "logistic")
    b <- c(0.2, 0.8)
    h <- 0.2
    p <- skmle:::.async_parse(y ~ x, d$y, d$x, "id", "time", TRUE)
    ac <- skmle:::.async_collapse(p, h, FALSE)

    # By hand over every within-subject pair inside the window.
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

test_that("an evaluation bandwidth that admits no pair is refused", {
    set.seed(61)
    d <- sim_async_data(n = 30)
    expect_error(
        kee_async_cv(d$y, d$x, y ~ x,
            id = id, time = time,
            h_grid = c(0.2, 0.3), K = 2, seed = 1, quiet = TRUE, h_eval = 1e-9
        ),
        "h_eval"
    )
})

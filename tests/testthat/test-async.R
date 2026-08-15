## Asynchronous longitudinal regression, Cao, Zeng and Fine (2015).
##
## The data generator sim_async_data() follows their Section 4: response and
## covariate observation times are independent Poisson streams on (0, 1), X is
## Gaussian with correlation exp(-|s - t|) -- an Ornstein-Uhlenbeck covariance
## whose diagonal has a KINK, which is exactly the one-sided-differentiable
## class their Section 6 leaves open -- and the error has covariance 2^(-|s-t|).
##
## Every fit below is checked against a plain-R transcription of the estimating
## equation as it is written in the paper, so the C++ collapse-onto-covariate-
## rows shortcut is never trusted on its own.

## ---------------------------------------------------------------- invariant

test_that("kee_async matches an independent transcription of the estimating equation", {
    set.seed(3)
    d <- sim_async_data(n = 80)
    h <- 0.25
    W <- function(z) pmax(0.75 * (1 - z^2), 0)

    fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h)

    # Accumulate A and b over (response, covariate) pairs within subject,
    # exactly as the estimating equation is written.
    A <- matrix(0, 2, 2)
    b <- numeric(2)
    for (j in seq_len(nrow(d$y))) {
        k <- which(d$x$id == d$y$id[j])
        if (!length(k)) next
        z <- (d$y$time[j] - d$x$time[k]) / h
        s <- abs(z) <= 1
        if (!any(s)) next
        w <- W(z[s]) / h
        xv <- cbind(1, d$x$x[k][s])
        A <- A + crossprod(xv, xv * w)
        b <- b + as.numeric((w * d$y$y[j]) %*% xv)
    }
    expect_equal(unname(coef(fit)), as.numeric(solve(A, b)), tolerance = 1e-10)
})

test_that("the half kernel drops every pair with the covariate after the response", {
    set.seed(21)
    d <- sim_async_data(n = 60)
    h <- 0.25

    full <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h)
    half <- kee_async(d$y, d$x, y ~ x, id = id, time = time,
        h = h, one_sided = TRUE
    )

    # Same kernel shape, so the half kernel keeps a strict subset of the pairs
    # and the two fits must differ.
    expect_lt(half$npair, full$npair)
    expect_false(isTRUE(all.equal(coef(full), coef(half))))

    # And it keeps exactly the pairs a strictly-past window admits.
    keep <- 0L
    for (j in seq_len(nrow(d$y))) {
        k <- which(d$x$id == d$y$id[j])
        z <- (d$y$time[j] - d$x$time[k]) / h
        keep <- keep + sum(z > 0 & z <= 1)
    }
    expect_equal(half$npair, keep)
})

test_that("kee_async recovers the truth, with a bias that grows in h", {
    skip_on_cran()
    set.seed(11)
    d <- sim_async_data(n = 400)
    hs <- c(0.05, 0.20, 0.45)

    bhat <- function(h) {
        coef(kee_async(d$y, d$x, y ~ x, id = id, time = time, h = h))[2]
    }
    b <- vapply(hs, bhat, numeric(1))

    expect_equal(unname(b[1]), 1.5, tolerance = 0.1)
    # The default covariate process is Ornstein-Uhlenbeck, whose covariance has
    # a kink on the diagonal.  That is the case Cao, Zeng and Fine's Section 6
    # leaves open, and it leaves an O(h) bias a symmetric kernel cannot remove,
    # so the error must grow with the bandwidth.
    expect_gt(abs(b[3] - 1.5), abs(b[1] - 1.5))
})

test_that("kee_async handles nonlinear links and rejects malformed input", {
    skip_on_cran()
    set.seed(5)
    d <- sim_async_data(n = 250, beta = c(0.3, 1.0), link = "logistic")

    fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time,
        h = 0.4, link = "logistic"
    )
    expect_true(all(is.finite(coef(fit))))
    expect_length(coef(fit), 2)
    expect_identical(fit$convergence, 0L)

    # Newton must land on a root of the estimating equation.
    U <- numeric(2)
    h <- 0.4
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    bb <- coef(fit)
    for (j in seq_len(nrow(d$y))) {
        k <- which(d$x$id == d$y$id[j])
        z <- (d$y$time[j] - d$x$time[k]) / h
        s <- abs(z) <= 1
        if (!any(s)) next
        w <- W(z[s]) / h
        xv <- cbind(1, d$x$x[k][s])
        U <- U + as.numeric((w * (d$y$y[j] - plogis(as.numeric(xv %*% bb)))) %*% xv)
    }
    expect_lt(max(abs(U)), 1e-6)

    expect_error(
        kee_async(d$y, d$x, y ~ x, id = id, time = time, h = -1),
        "positive"
    )
    expect_error(
        kee_async(d$y, d$x, y ~ x, id = nope, time = time, h = 0.4),
        "not found in data_y"
    )
    expect_error(
        kee_async(d$y, d$x, y ~ nope, id = id, time = time, h = 0.4),
        "not found"
    )
    expect_error(
        kee_async(as.list(d$y), d$x, y ~ x, id = id, time = time, h = 0.4),
        "data.frame"
    )
})

## ---------------------------------------------------------------- dependent

test_that("kee_async_td matches an independent transcription at each target time", {
    set.seed(13)
    d <- sim_async_data(n = 120, beta = function(tt) cbind(0.5, 1 + tt))
    times <- c(0.3, 0.5, 0.7)
    h <- 0.3
    W <- function(z) pmax(0.75 * (1 - z^2), 0)

    fit <- kee_async_td(d$y, d$x, y ~ x, id = id, time = time,
        times = times, h = h
    )

    for (m in seq_along(times)) {
        t0 <- times[m]
        A <- matrix(0, 2, 2)
        b <- numeric(2)
        for (j in seq_len(nrow(d$y))) {
            k <- which(d$x$id == d$y$id[j])
            if (!length(k)) next
            # The two occasions are weighted by their SEPARATE distances to
            # the target time; the lag between them never appears.  Both are
            # measured as "how far in the past", the package-wide convention.
            ww <- (W((t0 - d$y$time[j]) / h) / h) * (W((t0 - d$x$time[k]) / h) / h)
            if (all(ww == 0)) next
            xv <- cbind(1, d$x$x[k])
            A <- A + crossprod(xv, xv * ww)
            b <- b + as.numeric((ww * d$y$y[j]) %*% xv)
        }
        expect_equal(unname(fit$coefficients[m, ]), as.numeric(solve(A, b)),
            tolerance = 1e-10
        )
    }
})

test_that("kee_async_td tracks a coefficient curve", {
    skip_on_cran()
    set.seed(19)
    # A smooth covariance keeps the bias at O(h^2), so the curve is recovered
    # rather than attenuated by the kink term.
    d <- sim_async_data(
        n = 900, beta = function(tt) cbind(0.5, 1 + tt),
        lambda_y = 8, lambda_x = 8,
        x_cov = function(s, t) exp(-4 * (s - t)^2)
    )
    times <- seq(0.25, 0.75, by = 0.25)
    fit <- kee_async_td(d$y, d$x, y ~ x, id = id, time = time,
        times = times, h = 0.2
    )

    expect_true(all(fit$convergence == 0L))
    # Monotone increasing truth, so the estimate must increase too.
    expect_true(all(diff(fit$coefficients[, 2]) > 0))
    expect_equal(fit$coefficients[, 2], 1 + times, tolerance = 0.15, ignore_attr = TRUE)
    expect_equal(fit$coefficients[, 1], rep(0.5, 3), tolerance = 0.15, ignore_attr = TRUE)
    expect_true(all(fit$se > 0))
})

test_that("kee_async_td validates input and reports empty windows", {
    set.seed(23)
    d <- sim_async_data(n = 40)
    a <- list(d$y, d$x, y ~ x, id = "id", time = "time")

    expect_error(do.call(kee_async_td, c(a, list(times = 0.5, h = 0))), "positive")
    expect_error(do.call(kee_async_td, c(a, list(times = 0.5, h = 0.2, h2 = -1))), "positive")
    expect_error(do.call(kee_async_td, c(a, list(times = numeric(0), h = 0.2))), "non-empty")

    # A target time far outside the data has no weight on either side.
    expect_error(
        do.call(kee_async_td, c(a, list(times = 50, h = 0.05))),
        "no target time has data"
    )

    expect_warning(
        fit <- do.call(kee_async_td, c(a, list(times = c(0.5, 50), h = 0.05))),
        "no data in their window"
    )
    expect_identical(fit$convergence, c(0L, 3L))
    expect_true(all(is.na(fit$coefficients[2, ])))
})

test_that("print and plot methods work on a kee_td fit", {
    set.seed(29)
    d <- sim_async_data(n = 60)
    fit <- kee_async_td(d$y, d$x, y ~ x, id = id, time = time,
        times = c(0.4, 0.6), h = 0.3
    )
    expect_output(print(fit), "Time-dependent coefficients")
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    expect_invisible(plot(fit))
    expect_error(plot(fit, which = "nope"), "selects no coefficient")
})

test_that("one_sided restricts the covariate side of the time-dependent fit", {
    set.seed(31)
    d <- sim_async_data(n = 150)
    a <- list(d$y, d$x, y ~ x, id = "id", time = "time")

    full <- do.call(kee_async_td, c(a, list(times = c(0.4, 0.6), h = 0.3)))
    half <- do.call(kee_async_td, c(a, list(
        times = c(0.4, 0.6), h = 0.3, one_sided = TRUE
    )))
    expect_false(isTRUE(all.equal(coef(full), coef(half))))
    expect_true(all(is.finite(coef(half))))

    # Transcribe the half-kernel equation directly: the covariate side keeps
    # only observations strictly before the target time.
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    h <- 0.3
    t0 <- 0.4
    A <- matrix(0, 2, 2)
    b <- numeric(2)
    for (j in seq_len(nrow(d$y))) {
        k <- which(d$x$id == d$y$id[j])
        if (!length(k)) next
        lag_x <- t0 - d$x$time[k]
        vx <- W(lag_x / h) / h
        vx[lag_x <= 0] <- 0
        ww <- (W((t0 - d$y$time[j]) / h) / h) * vx
        if (all(ww == 0)) next
        xv <- cbind(1, d$x$x[k])
        A <- A + crossprod(xv, xv * ww)
        b <- b + as.numeric((ww * d$y$y[j]) %*% xv)
    }
    expect_equal(unname(coef(half)[1, ]), as.numeric(solve(A, b)), tolerance = 1e-10)
})

test_that("one_sided is validated", {
    set.seed(3)
    d <- sim_async_data(n = 30)
    a <- list(d$y, d$x, y ~ x, id = "id", time = "time")
    expect_error(do.call(kee_async, c(a, list(h = 0.3, one_sided = "yes"))), "TRUE or FALSE")
    expect_error(
        do.call(kee_async_td, c(a, list(times = 0.5, h = 0.3, one_sided = NA))),
        "TRUE or FALSE"
    )
})

## ------------------------------------------------------- bandwidth selection

test_that("kee_async_cv selects a bandwidth and refits there", {
    skip_on_cran()
    set.seed(41)
    d <- sim_async_data(n = 200)
    # The minimum may land on a grid endpoint here; that is exactly what the
    # endpoint warning is for, so accept it rather than letting it leak.
    cv <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = c(0.10, 0.20, 0.35), K = 4, seed = 1, quiet = TRUE
    ))
    expect_s3_class(cv, "cv.kee_async")
    expect_true(cv$h_cv %in% cv$h_grid)
    expect_equal(cv$fit$h, cv$h_cv)
    expect_s3_class(cv$fit, "kee")
    expect_true(all(cv$cv_results$nfold_used > 0))

    # The loss is a weighted MEAN of squared residuals, so it sits on the scale
    # of the error variance (1 in this generator) whatever the bandwidth.  An
    # unnormalised loss would instead grow with the number of pairs a wider
    # bandwidth admits, and would always select the smallest candidate.
    expect_true(all(cv$cv_results$cvloss > 0.1 & cv$cv_results$cvloss < 20))
    expect_lt(max(cv$cv_results$cvloss) / min(cv$cv_results$cvloss), 3)

    # The refit call must be readable, not the whole data frame.
    expect_lt(nchar(paste(deparse(cv$fit$call), collapse = " ")), 300)
    expect_lt(length(capture.output(print(cv))), 20)
})

test_that("kee_async_cv reproduces a hand-computed held-out loss", {
    skip_on_cran()
    set.seed(43)
    d <- sim_async_data(n = 60)
    h <- 0.3
    K <- 2
    cv <- kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        h_grid = h, K = K, seed = 7, quiet = TRUE
    )

    # Rebuild fold 1's contribution from scratch.
    set.seed(7)
    subj <- sort(unique(c(as.character(d$y$id), as.character(d$x$id))))
    fold <- sample(rep_len(seq_len(K), length(subj)))
    W <- function(z) pmax(0.75 * (1 - z^2), 0)
    loss_f <- function(f) {
        te <- subj[fold == f]
        tr <- subj[fold != f]
        b <- coef(kee_async(
            d$y[as.character(d$y$id) %in% tr, ],
            d$x[as.character(d$x$id) %in% tr, ],
            y ~ x,
            id = id, time = time, h = h
        ))
        num <- 0
        den <- 0
        yt <- d$y[as.character(d$y$id) %in% te, ]
        xt <- d$x[as.character(d$x$id) %in% te, ]
        for (j in seq_len(nrow(yt))) {
            k <- which(xt$id == yt$id[j])
            if (!length(k)) next
            w <- W((yt$time[j] - xt$time[k]) / h) / h
            mu <- b[1] + b[2] * xt$x[k]
            num <- num + sum(w * (yt$y[j] - mu)^2)
            den <- den + sum(w)
        }
        c(num, den)
    }
    hand <- vapply(seq_len(K), loss_f, numeric(2))
    expect_equal(cv$cv_results$cvloss, mean(hand[1, ] / hand[2, ]), tolerance = 1e-8)
})

test_that("kee_async_cv builds a default grid and flags endpoint minima", {
    skip_on_cran()
    set.seed(47)
    d <- sim_async_data(n = 120)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x,
        y ~ x, id = id, time = time,
        n_h = 6, K = 3, seed = 2, quiet = TRUE
    ))
    expect_length(cv$h_grid, 6L)
    expect_true(all(diff(cv$h_grid) > 0))
    expect_true(all(cv$h_grid > 0))

    # A grid with a single sensible value plus an absurd one must pick the
    # sensible one and warn that the minimum is on the boundary.
    expect_warning(
        kee_async_cv(d$y, d$x,
            y ~ x, id = id, time = time,
            h_grid = c(0.25, 5), K = 3, seed = 2, quiet = TRUE
        ),
        "endpoint"
    )
})

test_that("a bandwidth on the wrong scale is flagged rather than fitted quietly", {
    set.seed(53)
    d <- sim_async_data(n = 80)
    d$y$time <- d$y$time * 365
    d$x$time <- d$x$time * 365
    expect_error(
        kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.001),
        "same scale"
    )
    expect_warning(
        kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 1),
        "same scale as the observation times"
    )
})

test_that("the formula interface handles several covariates and missing values", {
    skip_on_cran()
    set.seed(59)
    d <- sim_async_data(n = 150, beta = c(0.5, 1.5, -0.8))
    fit <- kee_async(d$y, d$x,
        y ~ x1 + x2,
        id = id, time = time, h = 0.3
    )
    expect_named(coef(fit), c("(Intercept)", "x1", "x2"))
    expect_equal(unname(coef(fit)), c(0.5, 1.5, -0.8), tolerance = 0.25)

    # Rows with NA are dropped on the side they occur, keeping id/time aligned.
    dy <- d$y
    dx <- d$x
    dy$y[1:5] <- NA
    dx$x1[1:5] <- NA
    fit2 <- kee_async(dy, dx,
        y ~ x1 + x2,
        id = id, time = time, h = 0.3
    )
    expect_true(all(is.finite(coef(fit2))))
    expect_lt(fit2$npair, fit$npair)

    # id and time may also be given as strings.
    fit3 <- kee_async(d$y, d$x,
        y ~ x1 + x2,
        id = "id", time = "time", h = 0.3
    )
    expect_equal(coef(fit), coef(fit3))
})

## ------------------------------------------------------------------ methods

test_that("the standard model generics work on every fitted class", {
    skip_on_cran()
    set.seed(61)
    d <- sim_async_data(n = 150, beta = function(tt) cbind(0.5, 1 + tt))
    ti <- kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.3)
    td <- kee_async_td(d$y, d$x,
        y ~ x, id = id, time = time,
        times = c(0.3, 0.5, 0.7), h = 0.3
    )

    expect_equal(dim(vcov(ti)), c(2L, 2L))
    expect_equal(vcov(ti), ti$var)
    expect_identical(nobs(ti), 150L)

    ci <- confint(ti)
    expect_equal(dim(ci), c(2L, 2L))
    expect_true(all(ci[, 1] < coef(ti) & coef(ti) < ci[, 2]))

    expect_length(vcov(td), 3L)
    expect_equal(dim(vcov(td, time = 0.5)), c(2L, 2L))
    expect_equal(vcov(td, time = 0.49), vcov(td, time = 0.5))
    expect_identical(nobs(td), 150L)

    cid <- confint(td)
    expect_equal(nrow(cid), 6L)
    expect_true(all(cid$lower < cid$estimate | is.na(cid$lower)))
    expect_error(confint(td, parm = "nope"), "selects no coefficient")
})

## ----------------------------------------------------- defaults for novices

test_that("the asynchronous estimators run with only the data", {
    skip_on_cran()
    set.seed(71)
    d <- sim_async_data(n = 200, beta = c(0.5, 1.5))

    expect_message(
        fit <- kee_async(d$y, d$x, y ~ x, id = id, time = time),
        "rule of thumb"
    )
    expect_equal(unname(coef(fit)), c(0.5, 1.5), tolerance = 0.2)
    expect_true(fit$h > 0)

    # The default is the geometric midpoint of the grid kee_async_cv() searches,
    # which is the claim the message makes.
    pooled <- c(d$y$time, d$x$time)
    iqr <- diff(stats::quantile(pooled, c(0.25, 0.75), names = FALSE))
    expect_equal(fit$h, 2 * iqr * 200^(-0.5), tolerance = 1e-8)

    # kee_async_td() also defaults its target times, inside the data.
    expect_message(td <- kee_async_td(d$y, d$x, y ~ x, id = id, time = time))
    expect_length(td$times, 25L)
    expect_gte(min(td$times), min(d$y$time))
    expect_lte(max(td$times), max(d$y$time))
    expect_true(all(td$convergence == 0L))

    expect_no_message(
        kee_async(d$y, d$x, y ~ x, id = id, time = time, h = 0.25)
    )
})

test_that("the cross-validation curve plots", {
    skip_on_cran()
    set.seed(73)
    d <- sim_async_data(n = 120)
    cv <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x,
        id = id, time = time,
        h_grid = c(0.1, 0.2, 0.3, 0.5), K = 3, seed = 1, quiet = TRUE
    ))
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    expect_invisible(plot(cv))
})

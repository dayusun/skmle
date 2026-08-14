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

    fit <- kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x, h = h)

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
        xv <- cbind(1, d$x$x[k, 1][s])
        A <- A + crossprod(xv, xv * w)
        b <- b + as.numeric((w * d$y$y[j]) %*% xv)
    }
    expect_equal(unname(coef(fit)), as.numeric(solve(A, b)), tolerance = 1e-10)
})

test_that("the half kernel drops every pair with the covariate after the response", {
    set.seed(21)
    d <- sim_async_data(n = 60)
    h <- 0.25

    full <- kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x, h = h)
    half <- kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
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
        coef(kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x, h = h))[2]
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

    fit <- kee_async(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
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
        xv <- cbind(1, d$x$x[k, 1][s])
        U <- U + as.numeric((w * (d$y$y[j] - plogis(as.numeric(xv %*% bb)))) %*% xv)
    }
    expect_lt(max(abs(U)), 1e-6)

    yi <- d$y$id
    expect_error(kee_async(yi, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x, h = -1), "positive")
    expect_error(
        kee_async(yi, d$y$time, d$y$y, d$x$id[-1], d$x$time, d$x$x, h = 0.4),
        "one entry per row"
    )
    expect_error(
        kee_async(yi, d$y$time, d$y$y[-1], d$x$id, d$x$time, d$x$x, h = 0.4),
        "equal length"
    )
})

## ---------------------------------------------------------------- dependent

test_that("kee_async_td matches an independent transcription at each target time", {
    set.seed(13)
    d <- sim_async_data(n = 120, beta = function(tt) cbind(0.5, 1 + tt))
    times <- c(0.3, 0.5, 0.7)
    h <- 0.3
    W <- function(z) pmax(0.75 * (1 - z^2), 0)

    fit <- kee_async_td(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
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
            xv <- cbind(1, d$x$x[k, 1])
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
    fit <- kee_async_td(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
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
    a <- list(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x)

    expect_error(do.call(kee_async_td, c(a, list(times = 0.5, h = 0))), "positive")
    expect_error(do.call(kee_async_td, c(a, list(times = 0.5, h = 0.2, h2 = -1))), "positive")
    expect_error(do.call(kee_async_td, c(a, list(times = numeric(0), h = 0.2))), "non-empty")

    # A target time far outside the data has no weight on either side.
    expect_error(
        do.call(kee_async_td, c(a, list(times = 50, h = 0.05))),
        "no target time has data"
    )

    fit <- do.call(kee_async_td, c(a, list(times = c(0.5, 50), h = 0.05)))
    expect_identical(fit$convergence, c(0L, 3L))
    expect_true(all(is.na(fit$coefficients[2, ])))
})

test_that("print and plot methods work on a kee_td fit", {
    set.seed(29)
    d <- sim_async_data(n = 60)
    fit <- kee_async_td(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x,
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
    a <- list(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x)

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
        xv <- cbind(1, d$x$x[k, 1])
        A <- A + crossprod(xv, xv * ww)
        b <- b + as.numeric((ww * d$y$y[j]) %*% xv)
    }
    expect_equal(unname(coef(half)[1, ]), as.numeric(solve(A, b)), tolerance = 1e-10)
})

test_that("one_sided is validated", {
    set.seed(3)
    d <- sim_async_data(n = 30)
    a <- list(d$y$id, d$y$time, d$y$y, d$x$id, d$x$time, d$x$x)
    expect_error(do.call(kee_async, c(a, list(h = 0.3, one_sided = "yes"))), "TRUE or FALSE")
    expect_error(
        do.call(kee_async_td, c(a, list(times = 0.5, h = 0.3, one_sided = NA))),
        "TRUE or FALSE"
    )
})

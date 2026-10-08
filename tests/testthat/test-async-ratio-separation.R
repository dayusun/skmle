## The pairs correction reweights covariate occasions by an estimated
## response/covariate intensity ratio. When the covariate is collected in a
## window the responses barely reach (a marker assayed in one study phase only),
## that logistic fit separates; the ratio must stay finite and the fit must run.

sep_data <- function(n = 150, seed = 11) {
    set.seed(seed)
    # responses across [0, 17]; the covariate drawn in [18, 20] only, except one
    # stray early draw (the configuration of a late-phase assay)
    y <- do.call(rbind, lapply(seq_len(n), function(i) {
        data.frame(id = i, time = sort(runif(8, 0, 17)))
    }))
    x <- do.call(rbind, lapply(seq_len(n), function(i) {
        tt <- c(runif(2, 18, 20), if (i == 1L) 1)
        data.frame(id = i, time = tt, x = rnorm(1) + rnorm(length(tt), sd = 0.3))
    }))
    y$y <- 1 + 0.5 * x$x[match(y$id, x$id)] + rnorm(nrow(y))
    list(y = y, x = x)
}

test_that("separated observation schedules give a finite intensity ratio", {
    d <- sep_data()
    dd <- skmle:::.async_parse(y ~ x, d$y, d$x, "id", "time", TRUE)
    r <- suppressWarnings(skmle:::.async_intensity_ratio(dd))
    expect_true(all(is.finite(r)))
    expect_true(all(r > 0))
})

test_that("kee_async_cv with correct = 'pairs' runs when the schedules separate", {
    d <- sep_data()
    fit <- suppressWarnings(kee_async_cv(d$y, d$x, y ~ x, id = id, time = time,
                                         method = "variance", correct = "pairs",
                                         quiet = TRUE))
    expect_true(all(is.finite(coef(fit$fit))))
})

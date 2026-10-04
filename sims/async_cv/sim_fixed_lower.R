## Slim run for the undersmoothing study: per-h full-data slope and sandwich
## SE, and the bandwidth chosen by the fixed criterion at h_eval = lower end of
## the default grid (5-fold CV over 5 random splits, argmin). Output has the
## same structure as sim_acv4.R (picks holds only "fixed_lower") so that
## undersmooth.R reads either.
## Usage: Rscript sim_fixed_lower.R <n> <link> <first_rep> <n_reps> <out.rds> [grid_min n_grid]
## The default grid is 48 log-spaced points from 0.005 to 1.5.
args <- commandArgs(TRUE)
n <- as.integer(args[1]); link <- args[2]
first <- as.integer(args[3]); nrep <- as.integer(args[4]); out <- args[5]
suppressMessages(library(skmle))
ns <- asNamespace("skmle")
parse_ <- get(".async_parse", ns); subset_ <- get(".async_subset", ns)
fit_ <- get(".async_fit_ti", ns); collapse_ <- get(".async_collapse", ns)
lk <- c(identity = 0L, log = 1L, logistic = 2L)[[link]]
beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
gmin <- if (length(args) >= 7) as.numeric(args[6]) else 0.005
ng <- if (length(args) >= 7) as.integer(args[7]) else 48L
grid <- exp(seq(log(gmin), log(1.5), length.out = ng)); G <- length(grid)
K <- 5L; R <- 5L
kappa <- function(eta) if (lk == 0L) eta^2 else if (lk == 1L) exp(eta) else log1p(exp(-abs(eta))) + pmax(eta, 0)
cterm <- function(eta) if (lk == 0L) 2 * eta else eta

one <- function(rep_id) {
  set.seed(1e5 * lk + 1e3 * n + rep_id)
  d <- sim_async_data(n = n, beta = beta0, link = link)
  p <- parse_(y ~ x, d$y, d$x, "id", "time", TRUE)
  iqr <- diff(quantile(c(p$yt, p$xt), c(.25, .75), names = FALSE))
  lower <- 2 * iqr * n^(-0.7)
  b2 <- se2 <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    f <- tryCatch(suppressWarnings(kee_async(d$y, d$x, y ~ x, id = id, time = time,
                                             h = grid[g], link = link)),
                  error = function(e) NULL)
    if (is.null(f) || f$convergence == 1L) next
    b <- coef(f); s <- sqrt(diag(vcov(f)))
    if (all(is.finite(b)) && all(is.finite(s))) { b2[g] <- b[2]; se2[g] <- s[2] }
  }
  ac <- collapse_(p, lower, FALSE)
  om <- ac$omega / sum(ac$omega); cc <- ac$c / sum(ac$omega)
  folds <- matrix(replicate(R, sample(rep_len(seq_len(K), p$n))), nrow = p$n)
  acc <- matrix(0, G, R); bad <- matrix(FALSE, G, R)
  for (r in seq_len(R)) for (f in seq_len(K)) {
    te <- which(folds[, r] == f); tr <- subset_(p, which(folds[, r] != f))
    rows <- which(p$xi %in% te); Xr <- p$X[rows, , drop = FALSE]
    for (g in seq_len(G)) {
      ft <- fit_(tr, grid[g], FALSE, lk, 50L, 1e-8)
      if (is.null(ft) || ft$convergence == 1L || !all(is.finite(ft$beta))) { bad[g, r] <- TRUE; next }
      eta <- as.numeric(Xr %*% ft$beta)
      acc[g, r] <- acc[g, r] + sum(om[rows] * kappa(eta) - cc[rows] * cterm(eta))
    }
  }
  acc[bad] <- NA
  cv <- apply(acc, 1, function(v) if (all(is.na(v))) NA else mean(v, na.rm = TRUE))
  pick <- if (all(!is.finite(cv))) NA_integer_ else which.min(replace(cv, !is.finite(cv), Inf))
  list(rep = rep_id, lower = lower, b2 = b2, se2 = se2,
       picks = c(fixed_lower = pick), nfail_h = rowSums(bad))
}
res <- lapply(first:(first + nrep - 1L), function(r)
  tryCatch(one(r), error = function(e) list(rep = r, error = conditionMessage(e))))
saveRDS(list(n = n, link = link, grid = grid, beta0 = beta0, res = res), out)
cat("wrote", out, "\n")

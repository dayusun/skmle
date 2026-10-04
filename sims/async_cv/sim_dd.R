## Undersmoothing study with a data-driven bandwidth grid.
## Per data set: grid = 100 log-spaced points from h_min to h_max, the smallest
## and largest within-subject lag |T_ij - R_ik| (exact zeros excluded; the
## estimator here is two-sided, so all lags count). h_hat = argmin of the
## fixed criterion at h_eval = 2 IQR n^-0.7 (5-fold CV over 5 random splits).
## A candidate whose fit fails in any (split, fold) is excluded from CV.
## Full-data slope and sandwich SE are stored at every per-data-set grid
## point (for h_ci = c * h_hat, snapped there) and at a common reference grid
## (130 log-spaced points, 1e-8 to 1) for the yardstick and per-h diagnostics.
## A fit that errors, does not converge, or gives a non-finite estimate or SE
## is recorded as a failure (NA).
## Usage: Rscript sim_dd.R <n> <link> <first_rep> <n_reps> <out.rds>
args <- commandArgs(TRUE)
n <- as.integer(args[1]); link <- args[2]
first <- as.integer(args[3]); nrep <- as.integer(args[4]); out <- args[5]
suppressMessages(library(skmle))
ns <- asNamespace("skmle")
parse_ <- get(".async_parse", ns); subset_ <- get(".async_subset", ns)
fit_ <- get(".async_fit_ti", ns); collapse_ <- get(".async_collapse", ns)
sandwich_fit <- function(d, h) {
  f <- tryCatch(suppressWarnings(kee_async(d$y, d$x, y ~ x, id = id, time = time,
                                           h = h, link = link)), error = function(e) NULL)
  if (is.null(f) || f$convergence == 1L) return(c(NA_real_, NA_real_))
  b <- coef(f)[2]; s <- sqrt(diag(vcov(f)))[2]
  if (is.finite(b) && is.finite(s) && s > 0) c(b, s) else c(NA_real_, NA_real_)
}
lk <- c(identity = 0L, log = 1L, logistic = 2L)[[link]]
beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
NG <- 100L; ref <- exp(seq(log(1e-8), log(1), length.out = 130L))
K <- 5L; R <- 5L
kappa <- function(eta) if (lk == 0L) eta^2 else if (lk == 1L) exp(eta) else log1p(exp(-abs(eta))) + pmax(eta, 0)
cterm <- function(eta) if (lk == 0L) 2 * eta else eta

one <- function(rep_id) {
  set.seed(1e5 * lk + 1e3 * n + rep_id)
  d <- sim_async_data(n = n, beta = beta0, link = link)
  p <- parse_(y ~ x, d$y, d$x, "id", "time", TRUE)
  xs <- split(p$xt, p$xi); ys <- split(p$yt, p$yi)
  ids <- intersect(names(xs), names(ys))
  lg <- abs(unlist(lapply(ids, function(i) outer(ys[[i]], xs[[i]], "-")), use.names = FALSE))
  lg <- lg[lg > 0]
  hmin <- min(lg); hmax <- max(lg)
  grid <- exp(seq(log(hmin), log(hmax), length.out = NG))
  iqr <- diff(quantile(c(p$yt, p$xt), c(.25, .75), names = FALSE))
  lower <- 2 * iqr * n^(-0.7)
  own <- vapply(grid, function(h) sandwich_fit(d, h), numeric(2))
  refv <- vapply(ref, function(h) if (h < hmin / 2) c(NA_real_, NA_real_) else sandwich_fit(d, h), numeric(2))
  ac <- collapse_(p, lower, FALSE)
  om <- ac$omega / sum(ac$omega); cc <- ac$c / sum(ac$omega)
  folds <- matrix(replicate(R, sample(rep_len(seq_len(K), p$n))), nrow = p$n)
  acc <- matrix(0, NG, R); bad <- matrix(FALSE, NG, R)
  for (r in seq_len(R)) for (f in seq_len(K)) {
    te <- which(folds[, r] == f); tr <- subset_(p, which(folds[, r] != f))
    rows <- which(p$xi %in% te); Xr <- p$X[rows, , drop = FALSE]
    for (g in seq_len(NG)) {
      if (bad[g, r]) next
      ft <- tryCatch(fit_(tr, grid[g], FALSE, lk, 50L, 1e-8), error = function(e) NULL)
      if (is.null(ft) || ft$convergence == 1L || !all(is.finite(ft$beta))) { bad[g, r] <- TRUE; next }
      eta <- as.numeric(Xr %*% ft$beta)
      acc[g, r] <- acc[g, r] + sum(om[rows] * kappa(eta) - cc[rows] * cterm(eta))
    }
  }
  excluded <- rowSums(bad) > 0
  cv <- rowMeans(acc); cv[excluded] <- NA
  pick <- if (all(!is.finite(cv))) NA_integer_ else which.min(replace(cv, !is.finite(cv), Inf))
  list(rep = rep_id, hmin = hmin, hmax = hmax, npairs = length(lg), lower = lower,
       grid = grid, b_own = own[1, ], se_own = own[2, ], b_ref = refv[1, ], se_ref = refv[2, ],
       pick = pick, cv_excluded = excluded)
}
res <- lapply(first:(first + nrep - 1L), function(r)
  tryCatch(one(r), error = function(e) list(rep = r, error = conditionMessage(e))))
saveRDS(list(n = n, link = link, ref = ref, beta0 = beta0, res = res), out)
cat("wrote", out, "\n")

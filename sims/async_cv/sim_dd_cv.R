## Plain K-fold cross-validation over subjects on the data-driven grid: each
## candidate h is fitted on the training subjects and scored on the held-out
## pairs with its own kernel W_h, normalised by the total held-out weight at
## that h (pooled over folds, averaged over splits). No other bandwidth,
## rate or constant enters the loss.
## Loss per pair: squared error (identity) or negative Bernoulli
## log-likelihood (logistic). K = 5 folds x R = 5 splits, argmin. A candidate is
## excluded if its fit fails in any (split, fold), or if no held-out pair falls
## in its window in some (split, fold).
## The data set, its grid and the full-data fits are those of sim_dd.R (same
## seed, same grid construction); the aggregator joins on (link, n, rep) and
## checks the grid end points.
## Usage: Rscript sim_dd_cv.R <n> <link> <first_rep> <n_reps> <out.rds>
args <- commandArgs(TRUE)
n <- as.integer(args[1]); link <- args[2]
first <- as.integer(args[3]); nrep <- as.integer(args[4]); out <- args[5]
suppressMessages(library(skmle))
ns <- asNamespace("skmle")
parse_ <- get(".async_parse", ns); subset_ <- get(".async_subset", ns)
fit_ <- get(".async_fit_ti", ns); collapse_ <- get(".async_collapse", ns)
lk <- c(identity = 0L, log = 1L, logistic = 2L)[[link]]
beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
NG <- 100L; K <- 5L; R <- 5L
kappa <- function(eta) if (lk == 0L) eta^2 else if (lk == 1L) exp(eta) else log1p(exp(-abs(eta))) + pmax(eta, 0)
cterm <- function(eta) if (lk == 0L) 2 * eta else eta
# Loss of one fold at one evaluation window: numerator and denominator.
fold_loss <- function(ac, eta) c(sum(ac$omega * kappa(eta)) - sum(ac$c * cterm(eta)) +
                                   if (lk == 0L) ac$sy2 else 0, sum(ac$omega))

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
  folds <- matrix(replicate(R, sample(rep_len(seq_len(K), p$n))), nrow = p$n)
  fitbad <- matrix(FALSE, NG, R)
  numM <- denM <- matrix(0, NG, R); nopair <- matrix(FALSE, NG, R)
  for (r in seq_len(R)) for (f in seq_len(K)) {
    te <- which(folds[, r] == f); tr <- subset_(p, which(folds[, r] != f))
    dte <- subset_(p, te)
    for (g in seq_len(NG)) {
      if (fitbad[g, r]) next
      ft <- tryCatch(fit_(tr, grid[g], FALSE, lk, 50L, 1e-8), error = function(e) NULL)
      if (is.null(ft) || ft$convergence == 1L || !all(is.finite(ft$beta))) { fitbad[g, r] <- TRUE; next }
      ac <- collapse_(dte, grid[g], FALSE)
      if (is.null(ac)) { nopair[g, r] <- TRUE; next }
      v <- fold_loss(ac, as.numeric(dte$X %*% ft$beta)); numM[g, r] <- numM[g, r] + v[1]; denM[g, r] <- denM[g, r] + v[2]
    }
  }
  excl_fit <- rowSums(fitbad) > 0
  excl_m <- excl_fit | rowSums(nopair) > 0
  cvM <- rowMeans(numM / denM); cvM[excl_m] <- NA
  pickM <- if (all(!is.finite(cvM))) NA_integer_ else which.min(replace(cvM, !is.finite(cvM), Inf))
  list(rep = rep_id, hmin = hmin, hmax = hmax, pick_matched = pickM,
       excl_fit = excl_fit, excl_matched = excl_m)
}
res <- lapply(first:(first + nrep - 1L), function(r)
  tryCatch(one(r), error = function(e) list(rep = r, error = conditionMessage(e))))
saveRDS(list(n = n, link = link, beta0 = beta0, res = res), out)
cat("wrote", out, "\n")

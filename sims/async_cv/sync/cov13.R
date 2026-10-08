## Plug-in minimum-variance bandwidth for the pair-moment estimator:
##   h_V = argmin_grid tr{ Sigma_X V-hat(h) },  Sigma_X = sum r(S_k) X_k X_k' / sum r(S_k),
## against the sync-loss CV of cov12 on the same data and seeds. The whole grid
## of (beta-hat(h), SE(h)) is kept so the aggregate can locate the oracle h.
source("cov12.R")

## Grid floor on the full sample: at least ceil(sqrt(n)) subjects with a
## response-covariate pair and with a covariate-covariate pair inside the kernel.
floor13 <- function(P, Q, n) {
  f <- function(u, sid) { mi <- tapply(u, factor(sid, levels = seq_len(n)), min); mi[is.na(mi)] <- Inf
    m <- sort(mi); m[min(ceiling(sqrt(n)), length(m))] }
  max(f(abs(P[, "u"]), P[, "sid"]), f(abs(Q[, "v"]), Q[, "sid"])) * (1 + 1e-8)
}

select13 <- function(e, link, start, r) {
  P <- e$P; Q <- e$Q; n <- e$p$n
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), floor13(P, Q, n))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG))
  Sx <- crossprod(e$X * sqrt(r)) / sum(r)
  p <- ncol(e$X); EST <- SE <- matrix(NA_real_, NG, p); crit <- rep(Inf, NG); st <- start
  for (g in rev(seq_len(NG))) {
    ft <- tryCatch(fit_pmw(e, P, Q, grid[g], link, st, r), error = function(err) NULL)
    if (is.null(ft) || !ft$conv || !all(is.finite(ft$V))) next
    st <- ft$beta; EST[g, ] <- ft$beta; SE[g, ] <- sqrt(diag(ft$V)); crit[g] <- sum(Sx * ft$V)
  }
  g <- which.min(crit)
  list(h = grid[g], lo = lo, g = g, grid = grid, EST = EST, SE = SE, crit = crit)
}

one_rep13 <- function(n, link, rep_id, gen) {
  beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); r <- ratio_hat(e$p)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(length(nm)))$beta,
                 error = function(err) numeric(length(nm)))
  hs <- select12(e, link, st, r)
  f <- tryCatch(fit_pmw(e, e$P, e$Q, hs[["h"]], link, st, r), error = function(err) NULL)
  out <- row7("pm_sat", hs[["h"]], TRUE, f, nm); out$h_lo <- hs[["lo"]]
  s <- select13(e, link, st, r)
  ok <- is.finite(s$crit[s$g])
  out <- rbind(out, data.frame(h_type = "pm_var", h = s$h, corr = TRUE, coef = nm,
    est = if (ok) s$EST[s$g, ] else NA_real_, se = if (ok) s$SE[s$g, ] else NA_real_, conv = ok, h_lo = s$lo))
  list(tab = cbind(gen = gen, link = link, n = n, rep = rep_id, truth = beta0[match(out$coef, nm)], out),
       grid = list(gen = gen, link = link, n = n, rep = rep_id, truth = beta0, coef = nm,
                   grid = s$grid, EST = s$EST, SE = s$SE, crit = s$crit, g = s$g))
}

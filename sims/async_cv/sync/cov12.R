## Pair-moment with sync CV and the saturated-loss rule (the package's default
## path, verified identical in eqv12.R). Same seeds as job 8478458.
source("cov11.R")
sat_ref <- function(y, link) if (link == "log") ifelse(y > 0, y - y * log(y), 0) else 0 * y
select12 <- function(e, link, start, r) {
  P <- e$P; Q <- e$Q; pl <- loss_of(link); folds <- ((e$ids - 1L) %% K) + 1L
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), floor9(P, Q, e$p$n, folds))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG)); IC <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) { te <- which(folds == f)
    H <- P[P[, "sid"] %in% te, , drop = FALSE]; Ptr <- P[!P[, "sid"] %in% te, , drop = FALSE]; Qtr <- Q[!Q[, "sid"] %in% te, , drop = FALSE]
    HX <- e$X[H[, "k"], , drop = FALSE]; QR <- qr(cbind(1, e$gb[e$block_of(abs(H[, "u"])), , drop = FALSE])); s0 <- mean(sat_ref(H[, "y"], link)); st <- start
    for (g in rev(seq_len(NG))) { ft <- tryCatch(fit_pmw(e, Ptr, Qtr, grid[g], link, st, r), error = function(err) NULL)
      if (is.null(ft) || !ft$conv) next; st <- ft$beta
      ic <- qr.coef(QR, pl(H[, "y"], as.vector(HX %*% ft$beta)))[1]; if (is.finite(ic) && ic >= s0) IC[g, f] <- ic } }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; c(h = grid[which.min(cv)], lo = lo)
}
one_rep12 <- function(n, link, rep_id, gen) {
  beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); r <- ratio_hat(e$p)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(length(nm)))$beta,
                 error = function(err) numeric(length(nm)))
  hs <- select12(e, link, st, r)
  f <- tryCatch(fit_pmw(e, e$P, e$Q, hs[["h"]], link, st, r), error = function(err) NULL)
  out <- row7("pm_sat", hs[["h"]], TRUE, f, nm)
  cbind(gen = gen, link = link, n = n, rep = rep_id, truth = beta0[match(out$coef, nm)], out, h_lo = hs[["lo"]])
}

## Pair-moment estimator with sync-loss CV. No distributional model for x:
## fits use ordered within-subject covariate pairs weighted by the estimated
## intensity ratio lambda_T/lambda_S; Gamma-hat enters only the CV criterion.
source("cov8.R")

floor9 <- function(P, Q, n, folds) {
  f1 <- floor7(P, n, folds)
  Qm <- cbind(u = Q[, "v"], sid = Q[, "sid"]); f2 <- floor7(Qm, n, folds)
  max(f1, f2)
}

select9 <- function(e, link, start, r) {
  P <- e$P; Q <- e$Q; pl <- loss_of(link); folds <- ((e$ids - 1L) %% K) + 1L
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), floor9(P, Q, e$p$n, folds))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG))
  IC <- ICm <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) {
    te <- which(folds == f)
    H <- P[P[, "sid"] %in% te, , drop = FALSE]
    Ptr <- P[!P[, "sid"] %in% te, , drop = FALSE]; Qtr <- Q[!Q[, "sid"] %in% te, , drop = FALSE]
    CDtr <- cond_design(e, Ptr); HX <- e$X[H[, "k"], , drop = FALSE]
    QR <- qr(cbind(1, e$gb[e$block_of(abs(H[, "u"])), , drop = FALSE]))
    for (m in c("plain", "pm")) { st <- start
      for (g in rev(seq_len(NG))) {
        ft <- tryCatch(if (m == "plain") fit7(e, Ptr, grid[g], link, FALSE, st, CDtr) else
          fit_pmw(e, Ptr, Qtr, grid[g], link, st, r), error = function(err) NULL)
        if (is.null(ft) || !ft$conv) next
        st <- ft$beta
        v <- qr.coef(QR, pl(H[, "y"], as.vector(HX %*% ft$beta)))[1]
        if (m == "pm") ICm[g, f] <- v else IC[g, f] <- v } }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; cvm <- rowMeans(ICm); cvm[!is.finite(cvm)] <- Inf
  c(cv = grid[which.min(cv)], pm = grid[which.min(cvm)], lo = lo)
}

GENS <- list(gauss = function(n, b, l) sim_multi(n, b, l),
             skew = function(n, b, l) sim_skew2(n, b, l, early = FALSE),
             skew_early = function(n, b, l) sim_skew2(n, b, l, early = TRUE))

one_rep9 <- function(n, link, rep_id, gen) {
  beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); r <- ratio_hat(e$p)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(length(nm)))$beta,
                 error = function(err) numeric(length(nm)))
  hs <- select9(e, link, st, r); out <- NULL
  r1 <- tryCatch(fit7(e, e$P, hs[["cv"]], link, FALSE, st), error = function(err) NULL)
  out <- rbind(out, row7("cv", hs[["cv"]], FALSE, r1, nm))
  r2 <- tryCatch(fit_pmw(e, e$P, e$Q, hs[["pm"]], link, st, r), error = function(err) NULL)
  out <- rbind(out, row7("pm", hs[["pm"]], TRUE, r2, nm))
  r3 <- tryCatch(fit7(e, lvcf7(e), 1e6, link, FALSE, st), error = function(err) NULL)
  out <- rbind(out, row7("lvcf", NA, FALSE, r3, nm))
  cbind(gen = gen, link = link, n = n, rep = rep_id, truth = beta0[match(out$coef, nm)], out, h_lo = hs[["lo"]])
}

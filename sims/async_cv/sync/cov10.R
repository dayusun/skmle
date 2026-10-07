## Cross-fitted bandwidth for the pair-moment estimator. For each outer fold k,
## h_k is chosen by sync-loss CV on the subjects outside fold k; subjects in
## fold k enter the estimating equation at h_k. Given the other folds, h_k is
## independent of fold k, so every fold's moment has mean zero at beta_0.
source("cov9.R")
NG_IN <- 30L

## Sync-loss CV of the pair-moment fit using only subjects `sub`.
select_pm_sub <- function(e, link, start, r, sub, ng = NG_IN) {
  P <- e$P[e$P[, "sid"] %in% sub, , drop = FALSE]; Q <- e$Q[e$Q[, "sid"] %in% sub, , drop = FALSE]
  pl <- loss_of(link); inner <- integer(e$p$n); inner[sub] <- ((seq_along(sub) - 1L) %% K) + 1L
  ## floor on the inner folds, computed on the subset
  mi <- function(M, col) { m <- tapply(abs(M[, col]), factor(M[, "sid"], levels = sub), min); m[is.na(m)] <- Inf; m }
  m1 <- mi(P, "u"); m2 <- mi(Q, "v"); fi <- inner[sub]
  lo <- max(vapply(seq_len(K), function(f) { k <- ceiling(sqrt(sum(fi != f)))
    max(sort(m1[fi != f])[k], sort(m2[fi != f])[k]) }, 0)) * (1 + 1e-8)
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), lo)
  grid <- exp(seq(log(lo), log(max(lg)), length.out = ng))
  IC <- matrix(NA_real_, ng, K)
  for (f in seq_len(K)) {
    te <- sub[fi == f]
    H <- P[P[, "sid"] %in% te, , drop = FALSE]; Ptr <- P[!P[, "sid"] %in% te, , drop = FALSE]
    Qtr <- Q[!Q[, "sid"] %in% te, , drop = FALSE]; HX <- e$X[H[, "k"], , drop = FALSE]
    QR <- qr(cbind(1, e$gb[e$block_of(abs(H[, "u"])), , drop = FALSE])); st <- start
    for (g in rev(seq_len(ng))) {
      ft <- tryCatch(fit_pmw(e, Ptr, Qtr, grid[g], link, st, r), error = function(err) NULL)
      if (is.null(ft) || !ft$conv) next
      st <- ft$beta; IC[g, f] <- qr.coef(QR, pl(H[, "y"], as.vector(HX %*% ft$beta)))[1] }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf
  grid[which.min(cv)]
}

## Pair-moment fit where subjects in outer fold k use bandwidth hk[k].
fit_pm_cf <- function(e, hk, folds, link, start, r) {
  Kf <- length(hk); n <- e$p$n; Lk <- LINKS[[link]]
  full <- function(M, g) { out <- matrix(0, n, ncol(M)); rr <- rowsum(M, g); out[as.integer(rownames(rr)), ] <- rr; out }
  parts <- lapply(seq_len(Kf), function(k) {
    sub <- which(folds == k)
    P <- e$P[e$P[, "sid"] %in% sub, , drop = FALSE]; Q <- e$Q[e$Q[, "sid"] %in% sub, , drop = FALSE]
    w <- epan(P[, "u"], hk[k]); kp <- w > 0; P <- P[kp, , drop = FALSE]; w <- w[kp]
    v <- epan(Q[, "v"], hk[k]) * r[Q[, "b"]]; kq <- v > 0; Q <- Q[kq, , drop = FALSE]; v <- v[kq]
    list(P = P, Q = Q, w = w, v = v, WA = sum(w), WB = sum(v), share = length(sub) / n,
         Ai = full(w * e$X[P[, "k"], , drop = FALSE] * P[, "y"], P[, "sid"]),
         wa = full(cbind(w), P[, "sid"])[, 1], wb = full(cbind(v), Q[, "sid"])[, 1],
         Xa = e$X[Q[, "a"], , drop = FALSE], Xb = e$X[Q[, "b"], , drop = FALSE]) })
  ef <- function(b) {
    U <- 0; J <- 0; psi <- 0
    for (pt in parts) {
      eta <- as.vector(pt$Xb %*% b); Bi <- full(pt$v * pt$Xa * Lk$mu(eta), pt$Q[, "sid"])
      A <- colSums(pt$Ai) / pt$WA; B <- colSums(Bi) / pt$WB
      U <- U + pt$share * (A - B)
      J <- J - pt$share * crossprod(pt$Xa * pt$v * Lk$d(eta), pt$Xb) / pt$WB
      psi <- psi + pt$share * ((pt$Ai - outer(pt$wa, A)) / pt$WA - (Bi - outer(pt$wb, B)) / pt$WB) }
    list(U = U, J = J, psi = psi) }
  b <- start; st <- Inf
  for (it in 1:50) { ee <- ef(b); st <- tryCatch(solve(ee$J, ee$U), error = function(err) NULL)
    if (is.null(st)) return(NULL); b <- b - st
    if (!all(is.finite(b))) return(NULL); if (max(abs(st)) < 1e-10) break }
  ee <- ef(b); Ji <- solve(ee$J)
  list(beta = b, V = Ji %*% crossprod(ee$psi) %*% t(Ji), conv = max(abs(st)) < 1e-6)
}

one_rep10 <- function(n, link, rep_id, gen) {
  beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); r <- ratio_hat(e$p)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(length(nm)))$beta,
                 error = function(err) numeric(length(nm)))
  folds <- ((e$ids - 1L) %% K) + 1L
  ## outer folds use a different assignment from the inner ones: shift by a seed-free permutation
  folds <- folds[order(((e$ids * 7919L) %% e$p$n))]
  hk <- vapply(seq_len(K), function(k) select_pm_sub(e, link, st, r, which(folds != k)), 0)
  f <- tryCatch(fit_pm_cf(e, hk, folds, link, st, r), error = function(err) NULL)
  out <- row7("pm_cf", exp(mean(log(hk))), TRUE, f, nm)
  cbind(gen = gen, link = link, n = n, rep = rep_id, truth = beta0[match(out$coef, nm)], out,
        h_min = min(hk), h_max = max(hk))
}

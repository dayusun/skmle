## cov5 with a data-driven lower end for the h grid: the smallest h at which, in every
## training fold, at least ceiling(sqrt(n_train)) subjects have a pair with |u| < h.
source("cov5.R")

h_floor <- function(P, ids, folds) {
  mi <- tapply(abs(P[, "u"]), factor(P[, "sid"], levels = ids), min)   # closest pair per subject
  max(vapply(seq_len(K), function(f) {
    m <- sort(mi[folds != f]); k <- ceiling(sqrt(sum(folds != f))); m[min(k, length(m))] }, 0)) * (1 + 1e-8)
}

select_g <- function(e, link, beta0) {
  P <- e$P; pl <- loss_of(link)
  folds <- ((e$ids - 1L) %% K) + 1L
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), h_floor(P, e$ids, folds))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG))
  IC <- ICc <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) {
    H <- P[P[, "sid"] %in% which(folds == f), , drop = FALSE]
    Ptr <- P[P[, "sid"] %in% which(folds != f), , drop = FALSE]
    QR <- qr(cbind(1, e$gam(abs(H[, "u"]))))
    for (corr in c(FALSE, TRUE)) { st <- beta0
      for (g in rev(seq_len(NG))) {
        ft <- tryCatch(fit_g(Ptr, grid[g], link, corr, e$gam, e$m, e$s2, st), error = function(e) NULL)
        if (is.null(ft) || !ft$conv) next
        st <- ft$beta
        v <- qr.coef(QR, pl(H[, "y"], ft$beta[1] + ft$beta[2] * H[, "x"]))[1]
        if (corr) ICc[g, f] <- v else IC[g, f] <- v
      } }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; cvc <- rowMeans(ICc); cvc[!is.finite(cvc)] <- Inf
  .last_lo <<- lo
  c(cv = grid[which.min(cv)], cvc = grid[which.min(cvc)], lo = lo)
}

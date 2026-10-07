## Subject bootstrap at the cv+corr bandwidth.
## Per data set: h-hat by sync-loss CV on plain fold fits; corrected fit at h-hat with sandwich SE.
## Per resample (subjects drawn with replacement, relabelled as distinct subjects):
##   fix: gamma-hat, mean, variance re-estimated; h held at the original h-hat
##   sel: the same, plus h re-selected by CV on the resample
source("cov3.R")

select_cv <- function(e, lk) {
  p <- e$p; P <- e$P; gam <- e$gam
  pl <- function(y, eta) if (lk == 0L) (y - eta)^2 else log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta
  lg <- abs(P[, "u"]); lg <- lg[lg > 0]
  grid <- exp(seq(log(min(lg)), log(max(lg)), length.out = NG))
  folds <- ((e$ids - 1L) %% K) + 1L
  IC <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) {
    tr <- subset_(p, which(folds != f)); H <- P[P[, "sid"] %in% which(folds == f), , drop = FALSE]
    QR <- qr(cbind(1, gam(abs(H[, "u"]))))
    for (g in seq_len(NG)) {
      ft <- tryCatch(fit_(tr, grid[g], FALSE, lk, 50L, 1e-8), error = function(e) NULL)
      if (is.null(ft) || ft$convergence == 1L || !all(is.finite(ft$beta))) next
      IC[g, f] <- qr.coef(QR, pl(H[, "y"], ft$beta[1] + ft$beta[2] * H[, "x"]))[1]
    }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; grid[which.min(cv)]
}

## Same pairs, same order as pairs_of(), without the per-subject which() scans.
pairs_of <- function(p, xv, ids) {
  ys <- split(seq_along(p$yi), factor(p$yi, levels = ids)); xs <- split(seq_along(p$xi), factor(p$xi, levels = ids))
  do.call(rbind, lapply(seq_along(ids), function(k) {
    yi <- ys[[k]]; xi <- xs[[k]]; if (!length(yi) || !length(xi)) return(NULL)
    cbind(u = as.vector(outer(p$yt[yi], p$xt[xi], "-")), y = rep(p$y[yi], times = length(xi)),
          x = rep(xv[xi], each = length(yi)), sid = ids[k]) }))
}

resample <- function(d) {
  u <- sort(unique(d$y$id)); s <- sample(u, length(u), replace = TRUE)
  f <- function(df) { r <- split(seq_len(nrow(df)), factor(df$id, levels = u))[as.character(s)]
    out <- df[unlist(r), , drop = FALSE]; out$id <- rep(seq_along(s), lengths(r)); out }
  list(y = f(d$y), x = f(d$x))
}

gen <- function(n, link, rep_id, scen) {
  lk <- c(identity = 0L, logistic = 2L)[[link]]
  beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
  if (scen == "base") { set.seed(1e5 * lk + 1e3 * n + rep_id); d <- sim_async_data(n = n, beta = beta0, link = link) }
  else { set.seed(1e7 * match(scen, c("nonstat", "merr", "skew")) + 1e5 * lk + 1e3 * n + rep_id); d <- sim_rob(n, beta0, link, scen) }
  list(d = d, lk = lk, beta0 = beta0)
}

one_rep_boot <- function(n, link, rep_id, scen, B) {
  g <- gen(n, link, rep_id, scen); lk <- g$lk; beta0 <- g$beta0
  e <- prep(g$d); h <- select_cv(e, lk)
  r <- fit_pairs(e$P, h, lk, TRUE, e$gam, e$m, e$s2, beta0)
  set.seed(1e8 + 1e5 * lk + 1e3 * n + rep_id + 7e6 * match(scen, c("base", "nonstat", "merr", "skew")))
  bf <- bs <- hb <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    eb <- tryCatch(prep(resample(g$d)), error = function(e) NULL); if (is.null(eb)) next
    f1 <- tryCatch(fit_pairs(eb$P, h, lk, TRUE, eb$gam, eb$m, eb$s2, r$beta), error = function(e) NULL)
    if (!is.null(f1) && f1$conv) bf[b] <- f1$beta[2]
    hb[b] <- tryCatch(select_cv(eb, lk), error = function(e) NA)
    if (is.finite(hb[b])) {
      f2 <- tryCatch(fit_pairs(eb$P, hb[b], lk, TRUE, eb$gam, eb$m, eb$s2, r$beta), error = function(e) NULL)
      if (!is.null(f2) && f2$conv) bs[b] <- f2$beta[2]
    }
  }
  q <- function(v) quantile(v, c(.025, .975), na.rm = TRUE, names = FALSE)
  data.frame(scen = scen, link = link, n = n, rep = rep_id, h = h, b1 = r$beta[2], se_sand = sqrt(r$V[2, 2]),
    conv = r$conv, B_fix = sum(is.finite(bf)), B_sel = sum(is.finite(bs)),
    se_fix = sd(bf, na.rm = TRUE), se_sel = sd(bs, na.rm = TRUE),
    pf_lo = q(bf)[1], pf_hi = q(bf)[2], ps_lo = q(bs)[1], ps_hi = q(bs)[2],
    hb_med = median(hb, na.rm = TRUE))
}

## Aggregate results13/: sync-loss CV (pm_sat) vs plug-in minimum-variance h (pm_var),
## plus the oracle fixed h located from the stored grids.
L <- unlist(lapply(Sys.glob("results13/*.rds"), readRDS), recursive = FALSE)
err <- vapply(L, is.character, TRUE); cat("errors:", sum(err), "\n"); L <- L[!err]
x <- do.call(rbind, lapply(L, `[[`, "tab")); x <- x[x$coef != "(Intercept)", ]
x$arm <- paste0(x$h_type, "+corr")
summ <- function(d) {
  ok <- d$conv & is.finite(d$est) & is.finite(d$se); R <- sum(ok); b <- d$est[ok]; s <- d$se[ok]; tr <- d$truth[1]
  k <- sum(abs(b - tr) <= qnorm(.975) * s); bt <- binom.test(k, R, 0.95)
  data.frame(R = R, miss = nrow(d) - R, med_h = median(d$h[ok]), bias = mean(b) - tr, bias_se = sd(b) / sqrt(R),
    bias_ESE = (mean(b) - tr) / sd(b), ESE = sd(b), SE_ESE = median(s) / sd(b), RMSE = sqrt(mean((b - tr)^2)),
    CP = 100 * k / R, CP_lo = 100 * bt$conf.int[1], CP_hi = 100 * bt$conf.int[2], p = bt$p.value)
}
t <- do.call(rbind, lapply(split(x, paste(x$gen, x$link, x$n, x$arm, x$coef)), function(d)
  cbind(d[1, c("gen", "link", "n", "arm", "coef")], summ(d))))

## Oracle: common log-grid per cell, nearest stored h per rep; h* minimises the summed MSE of x1, x2, z.
G <- lapply(L, `[[`, "grid"); key <- vapply(G, function(g) paste(g$gen, g$link, g$n), "")
orc <- do.call(rbind, lapply(split(G, key), function(gs) {
  lo <- quantile(vapply(gs, function(g) min(g$grid), 0), .9); hi <- quantile(vapply(gs, function(g) max(g$grid), 0), .1)
  hc <- exp(seq(log(lo), log(hi), length.out = 40)); nm <- gs[[1]]$coef; tr <- gs[[1]]$truth; cc <- which(nm != "(Intercept)")
  pick <- function(g) { ix <- vapply(hc, function(h) which.min(abs(log(g$grid) - log(h))), 0L); list(E = g$EST[ix, cc, drop = FALSE], S = g$SE[ix, cc, drop = FALSE]) }
  P <- lapply(gs, pick); E <- simplify2array(lapply(P, `[[`, "E")); S <- simplify2array(lapply(P, `[[`, "S"))  # [h, coef, rep]
  mse <- apply(E, 1, function(m) sum(rowMeans((m - tr[cc])^2, na.rm = TRUE))); gstar <- which.min(mse)
  hv <- vapply(gs, function(g) g$grid[g$g], 0)
  do.call(rbind, lapply(seq_along(cc), function(j) { b <- E[gstar, j, ]; s <- S[gstar, j, ]; ok <- is.finite(b) & is.finite(s); b <- b[ok]; s <- s[ok]
    k <- sum(abs(b - tr[cc[j]]) <= qnorm(.975) * s); bt <- binom.test(k, length(b), .95)
    data.frame(gen = gs[[1]]$gen, link = gs[[1]]$link, n = gs[[1]]$n, arm = "oracle", coef = nm[cc[j]], R = length(b), miss = length(gs) - length(b),
      med_h = hc[gstar], bias = mean(b) - tr[cc[j]], bias_se = sd(b) / sqrt(length(b)), bias_ESE = (mean(b) - tr[cc[j]]) / sd(b), ESE = sd(b),
      SE_ESE = median(s) / sd(b), RMSE = sqrt(mean((b - tr[cc[j]])^2)), CP = 100 * k / length(b), CP_lo = 100 * bt$conf.int[1], CP_hi = 100 * bt$conf.int[2], p = bt$p.value,
      hV_q10 = quantile(hv, .1), hV_q50 = median(hv), hV_q90 = quantile(hv, .9), hV_within2 = mean(abs(log(hv / hc[gstar])) <= log(2))) }))
}))
t <- merge(t, orc[, c("gen", "link", "n", "coef", "hV_q10", "hV_q50", "hV_q90", "hV_within2")], all.x = TRUE)
t <- rbind(t, orc[, names(t)])
t$q <- p.adjust(t$p, "BH"); t$sig <- as.character(cut(t$q, c(-1, .001, .01, .05, 2), c("***", "**", "*", "")))
t <- t[order(t$gen, t$link, t$n, t$coef, t$arm), ]; rownames(t) <- NULL
saveRDS(t, "summary13.rds"); write.csv(t, "summary13.csv", row.names = FALSE)
options(width = 250); print(t[, c("gen", "link", "n", "arm", "coef", "R", "miss", "med_h", "bias", "bias_ESE", "ESE", "SE_ESE", "RMSE", "CP", "q", "sig")], digits = 3)
cat("\nh_V vs oracle h*:\n"); print(unique(t[t$arm == "oracle", c("gen", "link", "n", "med_h", "hV_q10", "hV_q50", "hV_q90", "hV_within2")]), digits = 3)

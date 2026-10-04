## Aggregate plain matched CV on the data-driven grid (sim_dd_cv.R), joined to
## the full-data slope and sandwich SE on each data set's own grid and the
## common reference grid from sim_dd.R (same seeds, same data, same grid; the
## grid end points are checked replicate by replicate).
## The best fixed h is the reference-grid h with the smallest simulated MSE
## among points that fit in at least 95% of replicates: a yardstick that uses
## the known beta, not a selector.
## Usage: Rscript agg_dd_cv.R <cv dir> <sim_dd dir>
a <- commandArgs(TRUE); cvdir <- a[1]; dddir <- a[2]
z <- qnorm(0.975)
key <- function(o) paste(o$link, o$n)
cvs <- lapply(Sys.glob(file.path(cvdir, "*_n*_r*.rds")), readRDS)
dds <- lapply(Sys.glob(file.path(dddir, "*_n*_r*.rds")), readRDS)
ref <- dds[[1]]$ref
tab <- list()
for (k in unique(sapply(cvs, key))) {
  C <- unlist(lapply(cvs[sapply(cvs, key) == k], `[[`, "res"), recursive = FALSE)
  D <- unlist(lapply(dds[sapply(dds, key) == k], `[[`, "res"), recursive = FALSE)
  o <- cvs[[which(sapply(cvs, key) == k)[1]]]; link <- o$link; n <- o$n; b0 <- o$beta0[2]
  C <- C[!duplicated(sapply(C, `[[`, "rep"))]; D <- D[!duplicated(sapply(D, `[[`, "rep"))]
  cerr <- sum(sapply(C, function(r) !is.null(r$error))); C <- C[sapply(C, function(r) is.null(r$error))]
  D <- D[match(sapply(C, `[[`, "rep"), sapply(D, `[[`, "rep"))]
  mismatch <- sum(mapply(function(c, d) is.null(d) || !isTRUE(all.equal(c(c$hmin, c$hmax), c(d$hmin, d$hmax))), C, D))
  B <- t(sapply(D, `[[`, "b_ref")); S <- t(sapply(D, `[[`, "se_ref"))
  okh <- colMeans(is.finite(B) & is.finite(S))
  mse_ref <- colMeans((B - b0)^2, na.rm = TRUE); mse_ref[okh < 0.95] <- NA
  g_best <- which.min(mse_ref)
  nr <- length(C)
  g <- sapply(C, `[[`, "pick_matched")
  hh <- est <- se <- rep(NA_real_, nr)
  for (i in seq_len(nr)) if (!is.na(g[i])) {
    hh[i] <- D[[i]]$grid[g[i]]; est[i] <- D[[i]]$b_own[g[i]]; se[i] <- D[[i]]$se_own[g[i]]
  }
  ok <- is.finite(est) & is.finite(se); e <- est[ok] - b0; m <- sum(ok)
  kc <- sum(abs(e) <= z * se[ok]); cp <- kc / m
  # MSE at the best fixed h, for the ratio, with its own MC SE.
  eb <- B[, g_best] - b0; eb <- eb[is.finite(eb)]
  tab[[length(tab) + 1]] <- data.frame(link = link, n = n, N = nr, rep_errors = cerr,
    grid_mismatch = mismatch, no_pick = sum(is.na(g)), fit_fail_at_hhat = sum(!is.na(g) & !ok),
    cand_excl_pct = 100 * mean(sapply(C, function(r) mean(r$excl_matched))),
    h_hat_q25 = quantile(hh, .25, na.rm = TRUE), h_hat_med = median(hh, na.rm = TRUE),
    h_hat_q75 = quantile(hh, .75, na.rm = TRUE), best_h = ref[g_best],
    pct_grid_min = 100 * mean(g == 1, na.rm = TRUE), pct_grid_max = 100 * mean(g == 100, na.rm = TRUE),
    bias = mean(e), bias_mcse = sd(e) / sqrt(m), ESE = sd(e),
    MSE = mean(e^2), MSE_mcse = sd(e^2) / sqrt(m),
    best_MSE = mean(eb^2), best_MSE_mcse = sd(eb^2) / sqrt(length(eb)),
    MSE_ratio = mean(e^2) / mean(eb^2),
    CP = 100 * cp, CP_mcse = 100 * sqrt(cp * (1 - cp) / m),
    p_vs95 = binom.test(kc, m, 0.95)$p.value, row.names = NULL)
}
tab <- do.call(rbind, tab)
tab <- tab[order(tab$link, tab$n), ]
tab$q_vs95 <- p.adjust(tab$p_vs95, "BH")
tab$sig <- ifelse(tab$q_vs95 < 0.001, "***", ifelse(tab$q_vs95 < 0.01, "**", ifelse(tab$q_vs95 < 0.05, "*", "")))
saveRDS(tab, file.path(cvdir, "agg_dd_cv.rds"))
fmt <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], signif, 3); df }
print(fmt(tab), row.names = FALSE)
# Trend of the MSE ratio with n, per link: slope of log ratio on log n.
for (L in unique(tab$link)) {
  t <- tab[tab$link == L, ]
  f <- lm(log(MSE_ratio) ~ log(n), t)
  cat(sprintf("%s: log(MSE ratio) ~ log(n) slope %.3f (SE %.3f)\n", L, coef(f)[2], sqrt(vcov(f)[2, 2])))
}

## Aggregate the loss-function comparison. Usage: Rscript agg_acv4.R <dir>
dir <- commandArgs(TRUE)[1]
fs <- Sys.glob(file.path(dir, "*_n*_c*.rds"))
objs <- lapply(fs, readRDS)
grid <- objs[[1]]$grid; G <- length(grid); z <- qnorm(0.975)
cells <- split(objs, sapply(objs, function(o) paste(o$link, o$n)))
star <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))
summ <- list(); tests <- list(); fails <- list()
for (cell in cells) {
  link <- cell[[1]]$link; n <- cell[[1]]$n; b0 <- cell[[1]]$beta0[2]
  R <- unlist(lapply(cell, `[[`, "res"), recursive = FALSE)
  R <- R[!duplicated(sapply(R, `[[`, "rep"))]
  err <- sapply(R, function(r) !is.null(r$error)); Rk <- R[!err]
  B <- t(sapply(Rk, `[[`, "b2")); S <- t(sapply(Rk, `[[`, "se2")); nr <- nrow(B)
  P <- t(sapply(Rk, `[[`, "picks"))
  mse_h <- colMeans((B - b0)^2, na.rm = TRUE); g_best <- which.min(mse_h)
  P <- cbind(P, best_fixed_h = g_best)
  nf <- t(sapply(Rk, `[[`, "nfail_h"))
  fails[[length(fails) + 1]] <- data.frame(link = link, n = n, N = nr, rep_errors = sum(err),
    fullfit_fail_pct = 100 * mean(!is.finite(B)), cvfit_fail_pct = 100 * mean(nf) / 5,
    best_h = grid[g_best], best_mse = mse_h[g_best])
  sq <- list()
  for (s in colnames(P)) {
    g <- P[, s]; est <- B[cbind(seq_len(nr), g)]; se <- S[cbind(seq_len(nr), g)]
    ok <- !is.na(g) & is.finite(est) & is.finite(se)
    e <- est[ok] - b0; m <- sum(ok); k <- sum(abs(e) <= z * se[ok])
    sq[[s]] <- ifelse(ok, (est - b0)^2, NA)
    summ[[length(summ) + 1]] <- data.frame(link = link, n = n, loss = s, N = nr,
      no_choice_pct = 100 * mean(!ok),
      h_q25 = unname(quantile(grid[g[ok]], .25)), h_med = median(grid[g[ok]]),
      h_q75 = unname(quantile(grid[g[ok]], .75)),
      pct_grid_min = 100 * mean(g[ok] == 1), pct_grid_max = 100 * mean(g[ok] == G),
      bias = mean(e), bias_mcse = sd(e) / sqrt(m), ESE = sd(e),
      MSE = mean(e^2), MSE_mcse = sd(e^2) / sqrt(m), MSE_ratio = mean(e^2) / mse_h[g_best],
      CP = 100 * k / m, CP_lo = 100 * binom.test(k, m)$conf.int[1],
      CP_hi = 100 * binom.test(k, m)$conf.int[2], p_CP95 = binom.test(k, m, .95)$p.value)
  }
  for (s in setdiff(colnames(P), c("matched", "best_fixed_h"))) {
    dd <- sq[[s]] - sq[["matched"]]; dd <- dd[is.finite(dd)]
    tests[[length(tests) + 1]] <- data.frame(link = link, n = n, loss = s,
      diff_vs_matched = mean(dd), se = sd(dd) / sqrt(length(dd)), N = length(dd),
      p = t.test(dd)$p.value)
  }
}
summ <- do.call(rbind, summ); tests <- do.call(rbind, tests); fails <- do.call(rbind, fails)
summ$q_CP95 <- p.adjust(summ$p_CP95, "BH"); summ$sig_CP <- star(summ$q_CP95)
tests$q <- p.adjust(tests$p, "BH"); tests$sig <- star(tests$q)
summ <- summ[order(summ$link, summ$n, summ$loss == "best_fixed_h", summ$MSE), ]
fmt <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], signif, 3); df }
print(fmt(fails), row.names = FALSE); cat("\n")
print(fmt(summ[, c("link","n","loss","N","no_choice_pct","MSE","MSE_mcse","MSE_ratio","bias","bias_mcse","ESE","CP","CP_lo","CP_hi","q_CP95","sig_CP")]), row.names = FALSE)
cat("\n")
print(fmt(summ[, c("link","n","loss","h_q25","h_med","h_q75","pct_grid_min","pct_grid_max")]), row.names = FALSE)
cat("\nPaired squared-error difference vs matched (negative = better than matched)\n")
print(fmt(tests), row.names = FALSE)
saveRDS(list(summary = summ, tests = tests, fails = fails), file.path(dir, "summary4.rds"))

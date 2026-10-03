## Per-h diagnostics from the sim4 results: bias/ESE, mean SE/ESE and coverage
## at every grid h, and at the h each loss selected.
## Usage: Rscript perh_diag.R <dir>
dir <- commandArgs(TRUE)[1]
objs <- lapply(Sys.glob(file.path(dir, "*_n*_c*.rds")), readRDS)
grid <- objs[[1]]$grid; G <- length(grid); z <- qnorm(0.975)
cells <- split(objs, sapply(objs, function(o) paste(o$link, o$n)))
perh <- list(); atsel <- list()
for (cell in cells) {
  link <- cell[[1]]$link; n <- cell[[1]]$n; b0 <- cell[[1]]$beta0[2]
  R <- unlist(lapply(cell, `[[`, "res"), recursive = FALSE)
  R <- R[!sapply(R, function(r) !is.null(r$error))]
  B <- t(sapply(R, `[[`, "b2")); S <- t(sapply(R, `[[`, "se2"))
  for (g in seq_len(G)) {
    ok <- is.finite(B[, g]) & is.finite(S[, g]); e <- B[ok, g] - b0
    perh[[length(perh) + 1]] <- data.frame(link = link, n = n, h = grid[g], N = sum(ok),
      bias = mean(e), ESE = sd(e), bias_over_ESE = mean(e) / sd(e),
      SE_over_ESE = mean(S[ok, g]) / sd(e), MSE = mean(e^2),
      CP = 100 * mean(abs(e) <= z * S[ok, g]),
      CP_mcse = 100 * sqrt(mean(abs(e) <= z * S[ok, g]) * (1 - mean(abs(e) <= z * S[ok, g])) / sum(ok)))
  }
  P <- t(sapply(R, `[[`, "picks"))
  mse <- colMeans((B - b0)^2, na.rm = TRUE); P <- cbind(P, best_fixed_h = which.min(mse))
  for (s in colnames(P)) {
    g <- P[, s]; est <- B[cbind(seq_along(g), g)]; se <- S[cbind(seq_along(g), g)]
    ok <- is.finite(est) & is.finite(se); e <- est[ok] - b0
    # Mean of the per-h bias/ESE at the chosen h: how far the chosen h sits
    # into the bias-dominated region.
    bE <- sapply(seq_len(G), function(k) { v <- B[, k] - b0; mean(v, na.rm = TRUE) / sd(v, na.rm = TRUE) })
    seE <- sapply(seq_len(G), function(k) mean(S[, k], na.rm = TRUE) / sd(B[, k], na.rm = TRUE))
    atsel[[length(atsel) + 1]] <- data.frame(link = link, n = n, loss = s, N = sum(ok),
      h_med = median(grid[g[ok]]),
      mean_bias_over_ESE_at_h = mean(bE[g[ok]]), mean_SE_over_ESE_at_h = mean(seE[g[ok]]),
      CP = 100 * mean(abs(e) <= z * se[ok]))
  }
}
perh <- do.call(rbind, perh); atsel <- do.call(rbind, atsel)
fmt <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], signif, 3); df }
saveRDS(list(perh = perh, atsel = atsel), file.path(dir, "perh_diag.rds"))
write.csv(perh, file.path(dir, "perh_diag.csv"), row.names = FALSE)
show <- perh[perh$h %in% grid[seq(1, G, by = 4)] | perh$h %in% grid[c(G)], ]
print(fmt(show), row.names = FALSE)
cat("\nAt the selected h\n")
print(fmt(atsel), row.names = FALSE)

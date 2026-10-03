## Undersmoothed confidence intervals for kee_async(), offline from the sim4
## results: point estimate at h_hat (fixed_lower CV), interval at
## h_ci = c * h_hat, mapped to the grid value nearest in log h. A tie between
## two grid values goes to the smaller one; h_ci below the smallest grid value
## is set to it (counted as a floor hit).
## Usage: Rscript undersmooth.R <sim4 results dir> [out.rds]
a <- commandArgs(TRUE); dir <- a[1]
out <- if (length(a) > 1) a[2] else file.path(dir, "undersmooth.rds")
objs <- lapply(Sys.glob(file.path(dir, "*_n*_c*.rds")), readRDS)
grid <- objs[[1]]$grid; G <- length(grid); lg <- log(grid); z <- qnorm(0.975)
to_grid <- function(h) {
  below <- h < grid[1]
  d <- abs(outer(log(pmax(h, grid[1])), lg, "-"))
  idx <- apply(d, 1, function(v) which(v <= min(v) + 1e-12)[1])  # ties -> smaller h
  list(idx = idx, floor = below)
}
rules <- c(list(c1 = 1, c075 = 0.75, c05 = 0.5, c035 = 0.35, c025 = 0.25),
           list(rate01 = -0.1, rate02 = -0.2, rate03 = -0.3))
cells <- split(objs, sapply(objs, function(o) paste(o$link, o$n)))
tab <- list()
for (cell in cells) {
  link <- cell[[1]]$link; n <- cell[[1]]$n; b0 <- cell[[1]]$beta0[2]
  R <- unlist(lapply(cell, `[[`, "res"), recursive = FALSE)
  R <- R[!sapply(R, function(r) !is.null(r$error))]
  B <- t(sapply(R, `[[`, "b2")); S <- t(sapply(R, `[[`, "se2")); nr <- nrow(B)
  ghat <- sapply(R, function(r) r$picks[["fixed_lower"]])
  mse <- colMeans((B - b0)^2, na.rm = TRUE); g_best <- which.min(mse)
  bE <- (colMeans(B, na.rm = TRUE) - b0) / apply(B, 2, sd, na.rm = TRUE)
  sE <- colMeans(S, na.rm = TRUE) / apply(B, 2, sd, na.rm = TRUE)
  spec <- c(lapply(names(rules), function(k) {
    v <- rules[[k]]; cc <- if (v > 0) v else n^v
    m <- to_grid(cc * grid[ghat]); list(rule = k, c = cc, idx = m$idx, floor = m$floor)
  }), list(list(rule = "best_fixed_h", c = NA, idx = rep(g_best, nr), floor = rep(FALSE, nr))))
  for (sp in spec) {
    g <- sp$idx; est <- B[cbind(seq_len(nr), g)]; se <- S[cbind(seq_len(nr), g)]
    ok <- is.finite(est) & is.finite(se); e <- est[ok] - b0
    k <- sum(abs(e) <= z * se[ok]); m <- sum(ok); cp <- k / m
    tab[[length(tab) + 1]] <- data.frame(link = link, n = n, rule = sp$rule, c = sp$c, N = m,
      h_ci_med = median(grid[g[ok]]), floor_pct = 100 * mean(sp$floor[ok] | g[ok] == 1),
      CP = 100 * cp, CP_mcse = 100 * sqrt(cp * (1 - cp) / m),
      p_vs95 = binom.test(k, m, 0.95)$p.value,
      width = mean(2 * z * se[ok]), width_mcse = sd(2 * z * se[ok]) / sqrt(m),
      bias_over_ESE = mean(bE[g[ok]]), SE_over_ESE = mean(sE[g[ok]]))
  }
}
tab <- do.call(rbind, tab)
tab$q_vs95 <- p.adjust(tab$p_vs95, "BH")
tab$sig <- ifelse(tab$q_vs95 < 0.001, "***", ifelse(tab$q_vs95 < 0.01, "**", ifelse(tab$q_vs95 < 0.05, "*", "")))
saveRDS(tab, out)
fmt <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], signif, 3); df }
for (key in unique(paste(tab$link, tab$n))) {
  t <- tab[paste(tab$link, tab$n) == key, ]
  t <- t[order(t$rule == "best_fixed_h", abs(t$CP - 95), t$width), ]
  print(fmt(t), row.names = FALSE); cat("\n")
}
# Ranking across cells: mean |CP - 95|, worst |CP - 95|, mean width.
rk <- aggregate(cbind(absdev = abs(tab$CP - 95), width = tab$width) ~ rule, tab, mean)
rk$worst_absdev <- aggregate(abs(tab$CP - 95) ~ rule, tab, max)[, 2]
print(fmt(rk[order(rk$absdev, rk$width), ]), row.names = FALSE)

## Aggregate the data-driven-grid undersmoothing study (sim_dd.R output).
## h_ci = c * h_hat is snapped to the data set's own grid (nearest in log h,
## ties to the smaller h); below the data set's h_min it is set to h_min and
## counted as a floor hit. A failed fit at h_ci is a failure: coverage is
## reported among non-failures and, conservatively, with failures counted as
## misses. bias/ESE and SE/ESE at h_ci are read off the common reference grid
## (130 points, 1e-8 to 1) at the reference point nearest h_ci. The
## yardstick is the reference-grid h with the smallest simulated MSE among
## points that fit in at least 95% of replicates.
## Usage: Rscript undersmooth_dd.R <dir>
dir <- commandArgs(TRUE)[1]
objs <- lapply(Sys.glob(file.path(dir, "*_n*_r*.rds")), readRDS)
ref <- objs[[1]]$ref; z <- qnorm(0.975); lref <- log(ref)
rules <- list(c1 = 1, c075 = 0.75, c05 = 0.5, c035 = 0.35, c025 = 0.25,
              rate01 = -0.1, rate02 = -0.2, rate03 = -0.3)
cells <- split(objs, sapply(objs, function(o) paste(o$link, o$n)))
tab <- list(); dist <- list(); perh <- list()
for (cell in cells) {
  link <- cell[[1]]$link; n <- cell[[1]]$n; b0 <- cell[[1]]$beta0[2]
  R <- unlist(lapply(cell, `[[`, "res"), recursive = FALSE)
  R <- R[!duplicated(sapply(R, `[[`, "rep"))]
  err <- sapply(R, function(r) !is.null(r$error)); R <- R[!err]; nr <- length(R)
  B <- t(sapply(R, `[[`, "b_ref")); S <- t(sapply(R, `[[`, "se_ref"))
  okh <- colMeans(is.finite(B) & is.finite(S))
  e <- B - b0
  ese <- apply(B, 2, sd, na.rm = TRUE); bias <- colMeans(e, na.rm = TRUE)
  mse <- colMeans(e^2, na.rm = TRUE); mse[okh < 0.95] <- NA
  bE <- bias / ese; sE <- colMeans(S, na.rm = TRUE) / ese
  cpr <- colMeans(abs(e) <= z * S, na.rm = TRUE)
  perh[[length(perh) + 1]] <- data.frame(link = link, n = n, h = ref, ok_pct = 100 * okh,
    bias_over_ESE = bE, SE_over_ESE = sE, MSE = mse, CP = 100 * cpr)
  g_best <- which.min(mse)
  hh <- sapply(R, function(r) if (is.na(r$pick)) NA else r$grid[r$pick])
  dist[[length(dist) + 1]] <- data.frame(link = link, n = n, N = nr, rep_errors = sum(err),
    h_min_med = median(sapply(R, `[[`, "hmin")), h_min_q99 = quantile(sapply(R, `[[`, "hmin"), .99),
    h_max_med = median(sapply(R, `[[`, "hmax")),
    h_hat_q25 = quantile(hh, .25, na.rm = TRUE), h_hat_med = median(hh, na.rm = TRUE),
    h_hat_q75 = quantile(hh, .75, na.rm = TRUE), no_pick = sum(is.na(hh)),
    own_fail_pct = 100 * mean(sapply(R, function(r) mean(is.na(r$b_own)))),
    cv_excluded_pct = 100 * mean(sapply(R, function(r) mean(r$cv_excluded))),
    best_h = ref[g_best], best_mse = mse[g_best], row.names = NULL)
  specs <- c(lapply(names(rules), function(k) { v <- rules[[k]]; list(rule = k, c = if (v > 0) v else n^v) }),
             list(list(rule = "best_fixed_h", c = NA)))
  for (sp in specs) {
    est <- se <- hci <- rep(NA_real_, nr); fl <- rep(FALSE, nr)
    for (i in seq_len(nr)) {
      r <- R[[i]]
      if (sp$rule == "best_fixed_h") { hci[i] <- ref[g_best]; est[i] <- r$b_ref[g_best]; se[i] <- r$se_ref[g_best]; next }
      if (is.na(r$pick)) next
      target <- sp$c * r$grid[r$pick]
      fl[i] <- target < r$grid[1]
      dd <- abs(log(max(target, r$grid[1])) - log(r$grid))
      k <- which(dd <= min(dd) + 1e-12)[1]
      hci[i] <- r$grid[k]; est[i] <- r$b_own[k]; se[i] <- r$se_own[k]
    }
    ok <- is.finite(est) & is.finite(se); m <- sum(ok); fail <- nr - m
    cov <- abs(est[ok] - b0) <= z * se[ok]; k <- sum(cov); cp <- k / m
    kr <- sapply(hci[ok], function(h) which.min(abs(lref - log(h))))
    tab[[length(tab) + 1]] <- data.frame(link = link, n = n, rule = sp$rule, c = sp$c, N = nr,
      h_ci_q25 = quantile(hci, .25, na.rm = TRUE), h_ci_med = median(hci, na.rm = TRUE),
      h_ci_q75 = quantile(hci, .75, na.rm = TRUE),
      floor_pct = 100 * mean(fl), fail_pct = 100 * fail / nr,
      CP = 100 * cp, CP_mcse = 100 * sqrt(cp * (1 - cp) / m),
      CP_fail_as_miss = 100 * k / nr, p_vs95 = binom.test(k, m, 0.95)$p.value,
      width = mean(2 * z * se[ok]), width_mcse = sd(2 * z * se[ok]) / sqrt(m),
      bias_over_ESE = mean(bE[kr], na.rm = TRUE), SE_over_ESE = mean(sE[kr], na.rm = TRUE),
      row.names = NULL)
  }
}
tab <- do.call(rbind, tab); dist <- do.call(rbind, dist); perh <- do.call(rbind, perh)
tab$q_vs95 <- p.adjust(tab$p_vs95, "BH")
tab$sig <- ifelse(tab$q_vs95 < 0.001, "***", ifelse(tab$q_vs95 < 0.01, "**", ifelse(tab$q_vs95 < 0.05, "*", "")))
saveRDS(list(table = tab, dist = dist, perh = perh), file.path(dir, "undersmooth_dd.rds"))
fmt <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], signif, 3); df }
print(fmt(dist), row.names = FALSE); cat("\n")
cpw <- reshape(tab[, c("link", "n", "rule", "CP")], idvar = c("link", "rule"), timevar = "n", direction = "wide")
print(fmt(cpw[order(cpw$link, cpw$rule), ]), row.names = FALSE); cat("\n")
print(fmt(tab[order(tab$link, tab$n, tab$rule == "best_fixed_h", abs(tab$CP - 95)),
              c("link","n","rule","c","CP","CP_mcse","CP_fail_as_miss","q_vs95","sig","width","bias_over_ESE","SE_over_ESE","floor_pct","fail_pct","h_ci_med")]), row.names = FALSE)
rk <- aggregate(cbind(absdev = abs(tab$CP - 95), width = tab$width) ~ rule, tab, mean)
rk$worst <- aggregate(abs(tab$CP - 95) ~ rule, tab, max)[, 2]
cat("\n"); print(fmt(rk[order(rk$absdev, rk$width), ]), row.names = FALSE)

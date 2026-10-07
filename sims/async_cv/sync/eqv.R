## Package (worktree) vs cov6.R on identical data and folds.
source("cov6.R")
WT <- "/mnt/fanxiangssd/files/coding/skmle/.claude/worktrees/agent-a7d2d695509039df3"
suppressMessages(pkgload::load_all(WT, quiet = TRUE, export_all = TRUE, helpers = FALSE))
ns <- asNamespace("skmle")
cmp <- function(link, n, rep_id) {
  beta0 <- BETA[[link]]; set.seed(2e7 + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  dat <- sim_link(n, beta0, link); e <- prep(dat)
  s6 <- select_g(e, link, beta0); f6 <- fit_g(e$P, s6[["cvc"]], link, TRUE, e$gam, e$m, e$s2, beta0)
  d <- ns$.async_parse(y ~ x, dat$y, dat$x, "id", "time", TRUE); lk <- ns$.async_link_code(link)
  vg <- ns$.async_variogram(d, 2L); P <- ns$.async_all_pairs(d, 2L, FALSE)
  folds <- matrix(((seq_len(d$n) - 1L) %% 5L) + 1L)
  g <- ns$.async_lag_grid(d, 100L, FALSE); lo <- ns$.async_grid_floor(P, d$n, folds)
  g <- exp(seq(log(lo), log(max(g)), length.out = 100))
  st <- as.numeric(ns$.async_fit_ti(d, max(g), FALSE, lk, 50L, 1e-8)$beta)
  cv <- ns$.async_cv_sync(d, P, vg, 2L, g, folds, FALSE, lk, st, 50L, 1e-10, TRUE)
  hp <- g[cv$index]
  fit <- kee_async(dat$y, dat$x, y ~ x, id = id, time = time, h = hp, link = link, correct = TRUE, tol = 1e-10)
  data.frame(link, n, rep = rep_id, floor_diff = abs(lo - s6[["lo"]]) / lo, h6 = s6[["cvc"]], hpkg = hp,
    b_diff = max(abs(coef(fit) - f6$beta)), se_diff = max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(f6$V)))))
}
J <- expand.grid(link = c("identity", "log", "logistic"), n = c(400, 2000), rep = 1:4, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) tryCatch(cmp(J$link[j], J$n[j], J$rep[j]), error = function(e) conditionMessage(e)), mc.cores = 16)
print(do.call(rbind, Filter(is.data.frame, r)), digits = 4); print(Filter(is.character, r))

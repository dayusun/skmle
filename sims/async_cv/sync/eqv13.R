## Package kee_async_cv(method = "variance", the default) vs the simulation reference select13.
source("cov13.R")
WT <- "/mnt/fanxiangssd/files/coding/skmle/.claude/worktrees/agent-a7d2d695509039df3"
cmp <- function(gen, link, n, rep_id) {
  suppressMessages(pkgload::load_all(WT, quiet = TRUE, export_all = TRUE, helpers = FALSE))
  beta0 <- BETA[[link]]; set.seed(5678 + rep_id + 1e3 * n)
  dat <- GENS[[gen]](n, beta0, link); e <- prep8(dat); r <- ratio_hat(e$p)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(ncol(e$X)))$beta, error = function(err) numeric(ncol(e$X)))
  ref <- select13(e, link, st, r)
  cv <- suppressWarnings(kee_async_cv(dat$y, dat$x, y ~ x1 + x2 + z, id = id, time = time, link = link, tol = 1e-10, quiet = TRUE))
  ok <- is.finite(ref$crit) & is.finite(cv$cv_results$trSV)
  data.frame(gen, link, n, rep = rep_id, grid_diff = max(abs(cv$h_grid - ref$grid) / ref$grid),
    crit_rel = max(abs(cv$cv_results$trSV[ok] - ref$crit[ok]) / ref$crit[ok]), n_fin_ref = sum(is.finite(ref$crit)), n_fin_pkg = sum(is.finite(cv$cv_results$trSV)),
    h_ref = ref$h, h_pkg = cv$h_cv, b_diff = max(abs(coef(cv$fit) - ref$EST[ref$g, ])), se_diff = max(abs(sqrt(diag(vcov(cv$fit))) - ref$SE[ref$g, ])))
}
J <- expand.grid(gen = c("gauss", "skew_early"), link = c("identity", "log", "logistic"), n = c(400, 2000), rep = 1:2, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) tryCatch(cmp(J$gen[j], J$link[j], J$n[j], J$rep[j]), error = function(e) conditionMessage(e)), mc.cores = 16)
options(width = 200); print(do.call(rbind, Filter(is.data.frame, r)), digits = 3, row.names = FALSE); print(Filter(is.character, r))

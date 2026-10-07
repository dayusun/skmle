## Package (worktree) vs scratch reference for the pair-moment method with sync CV.
source("cov11.R")
sat_ref <- function(y, link) if (link == "log") ifelse(y > 0, y - y * log(y), 0) else 0 * y
select_ref <- function(e, link, start, r) {
  P <- e$P; Q <- e$Q; pl <- loss_of(link); folds <- ((e$ids - 1L) %% K) + 1L
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), floor9(P, Q, e$p$n, folds))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG)); IC <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) { te <- which(folds == f)
    H <- P[P[, "sid"] %in% te, , drop = FALSE]; Ptr <- P[!P[, "sid"] %in% te, , drop = FALSE]; Qtr <- Q[!Q[, "sid"] %in% te, , drop = FALSE]
    HX <- e$X[H[, "k"], , drop = FALSE]; QR <- qr(cbind(1, e$gb[e$block_of(abs(H[, "u"])), , drop = FALSE])); s0 <- mean(sat_ref(H[, "y"], link)); st <- start
    for (g in rev(seq_len(NG))) { ft <- tryCatch(fit_pmw(e, Ptr, Qtr, grid[g], link, st, r), error = function(err) NULL)
      if (is.null(ft) || !ft$conv) next; st <- ft$beta
      ic <- qr.coef(QR, pl(H[, "y"], as.vector(HX %*% ft$beta)))[1]; if (is.finite(ic) && ic >= s0) IC[g, f] <- ic } }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; list(h = grid[which.min(cv)], grid = grid)
}
WT <- "/mnt/fanxiangssd/files/coding/skmle/.claude/worktrees/agent-a7d2d695509039df3"
cmp <- function(gen, link, n, rep_id) {
  suppressMessages(pkgload::load_all(WT, quiet = TRUE, export_all = TRUE, helpers = FALSE)); ns <- asNamespace("skmle")
  beta0 <- BETA[[link]]; set.seed(1234 + rep_id + 1e3 * n)
  dat <- GENS[[gen]](n, beta0, link); e <- prep8(dat); r <- ratio_hat(e$p); lk <- LCODE[[link]]
  st <- as.numeric(fit_(e$p, max(abs(e$P[, "u"])), FALSE, lk, 50L, 1e-8)$beta)
  ref <- select_ref(e, link, st, r); fr <- fit_pmw(e, e$P, e$Q, ref$h, link, st, r)
  d <- ns$.async_parse(y ~ x1 + x2 + z, dat$y, dat$x, "id", "time", TRUE)
  P <- ns$.async_ts_pairs(d, FALSE); Q <- ns$.async_cov_pairs(d, FALSE); rp <- ns$.async_intensity_ratio(d)
  tv <- ns$.async_time_varying(d); gdes <- ns$.async_gamma_design(d, tv); folds <- matrix(((seq_len(d$n) - 1L) %% 5L) + 1L)
  lo <- max(ns$.async_grid_floor(P, d$n, folds), ns$.async_grid_floor(cbind(u = Q[, "v"], sid = Q[, "sid"]), d$n, folds))
  g <- ns$.async_lag_grid(d, 100L, FALSE); g <- exp(seq(log(lo), log(max(g)), length.out = 100))
  fitter <- function(te, h, s) ns$.async_fit_pairs(P[!P[, "sid"] %in% te, , drop = FALSE], Q[!Q[, "sid"] %in% te, , drop = FALSE], h, d, lk, rp, s, 50L, 1e-10, FALSE)
  cv <- ns$.async_cv_sync(P, function(H) d$X[H[, "k"], , drop = FALSE], gdes, fitter, g, folds, lk, st, TRUE)
  hp <- g[cv$index]
  fit <- kee_async(dat$y, dat$x, y ~ x1 + x2 + z, id = id, time = time, h = hp, link = link, correct = "pairs", tol = 1e-10)
  data.frame(gen, link, n, rep = rep_id, grid_diff = max(abs(g - ref$grid) / g), r_diff = max(abs(rp - r) / r),
    h_ref = ref$h, h_pkg = hp, b_diff = max(abs(coef(fit) - fr$beta)), se_diff = max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(fr$V)))))
}
J <- expand.grid(gen = c("gauss", "skew_early"), link = c("identity", "log", "logistic"), n = c(400, 2000), rep = 1:2, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) tryCatch(cmp(J$gen[j], J$link[j], J$n[j], J$rep[j]), error = function(e) conditionMessage(e)), mc.cores = 16)
options(width = 200); print(do.call(rbind, Filter(is.data.frame, r)), digits = 4, row.names = FALSE); print(Filter(is.character, r))

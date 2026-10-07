## For wild log reps: is the selected (floor) candidate's fold intercept below the saturated loss?
source("cov9.R")
L <- unlist(lapply(Sys.glob("results9/*_log_*.rds"), readRDS), recursive = FALSE); x <- do.call(rbind, Filter(is.data.frame, L))
x <- x[x$h_type == "pm", ]; x$err <- abs(x$est - x$truth)
w <- aggregate(err ~ gen + n + rep + h + h_lo, x[x$coef != "(Intercept)", ], max)
bad <- head(w[w$err > 0.5 & w$n <= 2000, ], 24); good <- head(w[w$err < 0.1 & w$n == 2000, ], 8)
J <- rbind(cbind(bad, kind = "wild"), cbind(good, kind = "ok"))
sat_loss <- function(y) ifelse(y > 0, y - y * log(y), 0)   # min_eta exp(eta) - y eta
r <- parallel::mclapply(seq_len(nrow(J)), function(j) {
  gen <- J$gen[j]; n <- J$n[j]; link <- "log"; rep_id <- J$rep[j]; beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); rr <- ratio_hat(e$p)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(length(nm)))$beta, error = function(err) numeric(length(nm)))
  P <- e$P; Q <- e$Q; pl <- loss_of(link); folds <- ((e$ids - 1L) %% K) + 1L
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), floor9(P, Q, e$p$n, folds))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG)); IC <- SAT <- matrix(NA, NG, K)
  for (f in seq_len(K)) { te <- which(folds == f)
    H <- P[P[, "sid"] %in% te, , drop = FALSE]; Ptr <- P[!P[, "sid"] %in% te, , drop = FALSE]; Qtr <- Q[!Q[, "sid"] %in% te, , drop = FALSE]
    HX <- e$X[H[, "k"], , drop = FALSE]; QR <- qr(cbind(1, e$gb[e$block_of(abs(H[, "u"])), , drop = FALSE])); s0 <- mean(sat_loss(H[, "y"])); stt <- st
    for (g in rev(seq_len(NG))) { ft <- tryCatch(fit_pmw(e, Ptr, Qtr, grid[g], link, stt, rr), error = function(err) NULL)
      if (is.null(ft) || !ft$conv) next; stt <- ft$beta
      IC[g, f] <- qr.coef(QR, pl(H[, "y"], as.vector(HX %*% ft$beta)))[1]; SAT[g, f] <- s0 } }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; g0 <- which.min(cv)
  below <- rowSums(IC < SAT, na.rm = TRUE) > 0
  cv2 <- cv; cv2[below] <- Inf; g1 <- which.min(cv2)
  b1 <- fit_pmw(e, P, Q, grid[g1], link, st, rr)$beta
  data.frame(J[j, c("gen", "n", "rep", "kind")], h_old = grid[g0], old_err = J$err[j], below_at_old = below[g0],
             n_below = sum(below), h_new = grid[g1], new_err = max(abs(b1 - beta0)[-1])) }, mc.cores = 16)
d <- do.call(rbind, Filter(is.data.frame, r)); saveRDS(d, "sat9.rds"); options(width = 200); print(d, digits = 3, row.names = FALSE)

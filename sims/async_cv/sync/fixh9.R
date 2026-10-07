## Same seeds as job 8478458: pair-moment at fixed h vs the CV-selected h, n = 8000.
source("cov9.R")
L <- unlist(lapply(Sys.glob("results9/skew*_identity_8000_*.rds"), readRDS), recursive = FALSE)
job <- do.call(rbind, Filter(is.data.frame, L)); job <- job[job$h_type == "pm", ]
J <- expand.grid(gen = c("skew", "skew_early"), rep = 1:96, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) {
  gen <- J$gen[j]; link <- "identity"; n <- 8000; rep_id <- J$rep[j]; beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); r4 <- ratio_hat(e$p); r12 <- ratio_hat(e$p, 12)
  rt <- if (gen == "skew_early") 1 / dbeta(e$p$xt, 1, 3) else NULL
  out <- do.call(rbind, lapply(c(0.05, 0.1, 0.3, 0.6), function(h) do.call(rbind, lapply(c("r4", "r12", "rtrue"), function(m) {
    if (m == "rtrue" && is.null(rt)) return(NULL)
    f <- fit_pmw(e, e$P, e$Q, h, link, beta0, switch(m, r4 = r4, r12 = r12, rtrue = rt))
    data.frame(gen, rep = rep_id, h = h, m, coef = nm, err = f$beta - beta0, se = sqrt(diag(f$V))) }))))
  jb <- job[job$gen == gen & job$rep == rep_id, ]
  rbind(out, data.frame(gen, rep = rep_id, h = NA, m = "cv_job", coef = jb$coef, err = jb$est - jb$truth, se = jb$se)) }, mc.cores = 16)
x <- do.call(rbind, Filter(is.data.frame, r)); x <- x[x$coef != "(Intercept)", ]
f <- function(d) data.frame(R = nrow(d), bias = round(mean(d$err), 4), z = round(mean(d$err) / (sd(d$err) / sqrt(nrow(d))), 1),
  ESE = round(sd(d$err), 4), SE_ESE = round(median(d$se) / sd(d$err), 2), CP = round(100 * mean(abs(d$err) <= 1.96 * d$se), 1))
a <- do.call(rbind, lapply(split(x, paste(x$gen, x$coef, x$m, x$h)), function(d) cbind(d[1, c("gen", "coef", "m", "h")], f(d))))
options(width = 150); print(a, row.names = FALSE); saveRDS(x, "fixh9.rds")

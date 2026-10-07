source("cov10.R")
J <- expand.grid(gen = c("skew", "skew_early"), rep = 1:96, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) {
  gen <- J$gen[j]; link <- "identity"; n <- 8000; rep_id <- J$rep[j]; beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); r <- ratio_hat(e$p)
  folds <- ((e$ids - 1L) %% K) + 1L; folds <- folds[order(((e$ids * 7919L) %% e$p$n))]
  a <- fit_pm_cf(e, c(0.05, 0.1, 0.3, 0.6, 1), folds, link, beta0, r)
  b <- fit_pm_cf(e, rep(0.3, 5), folds, link, beta0, r); c0 <- fit_pmw(e, e$P, e$Q, 0.3, link, beta0, r)
  rbind(data.frame(gen, m = "cf_fixed_mixed", coef = nm, err = a$beta - beta0, se = sqrt(diag(a$V))),
        data.frame(gen, m = "cf_fixed_0.3", coef = nm, err = b$beta - beta0, se = sqrt(diag(b$V))),
        data.frame(gen, m = "pooled_0.3", coef = nm, err = c0$beta - beta0, se = sqrt(diag(c0$V)))) }, mc.cores = 16)
x <- do.call(rbind, Filter(is.data.frame, r)); x <- x[x$coef != "(Intercept)", ]
f <- function(d) data.frame(R = nrow(d), bias = round(mean(d$err), 4), ESE = round(sd(d$err), 4), SE_ESE = round(median(d$se) / sd(d$err), 2),
  CP = round(100 * mean(abs(d$err) <= 1.96 * d$se), 1))
print(do.call(rbind, lapply(split(x, paste(x$gen, x$m, x$coef)), function(d) cbind(d[1, c("gen", "m", "coef")], f(d)))), row.names = FALSE)

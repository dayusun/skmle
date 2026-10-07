source("cov11.R")
J <- expand.grid(gen = c("skew_early", "skew"), link = c("identity", "logistic"), rep = 1:48, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) {
  gen <- J$gen[j]; link <- J$link[j]; n <- 8000; rep_id <- J$rep[j]; beta0 <- BETA[[link]]
  set.seed(4e7 + 1e6 * match(gen, names(GENS)) + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep8(GENS[[gen]](n, beta0, link)); nm <- colnames(e$X); rf <- ratio_fit(e$p)
  do.call(rbind, lapply(c(0.1, 0.3), function(h) { f <- fit_pms(e, e$P, e$Q, h, link, beta0, rf)
    data.frame(gen, link, h, coef = nm, err = f$beta - beta0, se = sqrt(diag(f$V)), se0 = sqrt(diag(f$V0))) })) }, mc.cores = 16)
x <- do.call(rbind, Filter(is.data.frame, r)); x <- x[x$coef != "(Intercept)", ]
f <- function(d) data.frame(R = nrow(d), ESE = round(sd(d$err), 4), SE0_ESE = round(median(d$se0) / sd(d$err), 2), SE_ESE = round(median(d$se) / sd(d$err), 2),
  CP0 = round(100 * mean(abs(d$err) <= 1.96 * d$se0), 1), CP = round(100 * mean(abs(d$err) <= 1.96 * d$se), 1))
print(do.call(rbind, lapply(split(x, paste(x$gen, x$link, x$h, x$coef)), function(d) cbind(d[1, c("gen", "link", "h", "coef")], f(d)))), row.names = FALSE)

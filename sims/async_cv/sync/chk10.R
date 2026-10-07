source("cov10.R")
J <- expand.grid(gen = c("skew", "skew_early"), rep = 1:96, stringsAsFactors = FALSE)
r <- parallel::mclapply(seq_len(nrow(J)), function(j) tryCatch(one_rep10(8000, "identity", J$rep[j], J$gen[j]),
  error = function(e) conditionMessage(e)), mc.cores = 16)
cat("errors", sum(vapply(r, is.character, TRUE)), "\n")
x <- do.call(rbind, Filter(is.data.frame, r)); x <- x[x$coef != "(Intercept)", ]; x$err <- x$est - x$truth
f <- function(d) data.frame(R = nrow(d), bias = round(mean(d$err), 4), z = round(mean(d$err) / (sd(d$err) / sqrt(nrow(d))), 1),
  ESE = round(sd(d$err), 4), SE_ESE = round(median(d$se) / sd(d$err), 2), CP = round(100 * mean(abs(d$err) <= 1.96 * d$se), 1),
  h_med = signif(median(d$h), 3))
print(do.call(rbind, lapply(split(x, paste(x$gen, x$coef)), function(d) cbind(d[1, c("gen", "coef")], f(d)))), row.names = FALSE)
saveRDS(x, "chk10.rds")

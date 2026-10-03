## Replicate Cao, Zeng & Fine (2015) Table 1 (time-invariant coefficient)
## with kee_async(). Design as stated in their Section 4.1: 1000 data sets,
## n = 100, 400, 900; Poisson(5) response and covariate counts; Unif(0, 1)
## times; Gaussian X with corr exp(-|s - t|); Gaussian error with
## cov 2^-|s - t|; beta = (0.5, 1.5); logistic Y = I(U <= expit(eta));
## Epanechnikov kernel; fixed h = n^-0.5, n^-0.6, n^-0.8, n^-0.9.
## sim_async_data() matches this except that it draws max(1, Poisson(5))
## observations, so no subject has zero observations of either kind.
## Their "auto" rows are not replicated.
## Usage: Rscript replicate_czf_table1.R <cores> <out.rds>
args <- commandArgs(TRUE); cores <- as.integer(args[1]); out <- args[2]
suppressMessages(library(skmle))
gam <- c(0.5, 0.6, 0.8, 0.9); z <- qnorm(0.975)
one <- function(job) {
  n <- job$n; link <- job$link; rep_id <- job$rep
  set.seed(7e6 + 1e4 * (link == "logistic") + 10 * n + rep_id)
  d <- sim_async_data(n = n, beta = c(0.5, 1.5), link = link)
  t(vapply(gam, function(g) {
    f <- tryCatch(suppressWarnings(kee_async(d$y, d$x, y ~ x, id = id, time = time,
                                             h = n^(-g), link = link)),
                  error = function(e) NULL)
    if (is.null(f) || f$convergence == 1L) return(c(NA_real_, NA_real_))
    c(coef(f)[2], sqrt(diag(vcov(f)))[2])
  }, numeric(2)))
}
jobs <- expand.grid(rep = 1:1000, n = c(100L, 400L, 900L), link = c("identity", "logistic"),
                    stringsAsFactors = FALSE)
res <- parallel::mclapply(split(jobs, seq_len(nrow(jobs))), one, mc.cores = cores)
tab <- list()
for (key in unique(paste(jobs$link, jobs$n))) {
  idx <- which(paste(jobs$link, jobs$n) == key)
  for (k in seq_along(gam)) {
    est <- sapply(res[idx], function(m) m[k, 1]); se <- sapply(res[idx], function(m) m[k, 2])
    ok <- is.finite(est) & is.finite(se); e <- est[ok] - 1.5; cp <- mean(abs(e) <= z * se[ok])
    tab[[length(tab) + 1]] <- data.frame(link = jobs$link[idx[1]], n = jobs$n[idx[1]],
      BD = sprintf("n^-%.1f", gam[k]), h = jobs$n[idx[1]]^(-gam[k]), N = sum(ok),
      fail = sum(!ok), Bias = mean(e), Bias_mcse = sd(e) / sqrt(sum(ok)),
      SD = sd(e), SE = mean(se[ok]), CP = 100 * cp, CP_mcse = 100 * sqrt(cp * (1 - cp) / sum(ok)))
  }
}
tab <- do.call(rbind, tab)
saveRDS(list(table = tab, raw = res, jobs = jobs), out)
num <- sapply(tab, is.numeric); tab[num] <- lapply(tab[num], signif, 3)
print(tab, row.names = FALSE)

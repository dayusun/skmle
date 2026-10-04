## Pilot: range of within-subject response-covariate lags |T_ij - R_ik|
## (exact zeros excluded) in sim_async_data(), per n. Sets the density of the
## per-data-set grid and the common reference grid of the main study.
## Usage: Rscript pilot_lag_range.R <reps per n>
reps <- as.integer(commandArgs(TRUE)[1])
suppressMessages(library(skmle))
lag_range <- function(d) {
  xs <- split(d$x$time, d$x$id); ys <- split(d$y$time, d$y$id)
  ids <- intersect(names(xs), names(ys))
  lg <- abs(unlist(lapply(ids, function(i) outer(ys[[i]], xs[[i]], "-")), use.names = FALSE))
  lg <- lg[lg > 0]
  c(min = min(lg), max = max(lg), npairs = length(lg))
}
for (n in c(100, 200, 400, 900, 2000, 4000)) {
  r <- t(sapply(seq_len(reps), function(k) { set.seed(k + 10 * n); lag_range(sim_async_data(n = n)) }))
  cat(sprintf("n=%4d  h_min q01/med/q99 %.1e %.1e %.1e | h_max q01/med %.3f %.3f | pairs med %d | log-range med %.1f\n",
              n, quantile(r[, 1], .01), median(r[, 1]), quantile(r[, 1], .99),
              quantile(r[, 2], .01), median(r[, 2]), as.integer(median(r[, 3])),
              median(log(r[, 2] / r[, 1]))))
}

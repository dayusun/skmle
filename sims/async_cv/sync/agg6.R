## Aggregate results6/: grid floor sqrt(n_train), five links.
source("cov5.R", echo = FALSE)
L <- unlist(lapply(Sys.glob("results6/*.rds"), readRDS), recursive = FALSE)
cat("errors:", sum(vapply(L, is.character, TRUE)), "\n")
x <- do.call(rbind, Filter(is.data.frame, L))
x$truth <- vapply(x$link, function(l) BETA[[l]][2], 0)
x$arm <- paste0(x$h_type, ifelse(x$corr, "+corr", ""))
cat("max pkg_diff (log):", max(x$pkg_diff, na.rm = TRUE), "\n")
summ <- function(d) {
  ok <- d$conv & is.finite(d$b1) & is.finite(d$se1); R <- sum(ok); b <- d$b1[ok]; s <- d$se1[ok]; tr <- d$truth[1]
  k <- sum(abs(b - tr) <= qnorm(.975) * s); bt <- binom.test(k, R, 0.95)
  data.frame(R = R, miss = nrow(d) - R, med_h = median(d$h[ok]), bias = mean(b) - tr,
    bias_se = sd(b) / sqrt(R), ESE = sd(b), SEmed = median(s), SE_ESE = median(s) / sd(b),
    RMSE = sqrt(mean((b - tr)^2)), CP = 100 * k / R, CP_lo = 100 * bt$conf.int[1],
    CP_hi = 100 * bt$conf.int[2], p = bt$p.value)
}
t <- do.call(rbind, lapply(split(x, paste(x$link, x$n, x$arm)), function(d)
  cbind(d[1, c("link", "n", "arm")], summ(d))))
t$q <- p.adjust(t$p, "BH"); t$sig <- as.character(cut(t$q, c(-1, .001, .01, .05, 2), c("***", "**", "*", "")))
t <- t[order(t$link, t$n, t$arm), ]; rownames(t) <- NULL
saveRDS(t, "summary6.rds"); write.csv(t, "summary6.csv", row.names = FALSE)
print(t, digits = 3)
s <- x[x$h_type %in% c("cv", "cvc") & x$corr, ]
s$at_floor <- s$h < s$h_lo * 1.0001; s$wild <- abs(s$b1 - s$truth) > 0.3
print(aggregate(cbind(at_floor = 100 * s$at_floor, wild = 100 * s$wild, h_lo = s$h_lo) ~ link + n + h_type, s, mean), digits = 3)

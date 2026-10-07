## Aggregate results7/: multivariate covariates (x1, x2 time-varying, z time-invariant).
L <- unlist(lapply(Sys.glob("results7/*.rds"), readRDS), recursive = FALSE)
cat("errors:", sum(vapply(L, is.character, TRUE)), "\n")
x <- do.call(rbind, Filter(is.data.frame, L)); x <- x[x$coef != "(Intercept)", ]
x$arm <- paste0(x$h_type, ifelse(x$corr, "+corr", ""))
summ <- function(d) {
  ok <- d$conv & is.finite(d$est) & is.finite(d$se); R <- sum(ok); b <- d$est[ok]; s <- d$se[ok]; tr <- d$truth[1]
  k <- sum(abs(b - tr) <= qnorm(.975) * s); bt <- binom.test(k, R, 0.95)
  data.frame(R = R, miss = nrow(d) - R, med_h = median(d$h[ok]), bias = mean(b) - tr, bias_se = sd(b) / sqrt(R),
    ESE = sd(b), SE_ESE = median(s) / sd(b), RMSE = sqrt(mean((b - tr)^2)), CP = 100 * k / R,
    CP_lo = 100 * bt$conf.int[1], CP_hi = 100 * bt$conf.int[2], p = bt$p.value)
}
t <- do.call(rbind, lapply(split(x, paste(x$link, x$n, x$arm, x$coef)), function(d)
  cbind(d[1, c("link", "n", "arm", "coef")], summ(d))))
t$q <- p.adjust(t$p, "BH"); t$sig <- as.character(cut(t$q, c(-1, .001, .01, .05, 2), c("***", "**", "*", "")))
t <- t[order(t$link, t$n, t$coef, t$arm), ]; rownames(t) <- NULL
saveRDS(t, "summary7.rds"); write.csv(t, "summary7.csv", row.names = FALSE)
options(width = 200); print(t, digits = 3)
s <- x[x$h_type %in% c("cv", "cvc") & x$corr & x$coef == "x1", ]
print(aggregate(cbind(at_floor = 100 * (s$h < s$h_lo * 1.0001)) ~ link + n + h_type, s, mean))

## Aggregate results12/: pair-moment estimator, three covariate scenarios.
L <- unlist(lapply(Sys.glob("results12/*.rds"), readRDS), recursive = FALSE)
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
t <- do.call(rbind, lapply(split(x, paste(x$gen, x$link, x$n, x$arm, x$coef)), function(d)
  cbind(d[1, c("gen", "link", "n", "arm", "coef")], summ(d))))
t$q <- p.adjust(t$p, "BH"); t$sig <- as.character(cut(t$q, c(-1, .001, .01, .05, 2), c("***", "**", "*", "")))
t <- t[order(t$gen, t$link, t$n, t$coef, t$arm), ]; rownames(t) <- NULL
saveRDS(t, "summary12.rds"); write.csv(t, "summary12.csv", row.names = FALSE)
options(width = 200); print(t, digits = 3)
s <- x[x$h_type %in% c("cv", "pm") & x$coef == "x1", ]
print(aggregate(cbind(at_floor = 100 * (s$h < s$h_lo * 1.0001)) ~ gen + link + n + h_type, s, mean))

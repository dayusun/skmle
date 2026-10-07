## Pair-moment estimator: no distributional model for x.
## U_h(b) = sum_TS w X_S Y_T / sum_TS w  -  sum_{k != k'} w X_k g(X_k' b) / sum_SS w,
## with w = w_h(lag) on response-covariate pairs and on ordered covariate pairs.
source("cov7.R")

## Non-Gaussian covariates: c_j = (OU_j^2 - 1)/sqrt(2) (stationary, skewed);
## x = A z + L c(t).
sim_skew <- function(n, beta, link) {
  rgp <- function(tt, f) as.numeric(t(chol(outer(tt, tt, f) + diag(1e-8, length(tt)))) %*% rnorm(length(tt)))
  L <- matrix(c(1, 0.6, 0, 0.8), 2); A <- c(0.5, -0.3)
  yl <- xl <- vector("list", n)
  for (i in seq_len(n)) {
    nT <- max(1L, rpois(1L, 5)); nS <- max(1L, rpois(1L, 5))
    Tt <- sort(runif(nT)); Ss <- sort(runif(nS)); at <- c(Tt, Ss); z <- rbinom(1, 1, 0.5)
    cc <- cbind(rgp(at, function(s, t) exp(-abs(s - t))), rgp(at, function(s, t) exp(-3 * abs(s - t))))
    cc <- (cc^2 - 1) / sqrt(2)
    X <- sweep(cc %*% t(L), 2, A * z, "+")
    eta <- beta[1] + X[seq_len(nT), , drop = FALSE] %*% beta[2:3] + beta[4] * z
    yv <- switch(link, identity = eta + rgp(Tt, function(s, t) 2^(-abs(s - t))),
      log = rpois(nT, exp(eta)), rbinom(nT, 1L, LINKS[[link]]$mu(eta)))
    yl[[i]] <- data.frame(id = i, time = Tt, y = as.numeric(yv))
    xl[[i]] <- data.frame(id = i, time = Ss, x1 = X[nT + seq_len(nS), 1], x2 = X[nT + seq_len(nS), 2], z = z)
  }
  list(y = do.call(rbind, yl), x = do.call(rbind, xl))
}

## Ordered within-subject covariate pairs: lag S_k' - S_k, subject, rows k (instrument), k' (stand-in for T).
cov_pairs <- function(p) do.call(rbind, lapply(split(seq_along(p$xi), p$xi), function(ix) {
  if (length(ix) < 2) return(NULL); g <- expand.grid(a = ix, b = ix); g <- g[g$a != g$b, ]
  cbind(v = p$xt[g$b] - p$xt[g$a], sid = p$xi[g$a], a = g$a, b = g$b) }))

prep8 <- function(dat) { e <- prep7(dat); e$Q <- cov_pairs(e$p); e }

fit_pm <- function(e, P, Q, h, link, start) {
  w <- epan(P[, "u"], h); kp <- w > 0; P <- P[kp, , drop = FALSE]; w <- w[kp]
  v <- epan(Q[, "v"], h); kq <- v > 0; Q <- Q[kq, , drop = FALSE]; v <- v[kq]
  if (!nrow(P) || !nrow(Q)) return(NULL)
  Lk <- LINKS[[link]]; n <- e$p$n
  XS <- e$X[P[, "k"], , drop = FALSE]; Xa <- e$X[Q[, "a"], , drop = FALSE]; Xb <- e$X[Q[, "b"], , drop = FALSE]
  WA <- sum(w); WB <- sum(v)
  full <- function(M, g) { out <- matrix(0, n, ncol(M)); r <- rowsum(M, g); out[as.integer(rownames(r)), ] <- r; out }
  Ai <- full(w * XS * P[, "y"], P[, "sid"])
  wa <- tapply(w, factor(P[, "sid"], levels = seq_len(n)), sum); wa[is.na(wa)] <- 0
  wb <- tapply(v, factor(Q[, "sid"], levels = seq_len(n)), sum); wb[is.na(wb)] <- 0
  A <- colSums(Ai) / WA
  ef <- function(b) { eta <- as.vector(Xb %*% b)
    Bi <- full(v * Xa * Lk$mu(eta), Q[, "sid"])
    list(Bi = Bi, U = A - colSums(Bi) / WB, J = -crossprod(Xa * v * Lk$d(eta), Xb) / WB) }
  b <- start; st <- Inf
  for (it in 1:50) { ee <- ef(b); st <- tryCatch(solve(ee$J, ee$U), error = function(err) NULL)
    if (is.null(st)) return(NULL); b <- b - st
    if (!all(is.finite(b))) return(NULL); if (max(abs(st)) < 1e-10) break }
  ee <- ef(b); B <- colSums(ee$Bi) / WB
  ## per-subject influence of the two ratio estimators
  psi <- (Ai - outer(wa, A)) / WA - (ee$Bi - outer(wb, B)) / WB
  Ji <- solve(ee$J)
  list(beta = b, V = Ji %*% crossprod(psi) %*% t(Ji), conv = max(abs(st)) < 1e-6)
}

## Same as sim_skew, but covariate times can be front-loaded: S ~ Beta(1, 3).
sim_skew2 <- function(n, beta, link, early = FALSE) {
  rgp <- function(tt, f) as.numeric(t(chol(outer(tt, tt, f) + diag(1e-8, length(tt)))) %*% rnorm(length(tt)))
  L <- matrix(c(1, 0.6, 0, 0.8), 2); A <- c(0.5, -0.3)
  yl <- xl <- vector("list", n)
  for (i in seq_len(n)) {
    nT <- max(1L, rpois(1L, 5)); nS <- max(1L, rpois(1L, 5))
    Tt <- sort(runif(nT)); Ss <- sort(if (early) rbeta(nS, 1, 3) else runif(nS)); at <- c(Tt, Ss); z <- rbinom(1, 1, 0.5)
    cc <- cbind(rgp(at, function(s, t) exp(-abs(s - t))), rgp(at, function(s, t) exp(-3 * abs(s - t))))
    cc <- (cc^2 - 1) / sqrt(2)
    X <- sweep(cc %*% t(L), 2, A * z, "+")
    eta <- beta[1] + X[seq_len(nT), , drop = FALSE] %*% beta[2:3] + beta[4] * z
    yv <- switch(link, identity = eta + rgp(Tt, function(s, t) 2^(-abs(s - t))),
      log = rpois(nT, exp(eta)), rbinom(nT, 1L, LINKS[[link]]$mu(eta)))
    yl[[i]] <- data.frame(id = i, time = Tt, y = as.numeric(yv))
    xl[[i]] <- data.frame(id = i, time = Ss, x1 = X[nT + seq_len(nS), 1], x2 = X[nT + seq_len(nS), 2], z = z)
  }
  list(y = do.call(rbind, yl), x = do.call(rbind, xl))
}

## Intensity ratio lambda_T / lambda_S at each covariate time, from the pooled
## observation times alone: logistic regression of "is a response occasion" on a
## natural spline of time (4 df). Only times enter, never x or y.
ratio_hat <- function(p, df = 4) {
  tt <- c(p$yt, p$xt); lab <- rep(1:0, c(length(p$yt), length(p$xt)))
  B <- splines::ns(tt, df = df); f <- glm.fit(cbind(1, B), lab, family = binomial())
  pr <- plogis(cbind(1, predict(B, p$xt)) %*% f$coefficients)
  as.vector(pr / (1 - pr))
}

## fit_pm with covariate-pair weights r[b] (r indexed by covariate row).
fit_pmw <- function(e, P, Q, h, link, start, r = NULL) {
  if (!is.null(r)) { Q <- cbind(Q, rw = r[Q[, "b"]]) } else Q <- cbind(Q, rw = 1)
  w <- epan(P[, "u"], h); kp <- w > 0; P <- P[kp, , drop = FALSE]; w <- w[kp]
  v <- epan(Q[, "v"], h) * Q[, "rw"]; kq <- v > 0; Q <- Q[kq, , drop = FALSE]; v <- v[kq]
  if (!nrow(P) || !nrow(Q)) return(NULL)
  Lk <- LINKS[[link]]; n <- e$p$n
  XS <- e$X[P[, "k"], , drop = FALSE]; Xa <- e$X[Q[, "a"], , drop = FALSE]; Xb <- e$X[Q[, "b"], , drop = FALSE]
  WA <- sum(w); WB <- sum(v)
  full <- function(M, g) { out <- matrix(0, n, ncol(M)); r <- rowsum(M, g); out[as.integer(rownames(r)), ] <- r; out }
  Ai <- full(w * XS * P[, "y"], P[, "sid"])
  wa <- full(cbind(w), P[, "sid"])[, 1]; wb <- full(cbind(v), Q[, "sid"])[, 1]
  A <- colSums(Ai) / WA
  ef <- function(b) { eta <- as.vector(Xb %*% b); Bi <- full(v * Xa * Lk$mu(eta), Q[, "sid"])
    list(Bi = Bi, U = A - colSums(Bi) / WB, J = -crossprod(Xa * v * Lk$d(eta), Xb) / WB) }
  b <- start; st <- Inf
  for (it in 1:50) { ee <- ef(b); st <- tryCatch(solve(ee$J, ee$U), error = function(err) NULL)
    if (is.null(st)) return(NULL); b <- b - st
    if (!all(is.finite(b))) return(NULL); if (max(abs(st)) < 1e-10) break }
  ee <- ef(b); B <- colSums(ee$Bi) / WB
  psi <- (Ai - outer(wa, A)) / WA - (ee$Bi - outer(wb, B)) / WB
  Ji <- solve(ee$J)
  list(beta = b, V = Ji %*% crossprod(psi) %*% t(Ji), conv = max(abs(st)) < 1e-6)
}

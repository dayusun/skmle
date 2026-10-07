## Pair-moment with (i) the intensity-ratio model stacked into the sandwich and
## (ii) the saturated-loss rule in the CV.
source("cov10.R")

## Logistic model for "is a response occasion" on ns(time, 4); returns the
## ratio at covariate rows plus what the stacked sandwich needs.
ratio_fit <- function(p, df = 4) {
  tt <- c(p$yt, p$xt); lab <- rep(1:0, c(length(p$yt), length(p$xt))); sid <- c(p$yi, p$xi)
  B <- splines::ns(tt, df = df); D <- cbind(1, B); f <- glm.fit(D, lab, family = binomial())
  pr <- as.vector(plogis(D %*% f$coefficients))
  Dx <- D[-seq_along(p$yt), , drop = FALSE]; px <- pr[-seq_along(p$yt)]
  S <- matrix(0, p$n, ncol(D)); rr <- rowsum(D * (lab - pr), sid); S[as.integer(rownames(rr)), ] <- rr
  list(r = px / (1 - px), Dx = Dx, S = S, Iinv = solve(crossprod(D * (pr * (1 - pr)), D)))
}

fit_pms <- function(e, P, Q, h, link, start, rf) {
  r <- rf$r; Q <- cbind(Q, rw = r[Q[, "b"]])
  w <- epan(P[, "u"], h); kp <- w > 0; P <- P[kp, , drop = FALSE]; w <- w[kp]
  v <- epan(Q[, "v"], h) * Q[, "rw"]; kq <- v > 0; Q <- Q[kq, , drop = FALSE]; v <- v[kq]
  if (!nrow(P) || !nrow(Q)) return(NULL)
  Lk <- LINKS[[link]]; n <- e$p$n
  XS <- e$X[P[, "k"], , drop = FALSE]; Xa <- e$X[Q[, "a"], , drop = FALSE]; Xb <- e$X[Q[, "b"], , drop = FALSE]
  Db <- rf$Dx[Q[, "b"], , drop = FALSE]                     # d log r / d theta at the stand-in time
  WA <- sum(w); WB <- sum(v)
  full <- function(M, g) { out <- matrix(0, n, ncol(M)); rr <- rowsum(M, g); out[as.integer(rownames(rr)), ] <- rr; out }
  Ai <- full(w * XS * P[, "y"], P[, "sid"]); wa <- full(cbind(w), P[, "sid"])[, 1]; wb <- full(cbind(v), Q[, "sid"])[, 1]
  A <- colSums(Ai) / WA
  ef <- function(b) { eta <- as.vector(Xb %*% b); Bi <- full(v * Xa * Lk$mu(eta), Q[, "sid"])
    list(Bi = Bi, eta = eta, U = A - colSums(Bi) / WB, J = -crossprod(Xa * v * Lk$d(eta), Xb) / WB) }
  b <- start; st <- Inf
  for (it in 1:50) { ee <- ef(b); st <- tryCatch(solve(ee$J, ee$U), error = function(err) NULL)
    if (is.null(st)) return(NULL); b <- b - st
    if (!all(is.finite(b))) return(NULL); if (max(abs(st)) < 1e-10) break }
  ee <- ef(b); B <- colSums(ee$Bi) / WB
  psi <- (Ai - outer(wa, A)) / WA - (ee$Bi - outer(wb, B)) / WB
  ## dU/dtheta = -d(B/WB)/dtheta; then add (dU/dtheta) Iinv s_i to each subject
  G <- Xa * Lk$mu(ee$eta)
  dB <- (crossprod(G * v, Db) - outer(B, colSums(v * Db))) / WB
  psi2 <- psi - rf$S %*% t(rf$Iinv) %*% t(dB)
  Ji <- solve(ee$J)
  list(beta = b, V = Ji %*% crossprod(psi2) %*% t(Ji), V0 = Ji %*% crossprod(psi) %*% t(Ji), conv = max(abs(st)) < 1e-6)
}

## Other links, correct model: log (Poisson), probit, cloglog (Bernoulli).
## Same estimating equation as the package, U = sum w X {y - mu}, with the corrected
## version replacing mu(eta_S) by E[mu(eta_T) | x_S] under a stationary Gaussian x
## (30-node Gauss-Hermite). Plain and corrected fold fits both use fit_g, so h-hat_cv
## and h-hat_cvc come from the same code path.
source("cov4.R")

LINKS <- list(
  identity = list(mu = function(e) e, d = function(e) rep(1, length(e))),
  logistic = list(mu = plogis, d = function(e) { p <- plogis(e); p * (1 - p) }),
  log      = list(mu = exp, d = exp),
  probit   = list(mu = pnorm, d = dnorm),
  cloglog  = list(mu = function(e) -expm1(-exp(e)), d = function(e) exp(e - exp(e))))
BETA <- list(identity = c(0.5, 1.5), logistic = c(0.3, 1.0), log = c(0.3, 0.5),
             probit = c(0.2, 0.6), cloglog = c(-0.5, 0.6))
LCODE <- c(identity = 0L, logistic = 2L, log = 1L, probit = 3L, cloglog = 4L)

## Held-out loss: squared error, Poisson deviance (up to y-only terms), Bernoulli log loss.
loss_of <- function(link) switch(link,
  identity = function(y, eta) (y - eta)^2,
  log = function(y, eta) exp(eta) - y * eta,
  function(y, eta) { p <- pmin(pmax(LINKS[[link]]$mu(eta), 1e-12), 1 - 1e-12)
    -y * log(p) - (1 - y) * log1p(-p) })

fit_g <- function(P, h, link, corr, gam, m, s2, start) {
  L <- LINKS[[link]]
  w <- epan(P[, "u"], h); k <- w > 0; P <- P[k, , drop = FALSE]; w <- w[k]
  X <- cbind(1, P[, "x"]); y <- P[, "y"]
  if (corr) { rho <- pmin(pmax(1 - gam(abs(P[, "u"])) / (2 * s2), -1), 1)
    cm <- m + rho * (P[, "x"] - m); sd <- sqrt(s2 * (1 - rho^2)) }
  ef <- function(b) {
    if (!corr) { eta <- as.vector(X %*% b); mu <- L$mu(eta); D <- X * L$d(eta) }
    else { mu <- 0; D <- matrix(0, nrow(P), 2)
      for (q in seq_along(gh$z)) { xx <- cm + sd * gh$z[q]; eta <- b[1] + b[2] * xx
        mu <- mu + gh$w[q] * L$mu(eta); D <- D + gh$w[q] * L$d(eta) * cbind(1, xx) } }
    list(U = w * X * (y - mu), J = -crossprod(X * w, D))
  }
  b <- start; st <- Inf
  for (it in 1:50) { e <- ef(b); st <- solve(e$J, colSums(e$U)); b <- b - st
    if (!all(is.finite(b))) break; if (max(abs(st)) < 1e-10) break }
  e <- ef(b); Ui <- rowsum(e$U, P[, "sid"]); Ai <- solve(e$J)
  list(beta = as.vector(b), V = Ai %*% crossprod(Ui) %*% t(Ai), conv = all(is.finite(b)) && max(abs(st)) < 1e-6)
}

## sim_async_data with the response drawn under any of the links (same RNG order).
sim_link <- function(n, beta, link) {
  rgp <- function(tt, f) as.numeric(t(chol(outer(tt, tt, f) + diag(1e-8, length(tt)))) %*% rnorm(length(tt)))
  yl <- xl <- vector("list", n)
  for (i in seq_len(n)) {
    nT <- max(1L, rpois(1L, 5)); nS <- max(1L, rpois(1L, 5))
    Tt <- sort(runif(nT)); Ss <- sort(runif(nS)); at <- c(Tt, Ss)
    Xa <- rgp(at, function(s, t) exp(-abs(s - t)))
    eta <- beta[1] + beta[2] * Xa[seq_len(nT)]
    yv <- switch(link, identity = eta + rgp(Tt, function(s, t) 2^(-abs(s - t))),
      log = rpois(nT, exp(eta)), rbinom(nT, 1L, LINKS[[link]]$mu(eta)))
    yl[[i]] <- data.frame(id = i, time = Tt, y = as.numeric(yv))
    xl[[i]] <- data.frame(id = i, time = Ss, x = Xa[nT + seq_len(nS)])
  }
  list(y = do.call(rbind, yl), x = do.call(rbind, xl))
}

## Sync-loss CV with plain (cv) and corrected (cvc) fold fits.
select_g <- function(e, link, beta0) {
  P <- e$P; pl <- loss_of(link)
  lg <- abs(P[, "u"]); lg <- lg[lg > 0]
  grid <- exp(seq(log(min(lg)), log(max(lg)), length.out = NG))
  folds <- ((e$ids - 1L) %% K) + 1L
  IC <- ICc <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) {
    H <- P[P[, "sid"] %in% which(folds == f), , drop = FALSE]
    Ptr <- P[P[, "sid"] %in% which(folds != f), , drop = FALSE]
    QR <- qr(cbind(1, e$gam(abs(H[, "u"]))))
    for (corr in c(FALSE, TRUE)) { st <- beta0
      for (g in rev(seq_len(NG))) {          # warm start from large h down
        ft <- tryCatch(fit_g(Ptr, grid[g], link, corr, e$gam, e$m, e$s2, st), error = function(e) NULL)
        if (is.null(ft) || !ft$conv) next
        st <- ft$beta
        v <- qr.coef(QR, pl(H[, "y"], ft$beta[1] + ft$beta[2] * H[, "x"]))[1]
        if (corr) ICc[g, f] <- v else IC[g, f] <- v
      } }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; cvc <- rowMeans(ICc); cvc[!is.finite(cvc)] <- Inf
  c(cv = grid[which.min(cv)], cvc = grid[which.min(cvc)])
}

one_rep_link <- function(n, link, rep_id) {
  beta0 <- BETA[[link]]
  set.seed(2e7 + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep(sim_link(n, beta0, link))
  hs <- select_g(e, link, beta0)
  out <- NULL
  for (a in list(c("cv", FALSE), c("cv", TRUE), c("cvc", TRUE))) {
    h <- hs[[a[1]]]; corr <- as.logical(a[2])
    r <- tryCatch(fit_g(e$P, h, link, corr, e$gam, e$m, e$s2, beta0), error = function(e) NULL)
    out <- rbind(out, row_of(a[1], h, corr, r))
  }
  PL <- lvcf_pairs(e$p, e$xv)
  r <- tryCatch(fit_g(PL, 1e6, link, FALSE, e$gam, e$m, e$s2, beta0), error = function(e) NULL)
  out <- rbind(out, row_of("lvcf", NA, FALSE, r))
  ## the plain fit equals the package fit where the package has the link
  pd <- NA_real_
  if (link %in% c("identity", "log", "logistic")) {
    pk <- fit_(e$p, hs[["cv"]], FALSE, LCODE[[link]], 50L, 1e-10)$beta
    pd <- max(abs(pk - c(out$b0[1], out$b1[1]))) }
  cbind(link = link, n = n, rep = rep_id, out, pkg_diff = pd)
}

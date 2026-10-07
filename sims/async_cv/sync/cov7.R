## Multivariate covariates: q_v time-varying (jointly stationary, symmetric
## cross-covariance, Gaussian for nonlinear links) plus time-invariant ones.
## Correction: x_T | x_S, z ~ N(m(z) + R(u)(x_S - m(z)), Sigma - C Sigma^-1 C),
## C(u) = Sigma - Gamma(u)/2, Gamma by polarised isotonic variograms.
## CV: held-out pair losses regressed on (1, vech Gamma-hat(u)); intercept.
source("cov2.R")   # gh, epan, NG, K, skmle loaded

LINKS <- list(
  identity = list(mu = function(e) e, d = function(e) rep(1, length(e))),
  logistic = list(mu = plogis, d = function(e) { p <- plogis(e); p * (1 - p) }),
  log      = list(mu = exp, d = exp),
  probit   = list(mu = pnorm, d = dnorm))
LCODE <- c(identity = 0L, logistic = 2L, log = 1L, probit = 3L)
BETA <- list(identity = c(0.5, 1.0, -0.8, 0.5), logistic = c(0.3, 0.8, -0.6, 0.4),
             log = c(0.2, 0.4, -0.3, 0.3), probit = c(0.2, 0.5, -0.4, 0.3))
loss_of <- function(link) switch(link,
  identity = function(y, eta) (y - eta)^2,
  log = function(y, eta) exp(eta) - y * eta,
  function(y, eta) { p <- pmin(pmax(LINKS[[link]]$mu(eta), 1e-12), 1 - 1e-12)
    -y * log(p) - (1 - y) * log1p(-p) })

## x(t) = A z + L c(t), c = two independent OU paths with rates 1 and 3:
## C(u) = L diag(e^-u, e^-3u) L' is symmetric and not separable.
sim_multi <- function(n, beta, link) {
  rgp <- function(tt, f) as.numeric(t(chol(outer(tt, tt, f) + diag(1e-8, length(tt)))) %*% rnorm(length(tt)))
  L <- matrix(c(1, 0.6, 0, 0.8), 2); A <- c(0.5, -0.3)
  yl <- xl <- vector("list", n)
  for (i in seq_len(n)) {
    nT <- max(1L, rpois(1L, 5)); nS <- max(1L, rpois(1L, 5))
    Tt <- sort(runif(nT)); Ss <- sort(runif(nS)); at <- c(Tt, Ss); z <- rbinom(1, 1, 0.5)
    cc <- cbind(rgp(at, function(s, t) exp(-abs(s - t))), rgp(at, function(s, t) exp(-3 * abs(s - t))))
    X <- sweep(cc %*% t(L), 2, A * z, "+")
    eta <- beta[1] + X[seq_len(nT), , drop = FALSE] %*% beta[2:3] + beta[4] * z
    yv <- switch(link, identity = eta + rgp(Tt, function(s, t) 2^(-abs(s - t))),
      log = rpois(nT, exp(eta)), rbinom(nT, 1L, LINKS[[link]]$mu(eta)))
    yl[[i]] <- data.frame(id = i, time = Tt, y = as.numeric(yv))
    xl[[i]] <- data.frame(id = i, time = Ss, x1 = X[nT + seq_len(nS), 1], x2 = X[nT + seq_len(nS), 2], z = z)
  }
  list(y = do.call(rbind, yl), x = do.call(rbind, xl))
}

## Isotonic variogram of a scalar increment series; returns values at lags `at`.
iso_at <- function(L, d, at) { o <- order(L); yf <- stats::isoreg(L[o], d[o])$yf
  i <- findInterval(at, L[o]); out <- yf[pmax(i, 1L)]; out[i == 0L] <- 0; out }

prep7 <- function(dat) {
  p <- parse_(y ~ x1 + x2 + z, dat$y, dat$x, "id", "time", TRUE)
  X <- p$X; cols <- colnames(X)[-1]
  const <- vapply(cols, function(cn) all(tapply(X[, cn], p$xi, function(v) diff(range(v)) == 0)), TRUE)
  tv <- cols[!const]; ti <- cols[const]; qv <- length(tv)
  ## mean of the time-varying part given the time-invariant covariates
  Zm <- cbind(1, X[, ti, drop = FALSE]); Bm <- qr.solve(Zm, X[, tv, drop = FALSE])
  M <- Zm %*% Bm; E <- X[, tv, drop = FALSE] - M; Sig <- crossprod(E) / nrow(E)
  ## within-subject covariate increments
  cp <- do.call(rbind, lapply(split(seq_along(p$xi), p$xi), function(ix) {
    if (length(ix) < 2) return(NULL); cb <- utils::combn(ix, 2)
    cbind(abs(p$xt[cb[1, ]] - p$xt[cb[2, ]]), X[cb[1, ], tv, drop = FALSE] - X[cb[2, ], tv, drop = FALSE]) }))
  lag <- cp[, 1]; D <- cp[, -1, drop = FALSE]
  ## breakpoints: every covariate-pair lag; Gamma is constant between them
  br <- sort(unique(lag)); at <- c(0, br)
  jk <- which(upper.tri(diag(qv), diag = TRUE), arr.ind = TRUE)
  diagv <- lapply(seq_len(qv), function(j) iso_at(lag, D[, j]^2, at))
  G <- array(0, c(length(at), qv, qv))
  for (r in seq_len(nrow(jk))) { j <- jk[r, 1]; k <- jk[r, 2]
    g <- if (j == k) diagv[[j]] else (iso_at(lag, (D[, j] + D[, k])^2, at) - diagv[[j]] - diagv[[k]]) / 2
    G[, j, k] <- G[, k, j] <- g }
  ## collapse to distinct Gamma values, then R and V per block
  key <- apply(matrix(G, length(at)), 1, paste, collapse = ",")
  blk <- cumsum(c(TRUE, key[-1] != key[-length(key)]))
  first <- which(!duplicated(blk)); Si <- solve(Sig)
  RV <- lapply(first, function(r) { C <- Sig - G[r, , ] / 2; R <- C %*% Si
    V <- Sig - C %*% Si %*% C; ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
    list(R = R, V = ev$vectors %*% diag(pmax(ev$values, 0), qv) %*% t(ev$vectors),
         g = G[r, , ][cbind(jk)]) })
  Rb <- simplify2array(lapply(RV, `[[`, "R")); Vb <- simplify2array(lapply(RV, `[[`, "V"))
  gb <- do.call(rbind, lapply(RV, `[[`, "g"))
  block_of <- function(u) blk[findInterval(u, at)]
  ## all within-subject pairs: lag, y, subject, covariate row
  ys <- split(seq_along(p$yi), factor(p$yi, levels = seq_len(p$n)))
  xs <- split(seq_along(p$xi), factor(p$xi, levels = seq_len(p$n)))
  P <- do.call(rbind, lapply(seq_len(p$n), function(i) { a <- ys[[i]]; b <- xs[[i]]
    if (!length(a) || !length(b)) return(NULL)
    cbind(u = as.vector(outer(p$yt[a], p$xt[b], "-")), y = rep(p$y[a], times = length(b)),
          sid = i, k = rep(b, each = length(a)), j = rep(a, times = length(b))) }))
  list(p = p, X = X, tv = tv, ti = ti, M = M, Rb = Rb, Vb = Vb, gb = gb,
       block_of = block_of, P = P, ids = seq_len(p$n))
}

## Per-pair conditional mean of the full design row at T and covariance of x_T.
cond_design <- function(e, P) {
  X <- e$X[P[, "k"], , drop = FALSE]; b <- e$block_of(abs(P[, "u"]))
  tv <- e$tv; qv <- length(tv); M <- e$M[P[, "k"], , drop = FALSE]
  Z <- X; dev <- X[, tv, drop = FALSE] - M
  for (j in seq_len(qv)) { zj <- M[, j]
    for (l in seq_len(qv)) zj <- zj + e$Rb[j, l, b] * dev[, l]
    Z[, tv[j]] <- zj }
  list(X = X, Z = Z, b = b)
}

fit7 <- function(e, P, h, link, corr, start, CD = NULL) {
  w <- epan(P[, "u"], h); kk <- w > 0; if (!any(kk)) return(NULL)
  P <- P[kk, , drop = FALSE]; w <- w[kk]; y <- P[, "y"]
  cd <- if (is.null(CD)) cond_design(e, P) else lapply(CD, function(a) if (is.matrix(a)) a[kk, , drop = FALSE] else a[kk])
  X <- cd$X; Z <- cd$Z; tvi <- match(e$tv, colnames(X)); qv <- length(tvi); Lk <- LINKS[[link]]
  ef <- function(bb) {
    if (!corr) { eta <- as.vector(X %*% bb); mu <- Lk$mu(eta); D <- X * Lk$d(eta) }
    else {
      bv <- bb[tvi]; mean <- as.vector(Z %*% bb)
      Vb <- matrix(0, nrow(X), qv)                  # per pair V b_v
      for (j in seq_len(qv)) for (l in seq_len(qv)) Vb[, j] <- Vb[, j] + e$Vb[j, l, cd$b] * bv[l]
      s2 <- pmax(as.vector(Vb %*% bv), 0)
      if (link == "identity") { mu <- mean; D <- Z }
      else if (link == "log") { mu <- exp(mean + s2 / 2); D <- Z * mu; D[, tvi] <- mu * (Z[, tvi] + Vb) }
      else { s <- sqrt(s2); mu <- 0; D <- 0 * X; dS <- Vb / ifelse(s > 0, s, 1)
        for (q in seq_along(gh$z)) { et <- mean + s * gh$z[q]; mu <- mu + gh$w[q] * Lk$mu(et)
          dq <- gh$w[q] * Lk$d(et); D <- D + dq * Z; D[, tvi] <- D[, tvi] + dq * gh$z[q] * dS } }
    }
    list(U = w * X * (y - mu), J = crossprod(X * w, D))
  }
  bb <- start; st <- Inf
  for (it in 1:50) { ee <- ef(bb); st <- tryCatch(solve(ee$J, colSums(ee$U)), error = function(err) NULL)
    if (is.null(st)) return(NULL); bb <- bb + st
    if (!all(is.finite(bb))) return(NULL); if (max(abs(st)) < 1e-10) break }
  ee <- ef(bb); Ui <- rowsum(ee$U, P[, "sid"]); Ai <- solve(ee$J)
  list(beta = bb, V = Ai %*% crossprod(Ui) %*% t(Ai), conv = max(abs(st)) < 1e-6)
}

floor7 <- function(P, n, folds) {
  mi <- tapply(abs(P[, "u"]), factor(P[, "sid"], levels = seq_len(n)), min); mi[is.na(mi)] <- Inf
  max(vapply(seq_len(K), function(f) { m <- sort(mi[folds != f]); m[min(ceiling(sqrt(sum(folds != f))), length(m))] }, 0)) * (1 + 1e-8)
}

select7 <- function(e, link, start) {
  P <- e$P; pl <- loss_of(link); folds <- ((e$ids - 1L) %% K) + 1L
  lg <- abs(P[, "u"]); lo <- max(min(lg[lg > 0]), floor7(P, e$p$n, folds))
  grid <- exp(seq(log(lo), log(max(lg)), length.out = NG))
  IC <- ICc <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) {
    H <- P[P[, "sid"] %in% which(folds == f), , drop = FALSE]
    Ptr <- P[P[, "sid"] %in% which(folds != f), , drop = FALSE]
    CDtr <- cond_design(e, Ptr); HX <- e$X[H[, "k"], , drop = FALSE]
    QR <- qr(cbind(1, e$gb[e$block_of(abs(H[, "u"])), , drop = FALSE]))
    for (corr in c(FALSE, TRUE)) { st <- start
      for (g in rev(seq_len(NG))) {
        ft <- tryCatch(fit7(e, Ptr, grid[g], link, corr, st, CDtr), error = function(err) NULL)
        if (is.null(ft) || !ft$conv) next
        st <- ft$beta
        v <- qr.coef(QR, pl(H[, "y"], as.vector(HX %*% ft$beta)))[1]
        if (corr) ICc[g, f] <- v else IC[g, f] <- v } }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; cvc <- rowMeans(ICc); cvc[!is.finite(cvc)] <- Inf
  c(cv = grid[which.min(cv)], cvc = grid[which.min(cvc)], lo = lo)
}

## LVCF pairs: each response with the latest covariate row at or before it.
lvcf7 <- function(e) {
  P <- e$P[e$P[, "u"] >= 0, , drop = FALSE]
  P <- P[order(P[, "j"], P[, "u"]), , drop = FALSE]
  P[!duplicated(P[, "j"]), , drop = FALSE]
}

row7 <- function(type, h, corr, r, nm) {
  b <- if (is.null(r)) rep(NA, length(nm)) else r$beta
  s <- if (is.null(r)) rep(NA, length(nm)) else sqrt(diag(r$V))
  data.frame(h_type = type, h = h, corr = corr, coef = nm, est = b, se = s,
             conv = !is.null(r) && isTRUE(r$conv) && all(is.finite(r$V)))
}

one_rep7 <- function(n, link, rep_id) {
  beta0 <- BETA[[link]]
  set.seed(3e7 + 1e5 * LCODE[[link]] + 1e3 * n + rep_id)
  e <- prep7(sim_multi(n, beta0, link)); nm <- colnames(e$X)
  st <- tryCatch(fit7(e, e$P, max(abs(e$P[, "u"])) * 1.01, link, FALSE, numeric(length(nm)))$beta,
                 error = function(err) numeric(length(nm)))
  hs <- select7(e, link, st); out <- NULL
  for (a in list(c("cv", FALSE), c("cv", TRUE), c("cvc", TRUE))) {
    h <- hs[[a[1]]]; corr <- as.logical(a[2])
    r <- tryCatch(fit7(e, e$P, h, link, corr, st), error = function(err) NULL)
    out <- rbind(out, row7(a[1], h, corr, r, nm)) }
  r <- tryCatch(fit7(e, lvcf7(e), 1e6, link, FALSE, st), error = function(err) NULL)
  out <- rbind(out, row7("lvcf", NA, FALSE, r, nm))
  cbind(link = link, n = n, rep = rep_id, truth = beta0[match(out$coef, nm)], out, h_lo = hs[["lo"]])
}

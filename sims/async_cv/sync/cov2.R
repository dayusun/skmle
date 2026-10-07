## Coverage at the synchronous-loss CV bandwidth.
## h-hat: sync-loss CV (isotonic variogram, K = 5 fixed folds, data-driven grid).
## At h-hat: (P) the kernel GEE with its sandwich; (C) the same equation with the
## asynchrony bias removed, which is unbiased at every lag under stationarity.
suppressMessages(library(skmle))
ns <- asNamespace("skmle")
parse_ <- get(".async_parse", ns); subset_ <- get(".async_subset", ns); fit_ <- get(".async_fit_ti", ns)
NG <- 100L; K <- 5L
gh <- local({ m <- 30; b <- sqrt(seq_len(m - 1) / 2); J <- matrix(0, m, m)
  J[cbind(1:(m - 1), 2:m)] <- b; J[cbind(2:m, 1:(m - 1))] <- b; e <- eigen(J, symmetric = TRUE)
  list(z = e$values * sqrt(2), w = e$vectors[1, ]^2) })           # E f(Z), Z ~ N(0,1)
vario_iso <- function(L, d) { o <- order(L); fit <- stats::isoreg(L[o], d[o]); xs <- L[o]; ys <- fit$yf
  function(u) { i <- findInterval(u, xs); out <- ys[pmax(i, 1L)]; out[i == 0L] <- 0; out } }
epan <- function(u, h) { z <- u / h; ifelse(abs(z) < 1, 0.75 * (1 - z^2) / h, 0) }

## all within-subject (response, covariate) pairs
pairs_of <- function(p, xv, ids) do.call(rbind, lapply(ids, function(i) {
  yi <- which(p$yi == i); xi <- which(p$xi == i); if (!length(yi) || !length(xi)) return(NULL)
  cbind(u = as.vector(outer(p$yt[yi], p$xt[xi], "-")), y = rep(p$y[yi], times = length(xi)),
        x = rep(xv[xi], each = length(yi)), sid = i) }))

## Solve sum_p w_p X_p {y_p - g_p(beta)} + corr = 0 by Newton; return beta and sandwich.
fit_pairs <- function(P, h, lk, corr, gam, m, s2, start) {
  w <- epan(P[, "u"], h); k <- w > 0; P <- P[k, , drop = FALSE]; w <- w[k]
  X <- cbind(1, P[, "x"]); y <- P[, "y"]; g <- gam(abs(P[, "u"]))
  rho <- pmin(pmax(1 - g / (2 * s2), -1), 1)
  ef <- function(b) {           # per-pair contributions (rows) and the Jacobian
    if (lk == 0L) {
      r <- y - X %*% b; U <- w * X * as.vector(r)
      J <- -crossprod(X * w, X)
      if (corr) { U[, 2] <- U[, 2] + w * b[2] * g / 2; J[2, 2] <- J[2, 2] + sum(w * g) / 2 }
    } else if (!corr) {
      mu <- plogis(as.vector(X %*% b)); U <- w * X * (y - mu); J <- -crossprod(X * w * mu * (1 - mu), X)
    } else {                    # E[mu(eta_T) | x_S] under a stationary Gaussian x
      cm <- m + rho * (P[, "x"] - m); sd <- sqrt(s2 * (1 - rho^2))
      E <- matrix(0, nrow(P), 3)
      for (q in seq_along(gh$z)) { xx <- cm + sd * gh$z[q]; mu <- plogis(b[1] + b[2] * xx)
        E <- E + gh$w[q] * cbind(mu, mu * (1 - mu), mu * (1 - mu) * xx) }
      U <- w * X * (y - E[, 1]); J <- -crossprod(X * w, E[, 2:3])
    }
    list(U = U, J = J)
  }
  b <- start
  for (it in 1:50) { e <- ef(b); st <- solve(e$J, colSums(e$U)); b <- b - st; if (max(abs(st)) < 1e-10) break }
  e <- ef(b); Ui <- rowsum(e$U, P[, "sid"]); Ai <- solve(e$J)
  list(beta = as.vector(b), V = Ai %*% crossprod(Ui) %*% t(Ai), conv = max(abs(st)) < 1e-6)
}

one_rep <- function(n, link, rep_id, h_ref = NULL) {
  lk <- c(identity = 0L, logistic = 2L)[[link]]
  beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
  pl <- function(y, eta) if (lk == 0L) (y - eta)^2 else log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta
  set.seed(1e5 * lk + 1e3 * n + rep_id)
  d <- sim_async_data(n = n, beta = beta0, link = link)
  p <- parse_(y ~ x, d$y, d$x, "id", "time", TRUE); xv <- p$X[, 2]
  ids <- seq_len(p$n); P <- pairs_of(p, xv, ids)
  lg <- abs(P[, "u"]); lg <- lg[lg > 0]
  grid <- exp(seq(log(min(lg)), log(max(lg)), length.out = NG))
  folds <- ((ids - 1L) %% K) + 1L
  cp <- do.call(rbind, lapply(split(seq_along(p$xt), p$xi), function(ix) {
    if (length(ix) < 2) return(NULL); cb <- utils::combn(ix, 2)
    cbind(abs(p$xt[cb[1, ]] - p$xt[cb[2, ]]), (xv[cb[1, ]] - xv[cb[2, ]])^2) }))
  gam <- vario_iso(cp[, 1], cp[, 2]); m <- mean(xv); s2 <- mean((xv - m)^2)
  IC <- ICc <- matrix(NA_real_, NG, K)
  for (f in seq_len(K)) {
    tr <- subset_(p, which(folds != f)); H <- P[P[, "sid"] %in% which(folds == f), , drop = FALSE]
    QR <- qr(cbind(1, gam(abs(H[, "u"]))))
    for (g in seq_len(NG)) {
      ft <- tryCatch(fit_(tr, grid[g], FALSE, lk, 50L, 1e-8), error = function(e) NULL)
      if (is.null(ft) || ft$convergence == 1L || !all(is.finite(ft$beta))) next
      IC[g, f] <- qr.coef(QR, pl(H[, "y"], ft$beta[1] + ft$beta[2] * H[, "x"]))[1]
    }
    Ptr <- P[P[, "sid"] %in% which(folds != f), , drop = FALSE]; st <- beta0
    for (g in rev(seq_len(NG))) {   # corrected fits, warm-started from large h down
      fc <- tryCatch(fit_pairs(Ptr, grid[g], lk, TRUE, gam, m, s2, st), error = function(e) NULL)
      if (is.null(fc) || !fc$conv || !all(is.finite(fc$beta))) next
      st <- fc$beta
      ICc[g, f] <- qr.coef(QR, pl(H[, "y"], fc$beta[1] + fc$beta[2] * H[, "x"]))[1]
    }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; hh <- grid[which.min(cv)]
  cvc <- rowMeans(ICc); cvc[!is.finite(cvc)] <- Inf; hc <- grid[which.min(cvc)]
  out <- NULL
  for (hs in c(cv = hh, ref = h_ref)) for (corr in c(FALSE, TRUE)) {
    r <- tryCatch(fit_pairs(P, hs, lk, corr, gam, m, s2, beta0), error = function(e) NULL)
    if (is.null(r)) { out <- rbind(out, data.frame(h_type = names(c(cv = hh, ref = h_ref))[c(hh, h_ref) == hs][1], h = hs, corr = corr, b1 = NA, se1 = NA, b0 = NA, se0 = NA, conv = FALSE)); next }
    out <- rbind(out, data.frame(h_type = if (identical(hs, hh)) "cv" else "ref", h = hs, corr = corr,
      b1 = r$beta[2], se1 = sqrt(r$V[2, 2]), b0 = r$beta[1], se0 = sqrt(r$V[1, 1]), conv = r$conv))
  }
  ## undersmoothing: h = h-hat * n^-c
  for (cc in c(1/10, 1/6, 1/4)) for (corr in c(FALSE, TRUE)) {
    hu <- hh * n^(-cc)
    r <- tryCatch(fit_pairs(P, hu, lk, corr, gam, m, s2, beta0), error = function(e) NULL)
    ok <- !is.null(r) && all(is.finite(r$V))
    out <- rbind(out, data.frame(h_type = sprintf("us%.3f", cc), h = hu, corr = corr,
      b1 = if (ok) r$beta[2] else NA, se1 = if (ok) sqrt(r$V[2, 2]) else NA,
      b0 = if (ok) r$beta[1] else NA, se0 = if (ok) sqrt(r$V[1, 1]) else NA, conv = ok && r$conv))
  }
  rc <- tryCatch(fit_pairs(P, hc, lk, TRUE, gam, m, s2, beta0), error = function(e) NULL)
  out <- rbind(out, data.frame(h_type = "cvc", h = hc, corr = TRUE,
    b1 = if (is.null(rc)) NA else rc$beta[2], se1 = if (is.null(rc)) NA else sqrt(rc$V[2, 2]),
    b0 = if (is.null(rc)) NA else rc$beta[1], se0 = if (is.null(rc)) NA else sqrt(rc$V[1, 1]),
    conv = !is.null(rc) && rc$conv))
  ## sanity: the plain pair-level fit equals the package fit at h-hat
  pk <- fit_(p, hh, FALSE, lk, 50L, 1e-10)$beta
  pr <- out[out$h_type == "cv" & !out$corr, ][1, ]
  cbind(link = link, n = n, rep = rep_id, out, pkg_diff = max(abs(pk - c(pr$b0, pr$b1))))
}

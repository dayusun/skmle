## (1) Oracle curves: plain and corrected fits on a common fixed-h grid, plus LVCF,
##     on the same data as the coverage run (same generator, same seeds).
## (3) Robustness: misspecified covariate processes; CV picks, fixed-h curves, LVCF.
source("cov2.R")

HFIX <- exp(seq(log(0.005), log(1), length.out = 40))

## Last value carried forward: each response paired with the latest covariate at or before it.
lvcf_pairs <- function(p, xv) do.call(rbind, lapply(seq_len(p$n), function(i) {
  yi <- which(p$yi == i); xi <- which(p$xi == i); if (!length(yi) || !length(xi)) return(NULL)
  xt <- p$xt[xi]; o <- order(xt); xt <- xt[o]; xi <- xi[o]
  k <- findInterval(p$yt[yi], xt); ok <- k > 0; if (!any(ok)) return(NULL)
  cbind(u = p$yt[yi][ok] - xt[k[ok]], y = p$y[yi][ok], x = xv[xi[k[ok]]], sid = i) }))

## Covariate processes. Z is the stationary OU path of sim_async_data.
sim_rob <- function(n, beta, link, scen) {
  rgp <- function(tt, f) as.numeric(t(chol(outer(tt, tt, f) + diag(1e-8, length(tt)))) %*% rnorm(length(tt)))
  yl <- xl <- vector("list", n)
  for (i in seq_len(n)) {
    nT <- max(1L, rpois(1L, 5)); nS <- max(1L, rpois(1L, 5))
    Tt <- sort(runif(nT)); Ss <- sort(runif(nS)); at <- c(Tt, Ss)
    Z <- rgp(at, function(s, t) exp(-abs(s - t)))
    X <- switch(scen,
      nonstat = 2 * at - 1 + (0.6 + 0.8 * at) * Z,      # mean and variance drift in t
      merr = Z,                                         # error added to observed x below
      skew = (Z^2 - 1) / sqrt(2))                       # stationary, chi-square marginal
    xs <- X[nT + seq_len(nS)]
    if (scen == "merr") xs <- xs + rnorm(nS, 0, 0.5)    # reliability 0.8
    eta <- beta[1] + beta[2] * X[seq_len(nT)]
    yv <- if (link == "identity") eta + rgp(Tt, function(s, t) 2^(-abs(s - t))) else rbinom(nT, 1L, plogis(eta))
    yl[[i]] <- data.frame(id = i, time = Tt, y = yv); xl[[i]] <- data.frame(id = i, time = Ss, x = xs)
  }
  list(y = do.call(rbind, yl), x = do.call(rbind, xl))
}

prep <- function(d) {
  p <- parse_(y ~ x, d$y, d$x, "id", "time", TRUE); xv <- p$X[, 2]
  ids <- seq_len(p$n); P <- pairs_of(p, xv, ids)
  cp <- do.call(rbind, lapply(split(seq_along(p$xt), p$xi), function(ix) {
    if (length(ix) < 2) return(NULL); cb <- utils::combn(ix, 2)
    cbind(abs(p$xt[cb[1, ]] - p$xt[cb[2, ]]), (xv[cb[1, ]] - xv[cb[2, ]])^2) }))
  list(p = p, xv = xv, ids = ids, P = P, gam = vario_iso(cp[, 1], cp[, 2]),
       m = mean(xv), s2 = mean((xv - mean(xv))^2))
}

## Sync-loss CV on plain fold fits (cv) and on corrected fold fits (cvc); same as one_rep.
select_h <- function(e, lk, beta0) {
  p <- e$p; P <- e$P; gam <- e$gam
  pl <- function(y, eta) if (lk == 0L) (y - eta)^2 else log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta
  lg <- abs(P[, "u"]); lg <- lg[lg > 0]
  grid <- exp(seq(log(min(lg)), log(max(lg)), length.out = NG))
  folds <- ((e$ids - 1L) %% K) + 1L
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
    for (g in rev(seq_len(NG))) {
      fc <- tryCatch(fit_pairs(Ptr, grid[g], lk, TRUE, gam, e$m, e$s2, st), error = function(e) NULL)
      if (is.null(fc) || !fc$conv || !all(is.finite(fc$beta))) next
      st <- fc$beta
      ICc[g, f] <- qr.coef(QR, pl(H[, "y"], fc$beta[1] + fc$beta[2] * H[, "x"]))[1]
    }
  }
  cv <- rowMeans(IC); cv[!is.finite(cv)] <- Inf; cvc <- rowMeans(ICc); cvc[!is.finite(cvc)] <- Inf
  c(cv = grid[which.min(cv)], cvc = grid[which.min(cvc)])
}

row_of <- function(type, h, corr, r) data.frame(h_type = type, h = h, corr = corr,
  b1 = if (is.null(r)) NA else r$beta[2], se1 = if (is.null(r)) NA else sqrt(r$V[2, 2]),
  b0 = if (is.null(r)) NA else r$beta[1], se0 = if (is.null(r)) NA else sqrt(r$V[1, 1]),
  conv = !is.null(r) && isTRUE(r$conv) && all(is.finite(r$V)))

fits <- function(e, lk, beta0, hs, types, corrs) {
  out <- NULL
  for (j in seq_along(hs)) {
    r <- tryCatch(fit_pairs(e$P, hs[j], lk, corrs[j], e$gam, e$m, e$s2, beta0), error = function(e) NULL)
    out <- rbind(out, row_of(types[j], hs[j], corrs[j], r))
  }
  PL <- lvcf_pairs(e$p, e$xv)       # constant weight: h far beyond every lag
  r <- tryCatch(fit_pairs(PL, 1e6, lk, FALSE, e$gam, e$m, e$s2, beta0), error = function(e) NULL)
  rbind(out, row_of("lvcf", NA, FALSE, r))
}

## (1) Same seed and generator call as one_rep(), so the data match the coverage run.
one_rep_oracle <- function(n, link, rep_id) {
  lk <- c(identity = 0L, logistic = 2L)[[link]]
  beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
  set.seed(1e5 * lk + 1e3 * n + rep_id)
  e <- prep(sim_async_data(n = n, beta = beta0, link = link))
  hs <- rep(HFIX, 2); cr <- rep(c(FALSE, TRUE), each = length(HFIX))
  cbind(link = link, n = n, rep = rep_id, scen = "base",
        fits(e, lk, beta0, hs, rep("fix", length(hs)), cr))
}

## (3) Misspecified covariate processes.
HROB <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 1)
one_rep_rob <- function(n, link, rep_id, scen) {
  lk <- c(identity = 0L, logistic = 2L)[[link]]
  beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
  si <- match(scen, c("nonstat", "merr", "skew"))
  set.seed(1e7 * si + 1e5 * lk + 1e3 * n + rep_id)
  e <- prep(sim_rob(n, beta0, link, scen))
  hsel <- select_h(e, lk, beta0)
  hs <- c(hsel["cv"], hsel["cv"], hsel["cvc"], rep(HROB, 2))
  ty <- c("cv", "cv", "cvc", rep("fix", 2 * length(HROB)))
  cr <- c(FALSE, TRUE, TRUE, rep(c(FALSE, TRUE), each = length(HROB)))
  cbind(link = link, n = n, rep = rep_id, scen = scen, fits(e, lk, beta0, unname(hs), ty, cr))
}

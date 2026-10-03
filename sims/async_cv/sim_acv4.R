## Which held-out loss for async CV? Every loss uses the same repeated K-fold
## (K = 5, R = 5), the same training fits per (split, fold, h), pooling over
## folds, averaging over splits, and the same rule (argmin).
## Losses: matched; fixed at h_e = grid min / default lower end / average of
## grid points in [min, lower]; locf; eeq at h_e = lower / grid min.
## Usage: Rscript sim_acv4.R <n> <link> <first_rep> <n_reps> <out.rds>
args <- commandArgs(TRUE)
n <- as.integer(args[1]); link <- args[2]
first <- as.integer(args[3]); nrep <- as.integer(args[4]); out <- args[5]
suppressMessages(library(skmle))
ns <- asNamespace("skmle")
parse_ <- get(".async_parse", ns); subset_ <- get(".async_subset", ns)
fit_ <- get(".async_fit_ti", ns); collapse_ <- get(".async_collapse", ns)
locf_ <- get(".async_locf_pairs", ns); mu_ <- get(".async_mu", ns)
lk <- c(identity = 0L, log = 1L, logistic = 2L)[[link]]
beta0 <- if (link == "logistic") c(0.3, 1.0) else c(0.5, 1.5)
grid <- exp(seq(log(0.005), log(1.5), length.out = 48)); G <- length(grid)
K <- 5L; R <- 5L
kappa <- function(eta) if (lk == 0L) eta^2 else if (lk == 1L) exp(eta) else log1p(exp(-abs(eta))) + pmax(eta, 0)
cterm <- function(eta) if (lk == 0L) 2 * eta else eta   # loss row: om*kappa - c*cterm (+ y^2 for identity)

one <- function(rep_id) {
  set.seed(1e5 * lk + 1e3 * n + rep_id)
  d <- sim_async_data(n = n, beta = beta0, link = link)
  p <- parse_(y ~ x, d$y, d$x, "id", "time", TRUE)
  iqr <- diff(quantile(c(p$yt, p$xt), c(.25, .75), names = FALSE))
  lower <- 2 * iqr * n^(-0.7)
  b2 <- se2 <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    f <- tryCatch(suppressWarnings(kee_async(d$y, d$x, y ~ x, id = id, time = time,
                                             h = grid[g], link = link)),
                  error = function(e) NULL)
    if (is.null(f) || f$convergence == 1L) next
    b <- coef(f); s <- sqrt(diag(vcov(f)))
    if (all(is.finite(b)) && all(is.finite(s))) { b2[g] <- b[2]; se2[g] <- s[2] }
  }
  # Fixed-pair evaluation objects (whole sample; a fold uses its own rows).
  ev <- list(fixed_min = grid[1], fixed_lower = lower, fixed_avg = grid[grid <= lower])
  wts <- lapply(ev, function(hs) {
    acs <- lapply(hs, function(he) collapse_(p, he, FALSE))
    if (any(vapply(acs, is.null, logical(1)))) return(NULL)
    D <- vapply(acs, function(a) sum(a$omega), 0)
    list(om = rowMeans(sweep(vapply(acs, `[[`, numeric(nrow(p$X)), "omega"), 2, D, "/")),
         cc = rowMeans(sweep(vapply(acs, `[[`, numeric(nrow(p$X)), "c"), 2, D, "/")))
  })
  al <- locf_(p)
  wts$locf <- list(om = al$omega / sum(al$omega), cc = al$c / sum(al$omega))
  ee <- list(eeq_lower = collapse_(p, lower, FALSE), eeq_min = collapse_(p, grid[1], FALSE))
  folds <- matrix(replicate(R, sample(rep_len(seq_len(K), p$n))), nrow = p$n)
  names_l <- c("matched", names(wts), names(ee))
  acc <- array(0, c(length(names_l), G, R), dimnames = list(names_l, NULL, NULL))
  den_m <- matrix(0, G, R); bad <- matrix(FALSE, G, R)
  eecount <- matrix(0, G, R)
  for (r in seq_len(R)) for (f in seq_len(K)) {
    te <- which(folds[, r] == f); tr <- subset_(p, which(folds[, r] != f))
    dte <- subset_(p, te); rows <- which(p$xi %in% te); Xr <- p$X[rows, , drop = FALSE]
    for (g in seq_len(G)) {
      ft <- fit_(tr, grid[g], FALSE, lk, 50L, 1e-8)
      if (is.null(ft) || ft$convergence == 1L || !all(is.finite(ft$beta))) { bad[g, r] <- TRUE; next }
      b <- as.numeric(ft$beta)
      acm <- collapse_(dte, grid[g], FALSE)
      if (!is.null(acm)) {
        e <- as.numeric(dte$X %*% b)
        acc["matched", g, r] <- acc["matched", g, r] + sum(acm$omega * kappa(e)) -
          sum(acm$c * cterm(e)) + if (lk == 0L) acm$sy2 else 0
        den_m[g, r] <- den_m[g, r] + sum(acm$omega)
      }
      eta <- as.numeric(Xr %*% b); ka <- kappa(eta); ct <- cterm(eta)
      for (s in names(wts)) if (!is.null(wts[[s]]))
        acc[s, g, r] <- acc[s, g, r] + sum(wts[[s]]$om[rows] * ka - wts[[s]]$cc[rows] * ct)
      mu <- mu_(eta, lk)
      for (s in names(ee)) if (!is.null(ee[[s]])) {
        D <- sum(ee[[s]]$omega[rows]); if (D <= 0) next
        u <- crossprod(Xr, ee[[s]]$c[rows] - ee[[s]]$omega[rows] * mu) / D
        acc[s, g, r] <- acc[s, g, r] + sum(u^2) / K
      }
    }
  }
  acc["matched", , ] <- acc["matched", , ] / den_m
  dim(acc) <- c(length(names_l), G * R); acc[, as.vector(bad)] <- NA
  dim(acc) <- c(length(names_l), G, R); dimnames(acc) <- list(names_l, NULL, NULL)
  curves <- apply(acc, c(1, 2), function(v) if (all(is.na(v))) NA else mean(v, na.rm = TRUE))
  picks <- apply(curves, 1, function(v) if (all(!is.finite(v))) NA_integer_ else which.min(replace(v, !is.finite(v), Inf)))
  avail <- c(wts = !vapply(wts, is.null, logical(1)), ee = !vapply(ee, is.null, logical(1)))
  list(rep = rep_id, lower = lower, b2 = b2, se2 = se2, picks = picks,
       nfail_h = rowSums(bad), avail = avail)
}
res <- lapply(first:(first + nrep - 1L), function(r)
  tryCatch(one(r), error = function(e) list(rep = r, error = conditionMessage(e))))
saveRDS(list(n = n, link = link, grid = grid, beta0 = beta0, res = res), out)
cat("wrote", out, "\n")

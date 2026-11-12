# The held-out score --------------------------------------------------------
#
# The held-out score is an ordinary log-likelihood, with no kernel and no h, so
# it can be compared across bandwidths. The covariate path a held-out subject
# needs between its observation times is filled in by carrying the last
# observation forward (and the first one backward, before the first
# observation).
#
#' Pieces of the held-out log-likelihood that do not depend on the bandwidth
#'
#' @param X_time,obs_times_vec,id_vec,delta Event/censoring times, covariate
#'   observation times, integer subject codes and event indicators, one entry per
#'   row of the model frame.
#' @param knots,tau Interior knots and the end of follow-up, as used by the fit.
#' @param nq Nodes per Legendre panel. The panels break at the observation times,
#'   where the carried-forward path jumps, so the integrand is smooth inside each
#'   one and few nodes are needed.
#' @return A list of the quadrature nodes' spline basis, weights, the row each
#'   node reads its covariate from, the subject each node belongs to, and the same
#'   for the event term.
#' @noRd
locf_score <- function(X_time, obs_times_vec, id_vec, delta, knots, tau,
                       nq = 16) {
  lqr <- gaussquad::legendre.quadrature.rules(nq)[[nq]]
  # Callers code subjects as 1..n, so position m below is subject m and the fold
  # membership vector indexes straight into the result.
  by_subj <- split(seq_along(id_vec), id_vec)
  ns_at <- function(t) {
    as.matrix(splines::ns(t, knots = knots, intercept = TRUE,
                          Boundary.knots = c(0, tau)))
  }

  nodes <- wts <- rows <- owner <- vector("list", length(by_subj))
  ev_row <- integer(length(by_subj))
  ev_X <- numeric(length(by_subj))

  for (m in seq_along(by_subj)) {
    r <- by_subj[[m]]
    r <- r[order(obs_times_vec[r])]
    ot <- obs_times_vec[r]
    Xi <- X_time[r][1]

    # findInterval() gives the last observation at or before a time; 0 means the
    # time precedes every observation, and then the first one is carried back.
    ev_row[m] <- r[max(findInterval(Xi, ot), 1L)]
    ev_X[m] <- Xi

    brk <- unique(c(0, ot[ot > 0 & ot < Xi], Xi))
    if (Xi <= 0 || length(brk) < 2) next
    lo <- brk[-length(brk)]
    hi <- brk[-1]
    half <- (hi - lo) / 2
    nd <- as.numeric(outer(half, lqr$x) + (lo + hi) / 2)
    nodes[[m]] <- nd
    wts[[m]] <- as.numeric(outer(half, lqr$w))
    rows[[m]] <- r[pmax(findInterval(nd, ot), 1L)]
    owner[[m]] <- rep(m, length(nd))
  }

  list(
    node_bs = ns_at(unlist(nodes)),
    node_w = unlist(wts),
    node_row = unlist(rows),
    node_subj = unlist(owner),
    ev_bs = ns_at(ev_X),
    ev_row = ev_row,
    ev_delta = delta[ev_row],
    n_subj = length(by_subj)
  )
}

# The held-out score itself lives in C++, as locf_nll_cpp() in
# src/skmle_cpp.cpp.  It was written here once, and scoring in R meant a second
# implementation of the model's transformation: R/utils.R carried a trans_link()
# that mirrored trans_fun() by hand, floor and all, with nothing able to compare
# them because trans_fun is not exported.  Scoring in C++ means the criterion
# and the objective call the same function, which is a guarantee rather than a
# test someone has to remember to write.

#' Select the Bandwidth by Cross-Validation
#'
#' @description
#' Perform K-fold cross-validation to select the kernel bandwidth for `skmle()`.
#'
#' @inheritParams skmle
#' @param K Number of folds.
#' @param h_grid Optional numeric vector of candidate bandwidth values. If `NULL`,
#'   a grid is generated automatically from the observed time gaps.
#' @param n_h Number of candidate bandwidths to generate when `h_grid` is `NULL`.
#' @param seed Optional integer seed for the random subject-to-fold assignment.
#'   If `NULL`, the current RNG state is used and no explicit seed is set.
#'   A supplied seed is used locally: the global RNG state is restored on exit.
#' @param quiet Logical; if `TRUE`, suppress progress output.
#'
#' @section Kernel choice:
#' `one_sided` is used inside the fold loop as well as being passed through to
#' the refit, so the bandwidth is selected under the same kernel the final fit
#' uses.
#'
#' @details
#' # How the folds are formed
#'
#' `skmle_cv()` splits subjects across folds. Several rows belong to the same
#' subject in long format, so splitting by row would put one subject on both
#' sides of the split.
#'
#' After choosing the bandwidth with the smallest average validation loss, the
#' function refits `skmle()` on the full data set using the selected value.
#' Because the fold assignment is random, pass `seed` (or `set.seed()` before
#' calling) to make the grid selection reproducible.
#'
#' # The loss
#'
#' Each training fold is fitted at the candidate bandwidth, and the fit is scored
#' on the held-out subjects by the log-likelihood per subject,
#'
#' \deqn{-\frac{1}{n_{\mathrm{test}}} \sum_{i \in \mathrm{test}}
#'   \left[ \delta_i \log g\{\hat\alpha(X_i) + Z_i(X_i)^\top \hat\beta\}
#'   - \int_0^{X_i} g\{\hat\alpha(t) + Z_i(t)^\top \hat\beta\}\,dt \right],}
#'
#' where \eqn{Z_i(t)} is the covariate carried forward from the last observation
#' at or before \eqn{t}. This loss does not depend on the bandwidth itself, so
#' it compares candidates on an equal footing.
#'
#' # The grid
#'
#' When `h_grid` is `NULL` the grid is log-spaced over
#' \eqn{[\max\{\min_i (X_i - T_{ij})_+,\ \tau n^{-0.6}\},\
#' \min\{\max_i \max_j (X_i - T_{ij})_+,\ \tau n^{-0.3}\}]}, with
#' \eqn{\tau = \max_i X_i}, so it adapts to the scale of the times on its own.
#' If the selected value is an endpoint of the grid, `skmle_cv()` warns; supply
#' a wider `h_grid`.
#'
#' `cv_results` lists the loss and its standard error across folds for every
#' candidate.
#'
#' @return
#' An object of class `cv.skmle` with components:
#'
#' * `h_cv`: selected bandwidth,
#' * `fit`: `skmle` fit refit on the full data,
#' * `cv_results`: data frame of candidate bandwidths, their CV losses, and the
#'   standard error of each loss across the folds,
#' * `h_grid`: bandwidth grid used in the search,
#' * `fold_id`: the subject-to-fold assignment vector (length `n`),
#' * `seed`: the value of `seed` supplied by the user, or `NULL`,
#' * `call`: the matched call.
#'
#' @examples
#' \donttest{
#' library(survival)
#'
#' set.seed(123)
#' dat <- sim_skmle_data(
#'   n = 60,
#'   mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
#'   mu_bar = 8,
#'   alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
#'   beta = c(1, -0.5),
#'   s = 0,
#'   cen = 0.7
#' )
#'
#' cv_fit <- skmle_cv(
#'   Surv(X, delta) ~ covariates,
#'   data = dat,
#'   id = id,
#'   obs_times = obs_times,
#'   s = 0,
#'   K = 3,
#'   seed = 2026,
#'   quiet = TRUE
#' )
#'
#' cv_fit$h_cv
#' # Read the whole table, not just the selection: a minimum on the edge of the
#' # grid is a boundary artefact and warns.
#' cv_fit$cv_results
#' summary(cv_fit$fit)
#' }
#'
#' @export
skmle_cv <- function(formula, data, id, obs_times, s = 0, K = 5, h_grid = NULL,
                     n_h = 10, nknots = 3, lq_nodes = 64,
                     maxeval = 10000, xtol_rel = 1e-6, seed = NULL,
                     quiet = FALSE, one_sided = TRUE, weight = NULL) {
  weight <- resolve_weight(weight, one_sided)
  if (missing(formula) || missing(data) || missing(id) || missing(obs_times)) {
    stop("formula, data, id and obs_times must all be supplied")
  }

  # parse formula and extract data once
  call <- match.call()
  m <- match(c("formula", "data", "id", "obs_times"), names(call), 0L)
  mf_call <- call[c(1L, m)]
  mf_call[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_call, parent.frame())

  Y <- stats::model.response(mf)
  if (!inherits(Y, "Surv")) {
    stop("Response must be a survival object created with Surv()")
  }

  X_time <- Y[, 1]
  delta <- Y[, 2]
  if (anyNA(X_time) || anyNA(delta)) {
    stop("missing values in survival response not permitted")
  }

  Z <- stats::model.matrix(formula, data = mf)
  if (ncol(Z) > 0 && colnames(Z)[1] == "(Intercept)") {
    Z <- Z[, -1, drop = FALSE]
  }
  if (ncol(Z) == 0) stop("model must contain at least one covariate")

  id_raw <- stats::model.extract(mf, "id")
  obs_times_vec <- as.numeric(stats::model.extract(mf, "obs_times"))

  if (length(id_raw) != length(X_time) || length(obs_times_vec) != length(X_time)) {
    stop("Length of 'id' and 'obs_times' must match number of rows in data/formula")
  }

  check_time_scale(X_time, obs_times_vec)

  id_vec <- as.integer(factor(id_raw))
  unique_ids <- unique(id_vec)
  n <- length(unique_ids)

  # K = 1 leaves the training fold empty, which reaches the optimiser as a fit to
  # no data at all.
  if (!is.numeric(K) || length(K) != 1L || is.na(K) || K < 2) {
    stop("'K' must be at least 2", call. = FALSE)
  }

  if (K > n) {
    if (!quiet) message("Requested K = ", K, " exceeds n = ", n,
                        " subjects; using K = ", n, " (leave-one-out).")
    K <- n
  }

  # End of follow-up. Needed by the default grid as well as by the basis and the
  # quadrature below, so it is computed before either.
  tau <- max(X_time)

  # generate h_grid if not provided
  if (is.null(h_grid)) {
    pos_diffs <- X_time - obs_times_vec
    pos_diffs <- pos_diffs[pos_diffs > 0]

    if (length(pos_diffs) == 0) {
      stop("No observation times are strictly prior to failure/censoring times. Cannot determine default h_grid.")
    }

    hmin <- max(min(pos_diffs), tau * n^(-0.6))

    # max diff per subject
    max_diffs <- tapply(X_time - obs_times_vec, id_vec, function(x) {
      v <- x[x > 0]
      if (length(v) > 0) max(v) else NA
    })
    hmax <- min(max(max_diffs, na.rm = TRUE), tau * n^(-0.3))

    if (hmin >= hmax) {
      h_grid <- c(hmin)
    } else {
      h_grid <- exp(seq(log(hmin), log(hmax), length.out = n_h + 1)[-1])
    }
  }

  # split data into K folds by subject id (reproducible if `seed` given)
  # without changing the caller's RNG stream
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  fold_id_subj <- sample(rep(1:K, length.out = n))

  # pre-compute quadrature and knots
  knots <- tau * (1:nknots) / (nknots + 1)
  lqrule <- gaussquad::legendre.quadrature.rules(lq_nodes)[[lq_nodes]]
  lq_x <- lqrule$x
  lq_w <- lqrule$w
  tts <- 0.5 * tau * (lq_x + 1)
  bsmat <- as.matrix(splines::ns(X_time, knots = knots, intercept = TRUE, Boundary.knots = c(0, tau)))
  bsmat_tt_mat <- as.matrix(splines::ns(tts, knots = knots, intercept = TRUE, Boundary.knots = c(0, tau)))
  p <- ncol(Z)
  gammap <- ncol(bsmat_tt_mat)

  # The held-out score is a function of the data alone, so it is built once and
  # reused for every (h, fold) pair. See locf_score(). It stays in R because it
  # is set-up rather than a loop, it needs splines::ns() and findInterval(), and
  # it evaluates no part of the model: it returns a basis, weights and row
  # indices. The SCORING runs in C++, because that evaluates the transformation,
  # and a second implementation of the transformation is the thing worth
  # avoiding.
  score <- locf_score(X_time, obs_times_vec, id_vec, delta, knots, tau)

  cv_losses <- skmle_cv_cpp(
    p = p, gammap = gammap, s = as.numeric(s), tau = as.numeric(tau),
    h_grid = as.numeric(h_grid), K = as.integer(K),
    fold_id_subj = as.integer(fold_id_subj),
    id_vec = as.integer(id_vec) - 1L,
    covariates = as.matrix(Z), bsmat = bsmat,
    X = as.numeric(X_time), obs_times = as.numeric(obs_times_vec),
    delta = as.numeric(delta),
    lq_x = as.numeric(lq_x), lq_w = as.numeric(lq_w),
    bsmat_tt_all = bsmat_tt_mat, tts = as.numeric(tts),
    node_bs = score$node_bs, node_w = as.numeric(score$node_w),
    node_row = as.integer(score$node_row) - 1L,
    node_subj = as.integer(score$node_subj) - 1L,
    ev_bs = score$ev_bs, ev_row = as.integer(score$ev_row) - 1L,
    ev_delta = as.numeric(score$ev_delta), n_subj = as.integer(score$n_subj),
    maxeval = as.integer(maxeval), xtol_rel = as.numeric(xtol_rel),
    quiet = as.logical(quiet),
    weight = weight, support = weight_support(weight),
    one_sided = as.logical(one_sided)
  )

  cv_results <- tibble::tibble(
    h = h_grid,
    cvloss = as.numeric(cv_losses[1, ]),
    # standard error of each loss across the folds
    se = as.numeric(cv_losses[2, ])
  )
  # A candidate scores +Inf when a training fit failed or a held-out score
  # overflowed on any fold. which.min() would return the first of an all-Inf
  # grid as though it had won, so refuse instead and say why.
  ok <- is.finite(cv_results$cvloss)
  if (!any(ok)) {
    stop(
      "no bandwidth in 'h_grid' gave a finite cross-validation loss on every ",
      "fold: each candidate had a training fit fail or a held-out score ",
      "overflow. Supply a grid of larger bandwidths, or fewer folds.",
      call. = FALSE
    )
  }
  if (!all(ok) && !quiet) {
    message(sum(!ok), " of ", length(ok), " candidate bandwidths had a fold ",
            "fail or overflow and were not eligible.")
  }
  best_h <- cv_results$h[ok][which.min(cv_results$cvloss[ok])]

  if (length(h_grid) > 1 && best_h %in% range(h_grid)) {
    warning(
      "the selected bandwidth is an endpoint of 'h_grid', so it is the best of ",
      "the values offered rather than a minimum. Widen the grid.",
      call. = FALSE
    )
  }

  if (!quiet) cat(sprintf("Selected optimal h = %f\n", best_h))

  # Refit with full data. Freeze `data` into the call so that lazy
  # re-evaluation of the `data` expression cannot silently produce a
  # different data frame from the one used for CV.
  fit_call <- call
  fit_call[[1]] <- quote(skmle::skmle)
  fit_call$data <- data
  fit_call$h <- best_h
  fit_call$K <- NULL
  fit_call$h_grid <- NULL
  fit_call$n_h <- NULL
  fit_call$seed <- NULL
  fit_call$quiet <- NULL

  fit <- eval(fit_call, envir = parent.frame())

  # The frozen data frame has done its job by now.  Leaving it in the fitted
  # object's call means print() and summary() deparse the whole data set in
  # place of the call, so put the user's original expression back.  The fit
  # itself was still evaluated against the frozen frame.
  fit$call <- fit_call
  fit$call$data <- call$data

  out <- list(
    h_cv = best_h,
    fit = fit,
    cv_results = cv_results,
    h_grid = h_grid,
    fold_id = fold_id_subj,
    seed = seed,
    call = call
  )

  class(out) <- "cv.skmle"
  return(out)
}


#' @param x A `cv.skmle` object.
#' @param ... Ignored.
#' @rdname skmle_cv
#' @export
print.cv.skmle <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n", length(unique(x$fold_id)),
      "-fold subject-level cross-validation\n\n", sep = "")
  print(x$cv_results, row.names = FALSE, digits = 5)
  cat("\nSelected h = ", format(x$h_cv), "\n", sep = "")
  cat("\nCoefficients at the refit:\n")
  print(x$fit$coefficients)
  invisible(x)
}

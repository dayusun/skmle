# The held-out score --------------------------------------------------------
#
# The kernel-weighted log-likelihood that `skmle()` maximises cannot be compared
# across bandwidths. Its weights are W(u/h)/h, so the weight a subject
# contributes falls away as h grows, and both the event term and the
# cumulative-hazard term shrink with it whatever the fit is worth. The criterion
# then decreases monotonically in h and the largest candidate always wins --
# normalising by the admitted weight does not rescue it either, because the fits
# obtained at large h keep scoring better even against a bandwidth held fixed.
#
# So the held-out score is an ordinary log-likelihood: no kernel, no h. The
# covariate path a held-out subject needs between its observation times is filled
# in by carrying the last observation forward (and the first one backward, before
# the first observation). That approximation is the same for every candidate, so
# it cancels out of the comparison and what is left is the quality of the fit.
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

#' Held-out negative log-likelihood for one fold
#'
#' @param score Output of `locf_score()`.
#' @param Z Covariate matrix, all rows.
#' @param beta,gamma Coefficients from the training fit.
#' @param s Transformation parameter.
#' @param keep Integer subject codes held out in this fold.
#' @return The negative log-likelihood per held-out subject.
#' @noRd
locf_nll <- function(score, Z, beta, gamma, s, keep) {
  haz <- trans_link(
    as.numeric(score$node_bs %*% gamma) +
      as.numeric(Z[score$node_row, , drop = FALSE] %*% beta), s
  )
  cumhaz <- numeric(score$n_subj)
  agg <- rowsum(score$node_w * haz, score$node_subj, reorder = TRUE)
  cumhaz[as.integer(rownames(agg))] <- agg[, 1]

  ev <- trans_link(
    as.numeric(score$ev_bs %*% gamma) +
      as.numeric(Z[score$ev_row, , drop = FALSE] %*% beta), s
  )
  -sum(score$ev_delta[keep] * log(ev[keep]) - cumhaz[keep]) / length(keep)
}

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
#' # What the held-out loss is
#'
#' Each training fold is fitted at the candidate bandwidth, and the fit is scored
#' on the held-out subjects by an ordinary log-likelihood per subject,
#'
#' \deqn{-\frac{1}{n_{\mathrm{test}}} \sum_{i \in \mathrm{test}}
#'   \left[ \delta_i \log g\{\hat\alpha(X_i) + Z_i(X_i)^\top \hat\beta\}
#'   - \int_0^{X_i} g\{\hat\alpha(t) + Z_i(t)^\top \hat\beta\}\,dt \right],}
#'
#' where \eqn{Z_i(t)} is the covariate carried forward from the last observation
#' at or before \eqn{t} (and the first observation carried back, before the first
#' observation time). The integral is exact up to the Legendre rule applied
#' between consecutive observation times, where the carried-forward path jumps.
#'
#' No kernel and no bandwidth appear in that expression, and that is the point.
#' The kernel-weighted log-likelihood `skmle()` maximises **cannot** be compared
#' across bandwidths: its weights are \eqn{W(u/h)/h}, so the weight each subject
#' contributes falls away as `h` grows, and the criterion decreases
#' monotonically in `h` whatever the fit is worth. Scored that way the largest
#' candidate wins every grid on every data set, and coefficients further from the
#' truth are preferred to coefficients nearer it. Dividing by the admitted
#' weight, or scoring at one bandwidth held fixed across the grid, does not
#' rescue it. The carried-forward likelihood is one yardstick for every
#' candidate, so its approximation cancels out of the comparison and what is left
#' is the quality of the fit.
#'
#' # The default grid
#'
#' When `h_grid` is `NULL` the grid is log-spaced over
#' \eqn{[\max\{\min_i (X_i - T_{ij})_+,\ \tau n^{-0.6}\},\
#' \min\{\max_i \max_j (X_i - T_{ij})_+,\ \tau n^{-0.3}\}]}, with
#' \eqn{\tau = \max_i X_i}, so it adapts to the scale of the times on its own.
#'
#' Always look at `cv_results`. A minimum at an endpoint of the grid raises a
#' warning: the selected value is then the best of the values offered rather than
#' a minimum, and the grid should be widened. On the automatic grid that warning
#' is common, because \eqn{n^{-0.3}} is the rate the asymptotics assume while the
#' finite-sample minimum of the loss frequently lies above it. Widening `h_grid`
#' by hand shows where the curve turns.
#'
#' # How sharp the selection is
#'
#' Not very, and the `se` column says so: it is the standard error of each loss
#' across the folds, and over a wide middle range of `h` the losses sit inside one
#' standard error of each other. Read the curve, not only `h_cv`.
#'
#' The criterion scores prediction of the held-out hazard, in which the baseline
#' \eqn{\hat\alpha} can absorb attenuation in \eqn{\hat\beta}, so it leans towards
#' more smoothing than the coefficients on their own would want. Over 10
#' replicates at \eqn{n = 200} on a grid spanning `0.05` to `0.9`, the mean
#' squared error of \eqn{\hat\beta} at the selected bandwidth was `0.155`, against
#' `0.205` at the largest candidate and `0.065` at the bandwidth an oracle would
#' have picked. Resist the temptation to correct the lean by taking the smallest
#' bandwidth within one standard error of the minimum: that lands in the noisy
#' small-`h` end, and scored `0.205` over the same replicates -- no better than
#' taking the largest candidate. The rise at the left of the curve is real.
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
                     quiet = FALSE, one_sided = TRUE) {
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
  if (!is.null(seed)) set.seed(seed)
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
  # reused for every (h, fold) pair. See locf_score().
  score <- locf_score(X_time, obs_times_vec, id_vec, delta, knots, tau)
  Zbs <- cbind(Z, bsmat)

  cv_losses <- vapply(h_grid, function(h) {
    if (!quiet) cat(sprintf("Evaluating bandwidth h = %g\n", h))
    kerval <- kernel_weights(X_time - obs_times_vec, h, one_sided)
    kerval_tt <- kernel_weights(outer(tts, obs_times_vec, "-"), h, one_sided)
    # Same support the fit uses, so the constraint set matches the objective.
    feasible <- if (s == 0) NULL else {
      (abs(X_time - obs_times_vec) <= h) & (!one_sided | (X_time > obs_times_vec))
    }

    fold_losses <- vapply(seq_len(K), function(k) {
      tr <- fold_id_subj[id_vec] != k
      ineqmat <- matrix(numeric(0), nrow = 0, ncol = p + gammap)
      if (s != 0 && any(tr & feasible)) {
        ineqmat <- Zbs[tr & feasible, , drop = FALSE]
      }
      fit <- skmle_cpp_fit(
        n = length(unique(id_vec[tr])), p = p, gammap = gammap,
        s = as.numeric(s), h = h, tau = tau,
        covariates = Z[tr, , drop = FALSE],
        bsmat = bsmat[tr, , drop = FALSE],
        X = X_time[tr], obs_times = obs_times_vec[tr], delta = delta[tr],
        kerval = kerval[tr],
        lq_x = lq_x, lq_w = lq_w,
        bsmat_tt_all = bsmat_tt_mat,
        kerval_tt_all = kerval_tt[, tr, drop = FALSE],
        ineqmat = ineqmat,
        maxeval = as.integer(maxeval), xtol_rel = as.numeric(xtol_rel)
      )
      # A training fit that failed cannot produce a trustworthy held-out score:
      # let this bandwidth lose the grid outright.
      if (fit$status < 0) {
        return(Inf)
      }
      locf_nll(score, Z, fit$solution[seq_len(p)], fit$solution[-seq_len(p)],
               s, which(fold_id_subj == k))
    }, numeric(1))

    c(mean(fold_losses), stats::sd(fold_losses) / sqrt(K))
  }, numeric(2))

  cv_results <- tibble::tibble(
    h = h_grid,
    cvloss = as.numeric(cv_losses[1, ]),
    # The criterion is flat over a wide range of h, so the standard error across
    # folds is worth reporting: it says how much of the gap between two
    # candidates is signal. It is not an invitation to a one-standard-error rule
    # -- that rule was measured here and made beta worse, see ?skmle_cv.
    se = as.numeric(cv_losses[2, ])
  )
  best_h <- cv_results$h[which.min(cv_results$cvloss)]

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

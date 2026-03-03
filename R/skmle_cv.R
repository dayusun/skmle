#' Bandwidth Selection for skmle Model using Cross-Validation
#'
#' @description Performs K-fold cross-validation to select the optimal bandwidth `h` for the `skmle` model.
#'
#' @inheritParams skmle
#' @param K Number of folds for cross-validation. Default is 5.
#' @param h_grid Optional numeric vector of bandwidth candidates. If NULL, a grid is generated automatically based on observation spacings.
#' @param n_h Number of bandwidths to generate if `h_grid` is NULL. Default is 10.
#' @param quiet Whether to suppress progress messages. Default is FALSE.
#'
#' @return A list of class `cv.skmle` containing the optimal bandwidth `h`, the estimated model with optimal bandwidth (`fit`), a data.frame of `cv_results`, and `h_grid`.
#'
#' @export
skmle_cv <- function(formula, data, id, obs_times, s, K = 5, h_grid = NULL, n_h = 10, nknots = 3, norder = 3, lq_nodes = 64, maxeval = 10000, xtol_rel = 1e-6, quiet = FALSE) {
  if (missing(formula) || missing(data) || missing(id) || missing(obs_times) || missing(s)) {
    stop("formula, data, id, obs_times, and s must all be supplied")
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
  
  id_vec <- as.numeric(stats::model.extract(mf, "id"))
  obs_times_vec <- as.numeric(stats::model.extract(mf, "obs_times"))
  
  if (length(id_vec) != length(X_time) || length(obs_times_vec) != length(X_time)) {
    stop("Length of 'id' and 'obs_times' must match number of rows in data/formula")
  }
  
  unique_ids <- unique(id_vec)
  n <- length(unique_ids)
  
  if (K > n) K <- n
  
  # generate h_grid if not provided
  if (is.null(h_grid)) {
    pos_diffs <- X_time - obs_times_vec
    pos_diffs <- pos_diffs[pos_diffs > 0]
    
    if (length(pos_diffs) == 0) {
      stop("No observation times are strictly prior to failure/censoring times. Cannot determine default h_grid.")
    }
    
    hmin <- max(min(pos_diffs), n^(-0.6))
    
    # max diff per subject
    max_diffs <- tapply(X_time - obs_times_vec, id_vec, function(x) {
      v <- x[x > 0]
      if (length(v) > 0) max(v) else NA
    })
    hmax <- min(max(max_diffs, na.rm = TRUE), n^(-0.3))
    
    if (hmin >= hmax) {
      h_grid <- c(hmin)
    } else {
      h_grid <- exp(seq(log(hmin), log(hmax), length.out = n_h + 1)[-1])
    }
  }
  
  # split data into K folds by subject id
  # Generate fold mapping
  fold_id_subj <- sample(rep(1:K, length.out = n))
  
  # pre-compute quadrature and knots
  knots <- (1:nknots) / (nknots + 1)
  lqrule <- gaussquad::legendre.quadrature.rules(lq_nodes)[[lq_nodes]]
  lq_x <- lqrule$x
  lq_w <- lqrule$w
  tts <- 0.5 * lq_x + 0.5
  bsmat_tt_mat <- as.matrix(splines::ns(tts, knots = knots, intercept = TRUE, Boundary.knots = c(0, 1)))
  
  # C++ implementation handles the iteration
  cv_losses <- skmle_cv_cpp(
    n = n, p = ncol(Z), gammap = ncol(bsmat_tt_mat), 
    s = as.numeric(s), 
    h_grid = as.numeric(h_grid),
    K = as.integer(K), 
    fold_id = as.numeric(fold_id_subj), 
    id_vec = as.numeric(id_vec),
    covariates = as.matrix(Z), 
    bsmat = as.matrix(splines::ns(X_time, knots = knots, intercept = TRUE, Boundary.knots = c(0, 1))),
    X = as.numeric(X_time), 
    obs_times = as.numeric(obs_times_vec), 
    delta = as.numeric(delta),
    lq_x = as.numeric(lq_x), 
    lq_w = as.numeric(lq_w), 
    bsmat_tt_all = as.matrix(bsmat_tt_mat),
    maxeval = as.integer(maxeval), 
    xtol_rel = as.numeric(xtol_rel), 
    quiet = as.logical(quiet)
  )
  
  cv_results <- data.frame(h = h_grid, cvloss = as.numeric(cv_losses))

  
  best_h <- cv_results$h[which.min(cv_results$cvloss)]
  
  if (!quiet) cat(sprintf("Selected optimal h = %f\n", best_h))
  
  # Refit with full data
  fit_call <- call
  fit_call[[1]] <- quote(skmle::skmle)
  fit_call$h <- best_h
  fit_call$K <- NULL
  fit_call$h_grid <- NULL
  fit_call$n_h <- NULL
  fit_call$quiet <- NULL
  
  fit <- eval(fit_call, envir = parent.frame())
  
  out <- list(
    h_cv = best_h,
    fit = fit,
    cv_results = cv_results,
    h_grid = h_grid
  )
  
  class(out) <- "cv.skmle"
  return(out)
}

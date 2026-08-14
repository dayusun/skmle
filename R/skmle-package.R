#' @details
#' A longitudinal covariate is rarely measured at the time you need it. `skmle`
#' handles that mismatch directly, by weighting each observation according to
#' how far its measurement time sits from the time being modelled, rather than
#' carrying a value forward or smoothing the covariate and substituting it.
#' Two outcome types are covered.
#'
#' @section Survival outcomes:
#' [skmle()] fits the transformed hazards family by the Sieve Maximum
#' Kernel-weighted Log-likelihood Estimator of Sun, Sun, Zhao and Cao (2025).
#' The Box-Cox parameter `s` indexes the family: `s = 0` is proportional
#' hazards, `s = 1` is additive hazards, and other values interpolate.
#' [kee_cox()] and [kee_additive()] are faster specialised estimating equations
#' for the two named cases, and [skmle_cv()] selects the bandwidth by
#' subject-level cross-validation. Simulate with [sim_skmle_data()].
#'
#' @section Asynchronous longitudinal outcomes:
#' When the outcome is itself a sparsely observed longitudinal process, recorded
#' on a time grid that does not line up with the covariate's, [kee_async()] and
#' [kee_async_td()] fit generalised linear models by the kernel-weighted
#' estimating equations of Cao, Zeng and Fine (2015), with time-invariant
#' coefficients and with a coefficient curve \eqn{\beta(t)} respectively.
#' Simulate with [sim_async_data()].
#'
#' @section Half and full kernels:
#' Every estimator takes `one_sided`. A half kernel admits only measurements
#' strictly before the time being modelled, which is the risk-set restriction of
#' a hazard model and the causal reading of a covariate path; a full kernel
#' smooths from both sides. The survival estimators default to the half kernel,
#' as published; the asynchronous ones default to the full kernel, as published.
#'
#' @references
#' Sun, D., Sun, Z., Zhao, X. and Cao, H. (2025). Kernel meets sieve:
#' transformed hazards models with sparse longitudinal covariates.
#' \emph{Journal of the American Statistical Association} 120, 2580-2591.
#'
#' Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
#' asynchronous longitudinal data. \emph{Journal of the Royal Statistical
#' Society, Series B} 77, 755-776.
#'
#' @useDynLib skmle, .registration = TRUE
#' @import nloptr
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.extract rnorm
#' @importFrom utils globalVariables
#' @keywords internal
"_PACKAGE"

globalVariables(c("Time", "Baseline"))

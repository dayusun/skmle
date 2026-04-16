#' @useDynLib skmle, .registration = TRUE
#' @import nloptr
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.extract rnorm
NULL

utils::globalVariables(c("Time", "Baseline"))

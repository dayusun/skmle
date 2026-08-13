#' @useDynLib skmle, .registration = TRUE
#' @import nloptr
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.extract rnorm
#' @importFrom utils globalVariables
NULL

globalVariables(c("Time", "Baseline"))

.onLoad <- function(libname, pkgname) {
  # Ensure nloptr is loaded before compiled code looks up its C-callables.
  requireNamespace("nloptr", quietly = TRUE)
}

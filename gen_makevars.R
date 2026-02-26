writeLines(
  c(
    "PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)",
    paste0("PKG_CXXFLAGS = -I\"", system.file("include", package="nloptr"), "\"")
  ),
  con = "src/Makevars.win"
)

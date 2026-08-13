## Submission

This is the initial CRAN submission of `skmle`.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 -- 0 errors | 0 warnings | 0 notes
* macOS builder (mac.r-project.org), macOS Tahoe 26.6, aarch64-apple-darwin23,
  R 4.6.1 Patched -- Status: OK (errors: no, warnings: no, notes: no)
* win-builder R-devel and R-release -- submitted; TODO record results here
  before submitting to CRAN
* GitHub Actions: R-devel, R-release, R-oldrel on Linux / macOS / Windows

## R CMD check results

0 errors | 0 warnings | 1 note

* Note: new submission.

The benchmark vignette uses `SurvSparse`, `bench`, `dplyr` and `ggplot2`
from Suggests; all of its chunks are gated on `requireNamespace()` so the
vignette builds when those packages are unavailable.

## Method reference

The methodology is described in

Sun, D., Sun, Z., Zhao, X., and Cao, H. (2025).
"Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal
Covariates." *Journal of the American Statistical Association*.
doi:10.1080/01621459.2025.2476781

and is cited in the `Description` field.

## Compiled code

The package links to `Rcpp`, `RcppArmadillo`, and `nloptr` (via `LinkingTo`)
and uses the `nloptr` C API for optimization. `SystemRequirements: C++17`.

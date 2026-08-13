## Submission

This is the initial CRAN submission of `skmle`.

## Test environments

* local: Ubuntu 24.04, R 4.6.1
* GitHub Actions: R-devel, R-release, R-oldrel on Linux / macOS / Windows
* win-builder: TODO before submitting -- run devtools::check_win_devel() and
  devtools::check_win_release(), then record the results here.

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

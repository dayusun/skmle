## Submission

This is the initial CRAN submission of `skmle`.

## Test environments

* local: Ubuntu 24.04, R 4.5.3
* GitHub Actions (planned): R-devel, R-release, R-oldrel on
  Linux / macOS / Windows
* R-hub and win-builder (planned)

## R CMD check results

0 errors | 0 warnings | 1 note

* Note: new submission.

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

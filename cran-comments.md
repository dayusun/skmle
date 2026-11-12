## Submission

This is the initial CRAN submission of `skmle`.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 (R CMD check --as-cran)
* GitHub Actions: Windows (R release), macOS (R release), Ubuntu (R devel,
  release and oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* Note: new submission.

The benchmark vignette uses `SurvSparse`, `bench`, `dplyr` and `ggplot2`
from Suggests; all of its chunks are gated on `requireNamespace()` so the
vignette builds when those packages are unavailable.

## Method references

The package implements two published methods, both cited in the `Description`
field.

Survival outcomes:

Sun, D., Sun, Z., Zhao, X., and Cao, H. (2025).
"Kernel Meets Sieve: Transformed Hazards Models with Sparse Longitudinal
Covariates." *Journal of the American Statistical Association*.
doi:10.1080/01621459.2025.2476781

Asynchronous longitudinal outcomes:

Cao, H., Zeng, D., and Fine, J. P. (2015).
"Regression Analysis of Sparse Asynchronous Longitudinal Data."
*Journal of the Royal Statistical Society: Series B* 77, 755-776.
doi:10.1111/rssb.12086

## Compiled code

The package links to `Rcpp`, `RcppArmadillo`, and `nloptr` (via `LinkingTo`)
and uses the `nloptr` C API for optimization. `SystemRequirements: C++17`.

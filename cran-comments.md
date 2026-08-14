## Submission

This is the initial CRAN submission of `skmle`.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 -- 0 errors | 0 warnings | 0 notes
* macOS builder (mac.r-project.org), macOS Tahoe 26.6, aarch64-apple-darwin23,
  R 4.6.1 Patched -- Status: OK (errors: no, warnings: no, notes: no)
* GitHub Actions -- all Status: OK
  - Windows Server 2022 x64, x86_64-w64-mingw32, R 4.6.1 ucrt
  - macOS, aarch64-apple-darwin23, R 4.6.1
  - Ubuntu, x86_64-pc-linux-gnu, R Under development (2026-06-21 r90185)
  - Ubuntu, x86_64-pc-linux-gnu, R 4.6.1
  - Ubuntu, x86_64-pc-linux-gnu, R 4.5.3 (oldrel-1)
* win-builder, Windows Server 2022 x64, x86_64-w64-mingw32, gcc 14.3.0,
  under R-devel (r90394 ucrt) and R 4.6.1 ucrt -- Windows binaries built with
  no example or test failures. win-builder truncates its published
  00check.log at the CRAN incoming feasibility step, so no Status line is
  quoted for those two runs; the Windows GitHub Actions job above covers
  the same platform and reports Status: OK.

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

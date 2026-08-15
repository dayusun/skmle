# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package overview

`skmle` is an R package implementing the Sieve Maximum Kernel-weighted
Log-likelihood Estimator (SMKLE) from Sun et al. 2025 (JASA) for
transformed hazards models with sparse, intermittently observed
longitudinal covariates. The heavy lifting is in a C++ backend linked
via `Rcpp` + `RcppArmadillo`, with `nloptr`’s C API used as the
optimizer.

## Build / develop commands

Use `devtools` from the package root:

- `devtools::load_all()` — compile C++ and load the package for
  interactive work.
- `devtools::document()` — regenerate `man/*.Rd` and `NAMESPACE` from
  roxygen comments. Always run after editing roxygen blocks.
- [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
  — regenerate `src/RcppExports.cpp` and `R/RcppExports.R` after
  adding/changing `// [[Rcpp::export]]` signatures. Must be run before
  `document()` when exported C++ prototypes change.
- `devtools::install()` / `R CMD INSTALL .` — full install.
- `R CMD build . && R CMD check skmle_*.tar.gz` — full CRAN-style check.

### Makevars nloptr include path

`src/Makevars` resolves the `nloptr` include dir at build time via
`Rscript -e "cat(system.file('include', package='nloptr'))"`.
`src/Makevars.win` hardcodes a Windows path — regenerate it on a Windows
machine with `source("dev/gen_makevars.R")` if the local nloptr install
path changes. Never edit `src/Makevars.win` by hand when changing
versions.

## Running verification scripts

There is no `testthat` suite yet. Correctness is verified by hand via
scripts under `tests/` that compare the packaged implementation against
a pure-R prototype.

- Each script expects to be run from the project root (paths like
  `dev/Rprototype_funcs.R` are relative).
- Scripts source `dev/Rprototype_funcs.R` (ground-truth R prototypes)
  and `tests/multiple_Cao_*_funcs.R` (per-model references), then run
  both side-by-side and print `all.equal` comparisons.
- Typical flow, from the package root in R:
  `source("tests/verify_skmle.R")`, `source("tests/verify_kee_cox.R")`,
  `source("tests/verify_kee_additive.R")`,
  `source("tests/verify_skmle_cv.R")`, `source("tests/test_plot.R")`.
- `tests/verify_0.R` re-uses the fitted object from `verify_skmle.R` to
  sanity-check the NLL at `maxeval = 1`.
- `dev/` is in `.Rbuildignore` and `.gitignore` — the prototype files
  live outside the built package but are required for the verify scripts
  to run.

## Architecture

### Three estimators, one shared pipeline

All three estimator wrappers
([`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md),
[`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md),
[`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
in `R/skmle.R`, `R/kee.R`) share the same front-end pattern:

1.  Parse the formula via
    [`stats::model.frame`](https://rdrr.io/r/stats/model.frame.html)
    using a manual
    [`match.call()`](https://rdrr.io/r/base/match.call.html) + `mf_call`
    trick so `id` and `obs_times` are carried alongside the `Surv`
    response and NA-synchronised with it.
2.  Drop the `(Intercept)` column from `model.matrix`. Intercept-only
    models are explicitly rejected.
3.  Coerce `id` to integer codes via `as.integer(factor(id_vec))` —
    non-numeric identifiers are supported because of this.
4.  Compute the Epanechnikov kernel `kerval = (1 - u^2) * 0.75 / h` with
    `u = (X - obs_times)/h`, zeroed where `X <= obs_times` (only past
    observations inform each subject’s likelihood contribution).
5.  Delegate to a `// [[Rcpp::export]]` function in `src/` via the thin
    `R/RcppExports.R` bridge. Never call `.Call` directly.
6.  Post-process: assemble a sandwich variance from the `A`/`B` (or
    `W`/`Sigma`) matrices returned by C++, then wrap in an S3 class
    (`skmle`, `kee`, or `cv.skmle`).

### skmle sieve fit

[`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md)
precomputes a B-spline sieve
([`splines::ns`](https://rdrr.io/r/splines/ns.html) with equally spaced
interior knots `(1:nknots)/(nknots+1)` on `[0,1]`) for the baseline
hazard, plus Legendre-Gauss quadrature nodes from `gaussquad` for the
cumulative-hazard integral. Two basis matrices are built: `bsmat` at
event/censoring times `X`, and `bsmat_tt_all` at mapped quadrature nodes
`0.5*lq_x + 0.5`, along with `kerval_tt_all` = kernel values between
each quadrature node and every observation time.

When `s != 0` (non-Cox transformations), an `ineqmat` of nloptr
inequality constraints is built from rows where
`0 < X - obs_times <= h`; this enforces `s*(Z'β + B'γ) + 1 > 0` inside
the Box-Cox transform. `src/skmle_cpp.cpp` defines `trans_fun`,
`trans_fun_d`, etc. as the Box-Cox family with `s = 0` special-cased to
`exp`.

`calc_A` / `calc_B` (called back into C++) form the Godambe sandwich
`A^{-1} B A^{-1} / n`.

### KEE estimators

[`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md)
solves a partial-likelihood-style kernel estimating equation via
[`nleqslv::nleqslv`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html);
variance comes from `kee_cox_var` returning `W`, `Sigma` and is
assembled as `W^{-1} Σ W^{-1}`.

[`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
is closed-form — `kee_additive_est` returns the estimate plus `A_est`,
`B_est`, `Sigma_est` in one C++ call.

### Cross-validation

[`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
(R/skmle_cv.R, exported) generates an automatic log-spaced `h_grid` in
`[max(min(pos_diffs), n^-0.6), min(max(max_diffs), n^-0.3)]`, assigns
folds by **subject id** (not by row), and delegates the whole
K×\|h_grid\| loop to `skmle_cv_cpp` in one call — the C++ side refits on
each train split and evaluates `skmle_eval_nll_cpp` on the held-out
subjects. After selecting `h_cv`, it refits on the full data by
rewriting `call[[1]]` to
[`skmle::skmle`](https://dayusun.github.io/skmle/reference/skmle.md) and
stripping CV-only arguments.

### Simulation

[`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md)
(`R/simulate.R`) draws a bivariate covariate (step-function AR-like
process + binary indicator), inverts the transformed cumulative hazard
via Legendre-Gauss quadrature + `nleqslv`, and draws observation times
from a thinned NHPP (`NonHPP_gen`). **Only two covariates are
supported** — `length(beta)` must be 2.

### Long-running C++ loops

Outer loops in `src/skmle_cpp.cpp` and `src/kee.cpp` call
`Rcpp::checkUserInterrupt()` every ~100 iterations. Preserve this when
adding new loops so Ctrl-C still works.

# Benchmarking skmle Against SurvSparse

## Overview

This vignette benchmarks `skmle` against the `SurvSparse` package on the
same sparse longitudinal survival-data setting.

We compare:

1.  [`SurvSparse::add.haz()`](https://rdrr.io/pkg/SurvSparse/man/add.haz.html)
    against
    [`skmle::kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
    and `skmle::skmle(s = 1)`.
2.  [`SurvSparse::trans.haz()`](https://rdrr.io/pkg/SurvSparse/man/trans.haz.html)
    against `skmle::skmle(s = 0)`.

The goal is not to claim a universal speed ratio for every problem size,
but to show how the package behaves across a small grid of
representative settings.

``` r
library(nloptr)
library(skmle)
library(SurvSparse)
library(survival)
library(dplyr)
library(bench)
```

## Data Generation

We simulate sparse longitudinal survival data using
[`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md).
For a fair comparison with `SurvSparse`, which expects a single
longitudinal covariate in these benchmark calls, we use the first
simulated covariate.

``` r
make_benchmark_data <- function(n, s_val) {
  dat <- sim_skmle_data(
    n = n,
    mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
    mu_bar = 8,
    alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
    beta = c(1, -0.5),
    s = s_val,
    cen = 1.0
  )

  dat %>%
    dplyr::rename(covariates_old = covariates) %>%
    dplyr::mutate(
      covariates = covariates_old[, 1],
      X = X,
      delta = delta,
      obs_times = obs_times,
      id = as.integer(factor(id))
    ) %>%
    dplyr::select(id, X, covariates, obs_times, delta) %>%
    as.data.frame()
}
```

## Benchmark Helpers

``` r
run_additive_benchmark <- function(n, iterations = 5) {
  dat <- make_benchmark_data(n, s_val = 1)
  h_val <- n^(-0.5)

  res <- bench::mark(
    SurvSparse_add_haz = add.haz(
      data = dat, n = n, tau = 1, h = h_val, method = 3
    ),
    skmle_kee_additive = kee_additive(
      Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times, h = h_val
    ),
    skmle_spline = skmle(
      Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times,
      s = 1, h = h_val, nknots = 3, norder = 3
    ),
    iterations = iterations,
    check = FALSE
  )

  as.data.frame(res[, c("expression", "median", "itr/sec")]) %>%
    dplyr::mutate(expression = as.character(expression)) %>%
    dplyr::mutate(
      scenario = sprintf("Additive (n = %d)", n),
      median_ms = as.numeric(median) / 1e6
    ) %>%
    dplyr::select(scenario, expression, median_ms, `itr/sec`)
}

run_transformed_benchmark <- function(n, iterations = 5) {
  dat <- make_benchmark_data(n, s_val = 0)
  h_val <- n^(-0.5)

  res <- bench::mark(
    SurvSparse_trans_haz = trans.haz(
      data = dat, n = n, nknots = 3, norder = 3, tau = 1, s = 0, h = h_val
    ),
    skmle_spline = skmle(
      Surv(X, delta) ~ covariates,
      data = dat, id = id, obs_times = obs_times,
      s = 0, h = h_val, nknots = 3, norder = 3
    ),
    iterations = iterations,
    check = FALSE
  )

  as.data.frame(res[, c("expression", "median", "itr/sec")]) %>%
    dplyr::mutate(expression = as.character(expression)) %>%
    dplyr::mutate(
      scenario = sprintf("Transformed (n = %d)", n),
      median_ms = as.numeric(median) / 1e6
    ) %>%
    dplyr::select(scenario, expression, median_ms, `itr/sec`)
}
```

## Results Across Multiple Scenarios

``` r
set.seed(20260325)

benchmark_results <- dplyr::bind_rows(
  run_additive_benchmark(200, iterations = 5),
  run_additive_benchmark(500, iterations = 5),
  run_transformed_benchmark(100, iterations = 5),
  run_transformed_benchmark(200, iterations = 5)
)

benchmark_results
#>                 scenario           expression    median_ms   itr/sec
#> 1     Additive (n = 200)   SurvSparse_add_haz 3.232287e-07  3.306787
#> 2     Additive (n = 200)   skmle_kee_additive 8.970429e-08 11.005874
#> 3     Additive (n = 200)         skmle_spline 1.272635e-07  7.728497
#> 4     Additive (n = 500)   SurvSparse_add_haz 6.365453e-07  1.384430
#> 5     Additive (n = 500)   skmle_kee_additive 9.459881e-08 10.527855
#> 6     Additive (n = 500)         skmle_spline 1.949731e-07  5.069229
#> 7  Transformed (n = 100) SurvSparse_trans_haz 5.178577e-07  1.938386
#> 8  Transformed (n = 100)         skmle_spline 9.936801e-08  9.897646
#> 9  Transformed (n = 200) SurvSparse_trans_haz 1.002256e-06  1.000100
#> 10 Transformed (n = 200)         skmle_spline 1.155567e-07  8.479270
```

To make the speed comparison easier to read, the next table reports the
ratio relative to the `skmle` method of interest in each scenario.

``` r
baseline_rows <- benchmark_results %>%
  dplyr::filter(
    (grepl("^Additive", scenario) & expression == "skmle_kee_additive") |
      (grepl("^Transformed", scenario) & expression == "skmle_spline")
  ) %>%
  dplyr::transmute(scenario, baseline_ms = median_ms)

speed_summary <- benchmark_results %>%
  dplyr::left_join(baseline_rows, by = "scenario") %>%
  dplyr::mutate(speedup_vs_baseline = median_ms / baseline_ms)

speed_summary
#>                 scenario           expression    median_ms   itr/sec
#> 1     Additive (n = 200)   SurvSparse_add_haz 3.232287e-07  3.306787
#> 2     Additive (n = 200)   skmle_kee_additive 8.970429e-08 11.005874
#> 3     Additive (n = 200)         skmle_spline 1.272635e-07  7.728497
#> 4     Additive (n = 500)   SurvSparse_add_haz 6.365453e-07  1.384430
#> 5     Additive (n = 500)   skmle_kee_additive 9.459881e-08 10.527855
#> 6     Additive (n = 500)         skmle_spline 1.949731e-07  5.069229
#> 7  Transformed (n = 100) SurvSparse_trans_haz 5.178577e-07  1.938386
#> 8  Transformed (n = 100)         skmle_spline 9.936801e-08  9.897646
#> 9  Transformed (n = 200) SurvSparse_trans_haz 1.002256e-06  1.000100
#> 10 Transformed (n = 200)         skmle_spline 1.155567e-07  8.479270
#>     baseline_ms speedup_vs_baseline
#> 1  8.970429e-08            3.603269
#> 2  8.970429e-08            1.000000
#> 3  8.970429e-08            1.418700
#> 4  9.459881e-08            6.728893
#> 5  9.459881e-08            1.000000
#> 6  9.459881e-08            2.061053
#> 7  9.936801e-08            5.211513
#> 8  9.936801e-08            1.000000
#> 9  1.155567e-07            8.673286
#> 10 1.155567e-07            1.000000
```

## Interpretation

The benchmark shows a stable pattern:

- For the additive comparison,
  [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
  is substantially faster than
  [`SurvSparse::add.haz()`](https://rdrr.io/pkg/SurvSparse/man/add.haz.html).
- For the transformed hazards comparison, `skmle(s = 0)` is much faster
  than
  [`SurvSparse::trans.haz()`](https://rdrr.io/pkg/SurvSparse/man/trans.haz.html)
  in these runs.
- The general spline-based `skmle(s = 1)` fit is slower than
  [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md),
  which is expected because it solves the broader joint optimization
  problem rather than a specialized estimating equation.

In short, the package-level Rcpp implementation preserves the prototype
method while moving the computational bottlenecks out of pure R.

# Default bandwidth for the survival estimators

Geometric midpoint of the grid
[`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
searches, which is built from the observed lags between covariate
observation times and event times.

## Usage

``` r
default_bandwidth_surv(X_time, obs_times, n)
```

## Arguments

- X_time, obs_times:

  Event/censoring and covariate observation times.

- n:

  Number of subjects.

## Value

A positive bandwidth.

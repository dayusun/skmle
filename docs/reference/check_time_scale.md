# Validate the time scale

Times may be on any scale. The sieve basis and the cumulative-hazard
quadrature are built on `[0, max(X_time)]`, so nothing has to be
rescaled to the unit interval first; that was a convention of the
original prototype rather than a requirement of the method. What is
required is that times are finite and non-negative, and that follow-up
has positive length.

## Usage

``` r
check_time_scale(X_time, obs_times_vec)
```

## Arguments

- X_time:

  Event or censoring times.

- obs_times_vec:

  Covariate observation times.

## Value

`NULL`, invisibly. Called for the error.

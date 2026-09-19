# Validate the time scale

Times may be on any scale. The sieve basis and the cumulative-hazard
quadrature are built on `[0, max(X_time)]`, so nothing has to be
rescaled to the unit interval first; that was a convention of the
original prototype rather than a requirement of the method. What is
required of `X_time` is that it is finite and non-negative, and that
follow-up has positive length.

## Usage

``` r
check_time_scale(X_time, obs_times_vec)
```

## Arguments

- X_time:

  Event or censoring times.

- obs_times_vec:

  Covariate observation times. May be negative.

## Value

`NULL`, invisibly. Called for the error.

## Details

`obs_times_vec` is held to the weaker standard of being finite, and may
be negative. No basis and no quadrature rule is ever evaluated at a
covariate observation time: it reaches the score only through the lag
`X_time - obs_times_vec`, and the one-sided rule tests the sign of that
lag, not the sign of the time. The earlier non-negativity requirement
applied a constraint belonging to `X_time` to an argument that does not
carry it. A negative observation time means a covariate measured before
the time origin, which is a scientific statement about the data and not
a numerical one; see “Negative observation times” in
[`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md).

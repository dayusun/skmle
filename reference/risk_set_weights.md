# Risk-set weights at the event times

The risk-set averages \\S^{(k)}\\ use the same kernel as the main term
of each estimating equation, so they are built from the same weight
function and the same one-sided rule, here, rather than rewritten in
C++. Rows are the distinct event times among the rows that contribute
(`delta == 1` and a nonzero row weight), columns are the observation
rows, and `index` maps each data row to its row of `K` (0-based, `-1`
for a row that does not contribute).

## Usage

``` r
risk_set_weights(X_time, obs_times, delta, kerval, h, weight, one_sided)
```

## Arguments

- X_time, obs_times, delta:

  Event times, observation times, event indicator.

- kerval:

  Row weights from
  [`kernel_weights()`](https://www.sundayu.me/skmle/reference/epan_kernel.md).

- h, weight, one_sided:

  As for
  [`kernel_weights()`](https://www.sundayu.me/skmle/reference/epan_kernel.md).

## Value

A list with `index` (integer) and `K` (numeric matrix).

## Details

The matrix has one row per distinct event time, not per data row,
because every row of a subject shares the subject's event time.

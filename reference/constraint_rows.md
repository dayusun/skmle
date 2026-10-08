# Rows that carry a feasibility constraint

For `s != 0` the transformed hazard needs `1 + s * eta > 0` wherever the
objective evaluates it with nonzero weight. A row is constrained when
its lag falls in the weight's declared support or its weight is nonzero,
and, under `one_sided`, its lag is positive. For the Epanechnikov
kernels this is exactly the earlier `|lag| <= h` rule; for a weight
reaching further than `h` it follows the weight. `skmle_cv_cpp()`
applies the same rule.

## Usage

``` r
constraint_rows(lag, kerval, h, weight, one_sided)
```

## Arguments

- lag:

  `X - obs_times`.

- kerval:

  Row weights from
  [`kernel_weights()`](https://www.sundayu.me/skmle/reference/epan_kernel.md).

- h, weight, one_sided:

  As for
  [`kernel_weights()`](https://www.sundayu.me/skmle/reference/epan_kernel.md).

## Value

Logical vector.

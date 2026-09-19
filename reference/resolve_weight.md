# Resolve the weight argument

`NULL` means "the Epanechnikov implied by `one_sided`", which is what
every call made before `weight` existed was asking for. Resolving it
here rather than as a default argument keeps that backward compatibility
exact: a literal default of
[`w_epan_half()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)
would silently turn `one_sided = FALSE` into a one-sided fit, because
the support would still be `[0, 1]`.

## Usage

``` r
resolve_weight(weight, one_sided)
```

## Arguments

- weight:

  A weight function, or `NULL`.

- one_sided:

  Logical, used only when `weight` is `NULL`.

## Value

A weight function.

# Default bandwidth for the asynchronous estimators

Geometric midpoint of the grid
[`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
searches, \\2 (Q_3 - Q_1) n^{-1/2}\\, so it adapts to whatever units
`time` is in.

## Usage

``` r
default_bandwidth_async(times, n)
```

## Arguments

- times:

  Pooled observation times from both tables.

- n:

  Number of subjects.

## Value

A positive bandwidth.

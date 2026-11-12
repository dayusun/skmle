# Tell the user which bandwidth was chosen for them

Tell the user which bandwidth was chosen for them

## Usage

``` r
announce_bandwidth(h, cv_fun = NULL)
```

## Arguments

- h:

  The chosen bandwidth.

- cv_fun:

  Name of the cross-validation function to point at, or `NULL` when the
  estimator has none.

## Value

`h`, invisibly.

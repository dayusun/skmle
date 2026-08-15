# Simulate asynchronous longitudinal data

Draw a sparse longitudinal response and a sparse longitudinal covariate
whose observation times are **independent**, so the two are never seen
together. This is the design of Cao, Zeng and Fine (2015), Section 4:
response times and covariate times are independent Poisson streams on
\\(0, 1)\\, the covariate is a mean-zero Gaussian process, and the error
is a second, independent Gaussian process.

## Usage

``` r
sim_async_data(
  n,
  beta = c(0.5, 1.5),
  lambda_y = 5,
  lambda_x = 5,
  link = c("identity", "log", "logistic"),
  x_cov = function(s, t) exp(-abs(s - t)),
  e_cov = function(s, t) 2^(-abs(s - t))
)
```

## Arguments

- n:

  Number of subjects.

- beta:

  Either a numeric vector `c(intercept, slopes)` or a function of a
  numeric time vector returning a matrix with `length(beta(t))` columns
  and one row per time.

- lambda_y, lambda_x:

  Poisson rates for the response and covariate observation streams on
  `(0, 1)`. Each subject is forced to have at least one observation of
  each.

- link:

  One of `"identity"`, `"log"` or `"logistic"`, matching
  [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md).
  Under `"logistic"` the response is Bernoulli and no error process is
  added; under `"log"` the response is Poisson.

- x_cov, e_cov:

  Covariance functions of two times, vectorised, for the covariate and
  error processes.

## Value

A list with two data frames on different time grids:

- y:

  `id`, `time`, `y` – one row per response occasion.

- x:

  `id`, `time`, and one column per covariate coordinate – named `x` when
  there is a single covariate, `x1`, `x2`, ... otherwise. One row per
  covariate occasion, on a time grid independent of `y`'s.

## Details

The default covariate covariance \\\exp(-\|s - t\|)\\ is the Ornstein-
Uhlenbeck covariance. Its diagonal has a kink, so it is only one-sided
differentiable, which is exactly the case Cao, Zeng and Fine's Section 6
leaves open. Supply a smooth `x_cov`, for instance
`function(s, t) exp(-(s - t)^2)`, to generate data inside their stated
assumptions instead.

`beta` may be a numeric vector, giving time-invariant coefficients
suitable for
[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md),
or a function of time returning one row per supplied time, giving a
coefficient curve suitable for
[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md).
In both cases the first element is the intercept and the remaining `p`
elements multiply the `p` covariate coordinates.

## References

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## Examples

``` r
set.seed(1)
d <- sim_async_data(n = 50)
str(d)
#> List of 2
#>  $ y:'data.frame':   257 obs. of  3 variables:
#>   ..$ id  : int [1:257] 1 1 1 1 2 2 2 2 2 2 ...
#>   ..$ time: num [1:257] 0.202 0.573 0.898 0.908 0.108 ...
#>   ..$ y   : num [1:257] -2.945 -0.819 0.211 0.33 2.494 ...
#>  $ x:'data.frame':   258 obs. of  3 variables:
#>   ..$ id  : int [1:258] 1 1 1 1 2 2 2 2 2 2 ...
#>   ..$ time: num [1:258] 0.0618 0.6291 0.6608 0.9447 0.0233 ...
#>   ..$ x   : num [1:258] -0.864 0.345 0.437 0.253 0.292 ...

# The two tables share `id` and nothing else; their `time` grids are
# independent draws, which is what makes the data asynchronous.
head(d$y, 3)
#>   id      time          y
#> 1  1 0.2016819 -2.9454025
#> 2  1 0.5728534 -0.8186868
#> 3  1 0.8983897  0.2109702
head(d$x, 3)
#>   id       time          x
#> 1  1 0.06178627 -0.8642242
#> 2  1 0.62911404  0.3452775
#> 3  1 0.66079779  0.4373871

# A coefficient curve, for kee_async_td().
d2 <- sim_async_data(n = 50, beta = function(tt) cbind(0.5, 1 + tt))
```

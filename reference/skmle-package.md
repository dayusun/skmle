# skmle: Sieve Kernel Maximum Likelihood Estimation

Estimation and inference when a longitudinal covariate is observed
sparsely and intermittently, at times that do not line up with the times
at which the outcome is measured. For a survival outcome, fits the
transformed hazards family – proportional hazards, additive hazards, and
the Box-Cox transformations between them – by the Sieve Maximum
Kernel-weighted Log-likelihood Estimator (SMKLE) of Sun, Sun, Zhao and
Cao (2025)
[doi:10.1080/01621459.2025.2476781](https://doi.org/10.1080/01621459.2025.2476781)
, with specialised kernel estimating-equation alternatives for the
proportional and additive cases. For a longitudinal outcome recorded on
its own time grid, fits generalised linear models by the kernel-weighted
estimating equations of Cao, Zeng and Fine (2015)
[doi:10.1111/rssb.12086](https://doi.org/10.1111/rssb.12086) , with
either time-invariant or time-dependent coefficients. Half and full
kernels, data simulation, and bandwidth selection by cross-validation
are available throughout, and the numerical routines use 'Rcpp',
'RcppArmadillo' and the 'nloptr' C API.

## Details

A longitudinal covariate is rarely measured at the time you need it.
`skmle` handles that mismatch directly, by weighting each observation
according to how far its measurement time sits from the time being
modelled, rather than carrying a value forward or smoothing the
covariate and substituting it. Two outcome types are covered.

## Survival outcomes

[`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) fits the
transformed hazards family by the Sieve Maximum Kernel-weighted
Log-likelihood Estimator of Sun, Sun, Zhao and Cao (2025). The Box-Cox
parameter `s` indexes the family: `s = 0` is proportional hazards,
`s = 1` is additive hazards, and other values interpolate.
[`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md) and
[`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
are faster specialised estimating equations for the two named cases, and
[`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
selects the bandwidth by subject-level cross-validation. Simulate with
[`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md).

## Asynchronous longitudinal outcomes

When the outcome is itself a sparsely observed longitudinal process,
recorded on a time grid that does not line up with the covariate's,
[`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
and
[`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
fit generalised linear models by the kernel-weighted estimating
equations of Cao, Zeng and Fine (2015), with time-invariant coefficients
and with a coefficient curve \\\beta(t)\\ respectively. Simulate with
[`sim_async_data()`](https://dayusun.github.io/skmle/reference/sim_async_data.md).

## Half and full kernels

Every estimator takes `one_sided`. A half kernel admits only
measurements strictly before the time being modelled, which is the
risk-set restriction of a hazard model and the causal reading of a
covariate path; a full kernel smooths from both sides. The survival
estimators default to the half kernel, as published; the asynchronous
ones default to the full kernel, as published.

## References

Sun, D., Sun, Z., Zhao, X. and Cao, H. (2025). Kernel meets sieve:
transformed hazards models with sparse longitudinal covariates. *Journal
of the American Statistical Association* 120, 2580-2591.

Cao, H., Zeng, D. and Fine, J. P. (2015). Regression analysis of sparse
asynchronous longitudinal data. *Journal of the Royal Statistical
Society, Series B* 77, 755-776.

## See also

Useful links:

- <https://www.sundayu.me/skmle/>

- <https://github.com/dayusun/skmle>

- Report bugs at <https://github.com/dayusun/skmle/issues>

## Author

**Maintainer**: Dayu Sun <dayusun@iu.edu> \[copyright holder\]

Authors:

- Dayu Sun <dayusun@iu.edu> \[copyright holder\]

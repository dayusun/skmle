# skmle 0.1.0

* Initial CRAN release.
* `skmle()` fits transformed hazards models by sieve maximum kernel-weighted
  log-likelihood estimation (SMKLE).
* `kee_cox()` and `kee_additive()` fit Cox and additive hazards models via
  kernel-weighted estimating equations.
* `skmle_cv()` selects the kernel bandwidth by K-fold cross-validation.
* `sim_skmle_data()` simulates survival data with sparse, intermittently
  observed longitudinal covariates.
* `plot()`, `print()` and `summary()` methods for fitted objects.

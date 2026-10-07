Simulation for `kee_async_cv(method = "sync")` and `kee_async(correct = TRUE)`.

- `cov2.R` ... `cov6.R`: pair-level implementation (each sources the previous);
  `cov6.R` is the method as shipped (grid floor ceiling(sqrt(n_train))).
- `run6.R`, `tasks6.txt`, `sim6.sbatch`: BigRed200 job 8466882, five links,
  n = 400, 2000, 8000, 1000 replicates per cell.
- `agg6.R` -> `summary6.csv`: bias, ESE, SE/ESE, RMSE, CP with binomial CI,
  p and BH q for plain cv, cv+corr, cvc+corr (the shipped method) and LVCF.
- `eqv.R`: the package reproduces `cov6.R` (same h, estimates to 1e-15).
- `theory_note.md`: working note with the derivations and all runs.

source("tests/verify_skmle.R")

cat('R f(0):', f(rep(0, 5)), '\n')

estres_pkg_0 <- skmle(
  formula = X_Surv ~ cov1 + cov2,
  data = simdata,
  id = id,
  obs_times = obs_times,
  s = s_val,
  h = hh,
  nknots = 3,
  norder = 3,
  maxeval = 1
)
cat('C++ nll(0):', -estres_pkg_0$loglik, '\n')

# rely on the development helper in the project directory
source("dev/Rprototype_funcs.R")
library(skmle)

set.seed(123)
nn <- 100
true_coef <- c(1, -0.5)
s_val <- 0.5
cen_val <- 0.7002909

simdata <-
  replicate(
    nn,
    simAsytransdata(
      mu = function(tt) 8 * (0.75 + (0.5 - tt) ^ 2),
      mu_bar = 8,
      alpha = function(tt) (1+s_val)/2*0.75 + 0.75*(tt * (1 - sin(2 * pi * (tt - 0.25)))),
      beta =true_coef ,
      s = s_val,
      cen= cen_val,
      nstep = 20
    ),
    simplify = F
  ) %>% bind_rows(.id = "id")

hh <- nn^(-0.4)

# Original R version
cat("Running original R estimation...\n")
ptm <- proc.time()
estres_ori <- estproc_ori_d_mul(
  simdata,
  nn,
  nknots = 3,
  norder = 3,
  s = s_val,
  h = hh
)
ori_time <- proc.time() - ptm
cat("Original R time:\n")
print(ori_time)


# Package version
# break the matrix into separate covariate columns (pkg can also handle a
# matrix variable in the formula, but the old script did not).
simdata$cov1 <- simdata$covariates[,1]
simdata$cov2 <- simdata$covariates[,2]
simdata$X_Surv <- survival::Surv(simdata$X, simdata$delta)

cat("Running skmle package estimation...\n")
ptm <- proc.time()
estres_pkg <- skmle(
  formula = X_Surv ~ cov1 + cov2,
  data = simdata,
  id = id,
  obs_times = obs_times,
  s = s_val,
  h = hh,
  nknots = 3,
  norder = 3
)
pkg_time <- proc.time() - ptm
cat("skmle Package time:\n")
print(pkg_time)

cat("\n--- Comparison ---\n")
cat("Original beta & gamma estimates:\n")
print(estres_ori$est)

cat("Package beta estimates:\n")
print(estres_pkg$coefficients)
cat("Package gamma estimates:\n")
print(estres_pkg$gamma)

cat("\nOriginal standard errors (se):\n")
print(estres_ori$se)
cat("Package standard errors (se):\n")
print(sqrt(diag(estres_pkg$var)))

cat("Package var:\n")
print(estres_pkg$var)

cat("Orig A_est:\n")
print(estres_ori$A_est)
cat("Pkg A_est:\n")
print(estres_pkg$A)

cat("Orig B_est:\n")
print(estres_ori$B_est)
cat("Pkg B_est:\n")
print(estres_pkg$B)


cat("Are betas & gammas identical? ", all.equal(estres_ori$est, c(estres_pkg$coefficients, estres_pkg$gamma), tolerance=1e-4, check.attributes = FALSE), "\n")
cat("Are standard errors identical? ", all.equal(estres_ori$se, sqrt(diag(estres_pkg$var)), tolerance=1e-4, check.attributes = FALSE), "\n")

cat("\n--- Summary output ---\n")
summary(estres_pkg)

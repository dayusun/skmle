# rely on the development helper in the project directory
source("dev/Rprototype_funcs.R")
library(skmle)

set.seed(42)
nn <- 50
true_coef <- c(1, -0.5)
s_val <- 0.5
cen_val <- 0.7002909

# Simulate data
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
    simplify = FALSE
  )
simdata <- do.call(rbind, simdata)
simdata$id <- rep(1:nn, each = 20)[1:nrow(simdata)] # This is roughly right but let's use actual bind_rows to be safe if dplyr is available, but the original used bind_rows.

# We'll use dplyr's bind_rows since it's already used in the project
library(dplyr)
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
    simplify = FALSE
  ) %>% bind_rows(.id = "id")

simdata$id <- as.numeric(simdata$id)
simdata$cov1 <- simdata$covariates[,1]
simdata$cov2 <- simdata$covariates[,2]
simdata$X_Surv <- survival::Surv(simdata$X, simdata$delta)

cat("Running skmle_cv cross-validation estimation...\n")
ptm <- proc.time()

# Provide a small grid
h_grid_test <- c(0.1, 0.2, 0.3)

cv_res <- skmle_cv(
  formula = X_Surv ~ cov1 + cov2,
  data = simdata,
  id = id,
  obs_times = obs_times,
  s = s_val,
  K = 3,
  h_grid = h_grid_test,
  nknots = 3,
  norder = 3,
  quiet = FALSE
)

cv_time <- proc.time() - ptm
cat("skmle_cv Package time:\n")
print(cv_time)

cat("\n--- Comparison & Validation ---\n")
cat("Returned class:", class(cv_res), "\n")
stopifnot(inherits(cv_res, "cv.skmle"))

cat("CV Results:\n")
print(cv_res$cv_results)

best_h_idx <- which.min(cv_res$cv_results$cvloss)
expected_h <- cv_res$cv_results$h[best_h_idx]

cat("Expected optimal h:", expected_h, "\n")
cat("Returned h_cv:", cv_res$h_cv, "\n")

stopifnot(all.equal(expected_h, cv_res$h_cv))

cat("Does the refitted model use the optimal h? ")
fit_h <- cv_res$fit$h
cat(fit_h, "\n")
stopifnot(all.equal(fit_h, cv_res$h_cv))

cat("\nCross-validation functionality passes basic verification!\n")

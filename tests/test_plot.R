library(skmle)
library(dplyr)
source("dev/Rprototype_funcs.R")

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
simdata$cov1 <- simdata$covariates[,1]
simdata$cov2 <- simdata$covariates[,2]
simdata$X_Surv <- survival::Surv(simdata$X, simdata$delta)

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

p <- plot(estres_pkg)
ggplot2::ggsave("tests/test_plot.png", plot = p, width = 6, height = 4)
cat("Plot saved to tests/test_plot.png\n")

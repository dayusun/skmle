library(devtools)
load_all()

source("tests/multiple_Cao_Cox_funcs.R")

# Generate 1 dataset
set.seed(123)
nn <- 200
s_val <- 0
cenrate <- 20
cen_val <- 0.7066465 # From cen_mat 1,1
true_coef <- c(1, -0.5)

simdata <- replicate(
    nn,
    simAsytransdata(
        mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
        mu_bar = 8,
        alpha = function(tt) (1 + s_val) / 2 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
        beta = true_coef,
        s = s_val,
        cen = cen_val,
        nstep = 20
    ),
    simplify = FALSE
) %>% bind_rows(.id = "id")

simdata$id <- as.numeric(simdata$id)

hh <- nn^(-0.4)

# 1. Original R implementation
cat("Running original R implementation...\n")
res_ori <- estproc_Cox_Cao_mul(data = simdata, n = nn, h = hh)

# 2. New kee_cox C++ implementation
cat("Running kee_cox (C++ via R wrapper)...\n")

# Format for formula
# Need covariates as separate columns
simdata_formula <- simdata
simdata_formula$Z1 <- simdata$covariates[, 1]
simdata_formula$Z2 <- simdata$covariates[, 2]

res_kee <- kee_cox(
    formula = Surv(X, delta) ~ Z1 + Z2,
    data = simdata_formula,
    id = id,
    obs_times = obs_times,
    h = hh
)

cat("Original Estimate:\n")
print(res_ori$est)
cat("New Estimate:\n")
print(res_kee$coefficients)

cat("Original se:\n")
print(res_ori$se)
cat("New se:\n")
print(sqrt(diag(res_kee$var)))

cat("Is estimate identical? ", isTRUE(all.equal(as.vector(res_ori$est), as.vector(res_kee$coefficients))), "\n")
cat("Is SE identical? ", isTRUE(all.equal(as.vector(res_ori$se), as.vector(sqrt(diag(res_kee$var))))), "\n")

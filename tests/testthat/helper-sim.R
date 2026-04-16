make_sim <- function(n = 80, s = 0, seed = 1) {
  set.seed(seed)
  sim_skmle_data(
    n = n,
    mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
    mu_bar = 8,
    alpha = function(tt) {
      (1 + s) / 2 * 0.75 +
        0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25))))
    },
    beta = c(1, -0.5),
    s = s,
    cen = 0.7
  )
}

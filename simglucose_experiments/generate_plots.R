rm(list = ls())

source("./utils/taus_from_oracle.R")

eff_days <- c(1, 2, 4, 8, Inf)
gamma_denoms <- 24 * eff_days
gammas <- 1 - 1/(gamma_denoms)

result_150 <- compute_tau(dir = "./oracle",
                          thresholds = seq(120, 180, by=5),
                          eval_at = 150)

plot_oracle(result_150, eval_at = 150)

p <- plot_oracle_one_gamma(
  result_150,
  gamma_target = 1 - 1/96,
  spar = 0.45
)

p2 <- plot_oracle_days_max(
  result_150,
  spar = 0.45,
  n_grid = 50,
  gamma_text = c("0", "1 - 1/24", "1 - 1/48", "1 - 1/96", "1 - 1/192", "1"),
  eff_days_text = c("0", "0.5", "1", "2", "4", "Inf")
)

p
p2
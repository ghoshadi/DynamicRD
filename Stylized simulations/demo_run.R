rm(list=ls())

source("./autoreg_utils.R")

run_all_tables(
    gammas = c(0.5, 0.8, 1.0),
    n_grid = 100*c(1, 2, 4, 8),
    reps = 500,
    c0 = 110, h = NULL,
    Tn = 12, mu = 100, sigma = 4, rho = 0.9, tau = 0.1, mu_trend = 1,
    K = K_uniform, 
    time_fe = TRUE,
    NREP_oracle = 1e6,
    cores = parallel::detectCores(), 
    save_dir = "./tables",
    save_tables = TRUE
)

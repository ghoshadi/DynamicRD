rm(list = ls())

source("./utils/simglucose_cov_width.R")

source("./utils/taus_from_oracle.R")

set.seed(42)
eff_days <- c(1, 2, 4, 8, Inf)
gamma_denoms <- 24 * eff_days
gammas <- 1 - 1/(gamma_denoms)

# ── Oracle computation for c=150 ──────────────────────────────────────────────

result_150 <- compute_tau(dir = "./oracle",
                          thresholds = seq(120, 180, by=5),
                          eval_at = 150)

taus_150 <- subset(result_150$tau, days == max(result_150$tau$days))$tau
names(taus_150) <- c("gamma=0", 
                     paste0("gamma=1-1/", gamma_denoms)
                     )
cat("Oracle tau_RDs (4-day horizon):\n")
print(data.frame(eff_days = c(0, eff_days), 
                 gamma = c(0,round(gammas, 3)), 
                 tau_RD = taus_150)
      )

# ── Run coverage / width simulation ───────────────────────────────────────────
run_sim(
  tau_RDs  = taus_150[-1],
  dir      = "data_demo",
  c0       = 150,
  gammas   = gammas,
  mc.cores = parallel::detectCores()-1
)

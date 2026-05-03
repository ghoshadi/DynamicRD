# ---------------------------------------------------------------------------
#   Codes for the Autoregressive Simulator
#   This runs the stylized simulation experiments in the paper:
#   Non-parametric Causal Inference in Dynamic Thresholding Designs
#   Aditya Ghosh and Stefan Wager (Stanford University)
#   ArXiv link: https://arxiv.org/abs/2512.15244
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(pbmcapply)
})
n.cores = parallel::detectCores() # for parallel computation in R

# --------------------------------------------------
#   Generate data from the Autoregressive process
# --------------------------------------------------
simu_X <- function(n, Tn, mu, sigma, rho, tau, thresh, mu_trend) {
  mu_t <- mu + mu_trend * (0:Tn)
  X <- matrix(NA_real_, n, Tn + 1)
  sigma0 <- sigma / sqrt(1 - rho^2)
  X[, 1] <- mu_t[1] + sigma0 * rnorm(n)
  for (t in 1:Tn) {
    Xt <- X[, t]
    Wt <- as.numeric(Xt >= thresh)
    Xc <- Xt - mu
    X[, t + 1] <- mu_t[t + 1] + rho * Xc - 
      Wt * tau * pmax(Xc, 0) + 
      sigma * rnorm(n)
  }
  X
}

# --------------------------------------------------
#   Kernels (helper functions)
# --------------------------------------------------
K_triangular <- function(u) pmax(1 - abs(u), 0)
K_uniform <- function(u) 0.5 * (abs(u) <= 1)
K_epanechnikov <- function(u) 0.75 * pmax(1 - u**2, 0)
# Or implement any K(u) that integrates to 1 on [-1,1].

# --------------------------------------------------
#   Imbens-Kalyanaraman bandwidth
# --------------------------------------------------
IK_bandwidth <- function(Y, X, threshold, kernel = c("triangular", "uniform")) {
  if (threshold >= max(X) || threshold <= min(X)) {
    stop("RD threshold is outside the running variable range.")
  }
  kernel <- match.arg(kernel)
  
  x <- X - threshold
  n <- length(Y)
  
  h.silverman <- 1.84 * stats::sd(x) * n^(-1/5)
  
  i.plus.1 <- x >= 0 & x <= h.silverman
  i.min.1  <- x < 0 & x >= -h.silverman
  
  y.ave.plus.1 <- mean(Y[i.plus.1])
  y.ave.min.1  <- mean(Y[i.min.1])
  
  dy.plus.1 <- Y[i.plus.1] - y.ave.plus.1
  dy.min.1  <- Y[i.min.1]  - y.ave.min.1
  
  sigmas <- (sum(dy.plus.1^2) + sum(dy.min.1^2)) / (sum(i.plus.1) + sum(i.min.1))
  fc <- (sum(i.plus.1) + sum(i.min.1)) / (2 * n * h.silverman)
  
  # third derivative
  median.plus <- stats::median(x[x >= 0])
  median.min  <- stats::median(x[x < 0])
  
  middle <- x >= median.min & x <= median.plus
  x.mid <- x[middle]
  y.mid <- Y[middle]
  
  tt <- cbind(1, x.mid >= 0, x.mid, x.mid^2, x.mid^3)
  gamma <- solve(t(tt) %*% tt, t(tt) %*% y.mid)
  third.der <- 6 * gamma[5]
  
  # second derivatives
  h.plus.2 <- 3.56 * (sigmas / (fc * max(third.der^2, 0.01)))^(1/7) * sum(x >= 0)^(-1/7)
  h.min.2  <- 3.56 * (sigmas / (fc * max(third.der^2, 0.01)))^(1/7) * sum(x < 0)^(-1/7)
  
  i.min.3  <- x < 0 & x >= -h.min.2
  i.plus.3 <- x >= 0 & x <= h.plus.2
  if (sum(i.min.1) <= 2 || sum(i.plus.3) <= 2) {
    stop("Insufficient observations near discontinuity.")
  }
  
  # left second derivative
  x.min <- x[i.min.3]
  y.min <- Y[i.min.3]
  t.min <- cbind(1, x.min, x.min^2)
  beta.min <- solve(t(t.min) %*% t.min, t(t.min) %*% y.min)
  second.der.min <- 2 * beta.min[3]
  
  # right second derivative
  x.plus <- x[i.plus.3]
  y.plus <- Y[i.plus.3]
  t.plus <- cbind(1, x.plus, x.plus^2)
  beta.plus <- solve(t(t.plus) %*% t.plus, t(t.plus) %*% y.plus)
  second.der.plus <- 2 * beta.plus[3]
  
  # regularization terms
  r.plus <- 2160 * sigmas / (sum(i.plus.3) * h.plus.2^4)
  r.min  <- 2160 * sigmas / (sum(i.min.3)  * h.min.2^4)
  
  CK <- 3.4375
  denom <- (second.der.plus - second.der.min)^2 + r.plus + r.min
  h.opt <- CK * ((2 * sigmas) / (fc * denom))^(1/5) * n^(-1/5)
  
  # Compute weights
  dist <- abs((X - threshold) / h.opt)
  if (kernel == "triangular") {
    weights <- (1 - dist) * (dist <= 1) / h.opt
    weights <- weights / sum(weights)
  } else {
    weights <- (1 - dist) * (dist <= 1)
    weights[weights > 0 ] <- 1
  }
  
  list(
    bandwidth = h.opt,
    weights = weights
  )
}


# --------------------------------------------------
# helper: build design matrix
# --------------------------------------------------
.build_X <- function(Zk, Wk, tk, Tn, time_fe = FALSE,
                     time_fe_Z = FALSE, time_fe_ZxW = FALSE) {
  X <- cbind(Intercept = 1.0, Z = Zk, W = Wk, ZxW = Zk * Wk)
  if (Tn > 1L && (time_fe || time_fe_Z || time_fe_ZxW)) {
    D <- model.matrix(~ factor(tk, levels = seq_len(Tn)))[, -1, drop = FALSE]
    colnames(D) <- paste0("tFE_", seq_len(Tn - 1L))
    if (time_fe) X <- cbind(X, D)
    if (time_fe_Z) {
      DZ <- sweep(D, 1, Zk, "*")
      colnames(DZ) <- paste0("Z:tFE_", seq_len(Tn - 1L))
      X <- cbind(X, DZ)
    }
    if (time_fe_ZxW) {
      DZxW <- sweep(D, 1, Zk * Wk, "*")
      colnames(DZxW) <- paste0("ZxW:tFE_", seq_len(Tn - 1L))
      X <- cbind(X, DZxW)
    }
  }
  X
}

# --------------------------------------------------
#   helper: Cholesky solver
# --------------------------------------------------
chol_solve <- function(A, B = diag(nrow(A))) {
  R <- chol(A)
  backsolve(R, forwardsolve(t(R), B))
}

# --------------------------------------------------
#   helper: Fast WLS + cluster VCV
# --------------------------------------------------
.fast_wls_cr <- function(Xw, yw, ck) {
  A <- crossprod(Xw)
  A <- A + 1e-6 * max(diag(A), 1e-6) * diag(ncol(Xw))  # relative ridge for stability
  Ainv <- chol_solve(A)
  beta <- Ainv %*% crossprod(Xw, yw)
  resid <- yw - Xw %*% beta
  U <- rowsum(Xw * as.numeric(resid), group = ck, reorder = FALSE)  # cluster sums of X'e
  meat <- crossprod(U)
  vcov <- Ainv %*% meat %*% Ainv
  list(beta = drop(beta), vcov = vcov, Ainv = Ainv, U = U)
}

# --------------------------------------------------
#   Standard (static) Local Linear Regression
# --------------------------------------------------
rdd_static <- function(Xmat, thresh, h, K = K_uniform, 
                       time_fe = FALSE, time_fe_Z = FALSE, time_fe_ZxW = FALSE) {
  stopifnot(is.matrix(Xmat), ncol(Xmat) >= 2L, h > 0)
  n  <- nrow(Xmat); Tn <- ncol(Xmat) - 1
  X_t  <- Xmat[, 1:Tn, drop = FALSE]
  Ymat <- -Xmat[, 2:(Tn + 1), drop = FALSE]
  
  Z <- as.vector(X_t) - thresh
  W <- as.numeric(Z >= 0)
  wK <- K(Z / h)
  keep <- wK > 0
  if (!any(keep)) stop("No obs in bandwidth.")
  
  Zk <- Z[keep]; Wk <- W[keep]
  Yk <- as.vector(Ymat)[keep]
  ck <- rep.int(seq_len(n), times = Tn)[keep]
  tk <- rep(seq_len(Tn), each = n)[keep]
  
  X   <- .build_X(Zk, Wk, tk, Tn, 
                  time_fe = time_fe, 
                  time_fe_Z = time_fe_Z, 
                  time_fe_ZxW = time_fe_ZxW)
  sw  <- sqrt(wK[keep])
  Xw  <- X * sw
  yw  <- Yk * sw
  
  keep_cols <- colSums(abs(Xw)) > 0
  if (!any(keep_cols)) {
    stop("All regressors have zero variation in the bandwidth.")
  }
  if (!keep_cols[which(colnames(X) == "W")]) {
    stop("No variation in W within the bandwidth (cannot estimate static RD).")
  }
  
  Xw_red <- Xw[, keep_cols, drop = FALSE]
  X_red  <- X[,  keep_cols, drop = FALSE]
  iW_red <- which(colnames(X_red) == "W")
  fit <- .fast_wls_cr(Xw_red, yw, ck)
  list(est = fit$beta[iW_red], se = sqrt(fit$vcov[iW_red, iW_red]))
  
  # fit <- .fast_wls_cr(Xw, yw, ck)
  # iW  <- which(colnames(X) == "W")
  # list(est = fit$beta[iW], se = sqrt(fit$vcov[iW, iW]))
}

# --------------------------------------------------
#   Dynamic RD ratio (proposed method)
# --------------------------------------------------
rdd_dynamic <- function(Xmat, thresh, h, disc_factor = 1, 
                        weights_by_time = TRUE,
                        K = K_uniform, time_fe = FALSE, 
                        time_fe_Z = FALSE, time_fe_ZxW = FALSE) {
  stopifnot(is.matrix(Xmat), ncol(Xmat) >= 2L, h > 0)
  n  <- nrow(Xmat); Tn <- ncol(Xmat) - 1
  
  X_t  <- Xmat[, 1:Tn, drop = FALSE]
  Ymat <- -Xmat[, 2:(Tn + 1), drop = FALSE]
  
  # discounted reverse cumsums
  G <- Ymat
  H <- (X_t >= thresh) * 1.0
  if (Tn > 1 && disc_factor != 0) {
    for (t in (Tn - 1):1) {
      G[, t] <- Ymat[, t] + disc_factor * G[, t + 1]
      H[, t] <- H[, t] + disc_factor * H[, t + 1]
    }
  }
  
  Z <- as.vector(X_t) - thresh
  W <- as.numeric(Z >= 0)
  wK <- K(Z / h)
  keep <- wK > 0
  if (!any(keep)) stop("No obs in bandwidth.")
  
  Zk <- Z[keep]; Wk <- W[keep]
  Gk <- as.vector(G)[keep]
  Hk <- as.vector(H)[keep]
  ck <- rep.int(seq_len(n), times = Tn)[keep]
  tk <- rep(seq_len(Tn), each = n)[keep]
  
  # time weights
  if (weights_by_time) {
    pow <- disc_factor^(0:(Tn - 1))
    wt_time <- rep(pow, each = n)[keep]
  } else {
    wt_time <- 1.0
  }
  wt <- wK[keep] * wt_time
  sw <- sqrt(wt)
  
  X  <- .build_X(Zk, Wk, tk, Tn, 
                 time_fe = time_fe, 
                 time_fe_Z = time_fe_Z, 
                 time_fe_ZxW = time_fe_ZxW)
  Xw <- X * sw
  
  keep_cols <- colSums(abs(Xw)) > 0
  if (!any(keep_cols)) {
    stop("All regressors have zero variation in the bandwidth.")
  }
  if (!keep_cols[which(colnames(X) == "W")]) {
    stop("No variation in W within the bandwidth (cannot estimate dynamic RD).")
  }
  Xw <- Xw[, keep_cols, drop = FALSE]
  X  <- X[,  keep_cols, drop = FALSE]

  fitG <- .fast_wls_cr(Xw, Gk * sw, ck)
  fitH <- .fast_wls_cr(Xw, Hk * sw, ck)
  
  iW <- which(colnames(X) == "W")
  bG <- fitG$beta[iW]; bH <- fitH$beta[iW]
  if (!is.finite(bH) || !is.finite(bG)) stop("Non-finite jump estimate (dynamic).")
  if (abs(bH) < 0.01) warning(sprintf("Dynamic: |bH| = %.4f is small; ratio may be unstable.", abs(bH)))
  
  VG <- fitG$vcov[iW, iW]
  VH <- fitH$vcov[iW, iW]
  
  # cross-cov via cluster influences (unchanged)
  inflG <- fitG$U %*% t(fitG$Ainv)
  inflH <- fitH$U %*% t(fitH$Ainv)
  Cov_GH <- sum(inflG[, iW] * inflH[, iW])
  
  tau_hat <- bG / bH
  var_tau <- (VG + tau_hat^2 * VH - 2 * tau_hat * Cov_GH) / (bH^2)
  se_tau  <- sqrt(max(var_tau, 0))
  
  list(est = tau_hat, se = se_tau, VG = VG, VH = VH, Cov_GH = Cov_GH,
       jump_G = bG, jump_H = bH)
}

# --------------------------------------------------
#   Naive RD (cross-sectional at t = 1)
#   LLR of G_{i,1} = sum_{s=1}^T gamma^{s-1} Y_{is}
#       and H_{i,1} = sum_{s=1}^T gamma^{s-1} A_{is}  on  Z_{i,1}
# --------------------------------------------------
rdd_naive <- function(Xmat, thresh, h, disc_factor = 1, K = K_uniform) {
  stopifnot(is.matrix(Xmat), ncol(Xmat) >= 2L, h > 0)
  n  <- nrow(Xmat); Tn <- ncol(Xmat) - 1L

  X_t  <- Xmat[, 1:Tn,        drop = FALSE]
  Ymat <- -Xmat[, 2:(Tn + 1), drop = FALSE]

  # discounted forward cumsums (same recursion as rdd_dynamic)
  Gmat <- Ymat
  Hmat <- (X_t >= thresh) * 1.0
  if (Tn > 1L && disc_factor != 0) {
    for (t in (Tn - 1L):1L) {
      Gmat[, t] <- Ymat[, t] + disc_factor * Gmat[, t + 1L]
      Hmat[, t] <- Hmat[, t] + disc_factor * Hmat[, t + 1L]
    }
  }

  # cross-sectional slice at t = 1
  Z1   <- X_t[, 1L] - thresh
  W1   <- as.numeric(Z1 >= 0)
  wK   <- K(Z1 / h)
  keep <- wK > 0
  if (!any(keep)) stop("No obs in bandwidth (naive).")

  Zk <- Z1[keep];      Wk <- W1[keep]
  Gk <- Gmat[keep, 1L]; Hk <- Hmat[keep, 1L]
  ck <- seq_len(n)[keep]          # one cluster per unit (cross-sectional)
  sw <- sqrt(wK[keep])

  X  <- cbind(Intercept = 1.0, Z = Zk, W = Wk, ZxW = Zk * Wk)
  Xw <- X * sw

  keep_cols <- colSums(abs(Xw)) > 0
  if (!keep_cols[which(colnames(X) == "W")]) stop("No variation in W (naive).")
  Xw <- Xw[, keep_cols, drop = FALSE]
  X  <- X[,  keep_cols, drop = FALSE]
  iW <- which(colnames(X) == "W")

  fitG <- .fast_wls_cr(Xw, Gk * sw, ck)
  fitH <- .fast_wls_cr(Xw, Hk * sw, ck)

  bG <- fitG$beta[iW]; bH <- fitH$beta[iW]
  if (!is.finite(bH) || !is.finite(bG)) stop("Non-finite jump estimate (naive).")
  if (abs(bH) < 0.1) warning(sprintf("Naive: |bH| = %.4f is small; ratio may be unstable.", abs(bH)))

  VG  <- fitG$vcov[iW, iW]; VH <- fitH$vcov[iW, iW]
  inflG  <- fitG$U %*% t(fitG$Ainv)
  inflH  <- fitH$U %*% t(fitH$Ainv)
  Cov_GH <- sum(inflG[, iW] * inflH[, iW])

  tau_hat <- bG / bH
  var_tau  <- (VG + tau_hat^2 * VH - 2 * tau_hat * Cov_GH) / bH^2
  se_tau   <- sqrt(max(var_tau, 0))

  list(est = tau_hat, se = se_tau)
}

# --------------------------------------------------
#   Oracle value function computation
#   Returns E[ sum_t γ^t Y_{t+1} ] and E[ sum_t γ^t 1{Z_t ≥ c} ]
# --------------------------------------------------
oracle_values <- function(thresholds, 
                          NREP_oracle, 
                          disc_factor, 
                          Tn, mu, rho, tau, sigma,
                          batch_oracle = 2e5,
                          mc.cores = n.cores,
                          base_seed = 42,
                          mu_trend = 0) {
  w <- disc_factor^(0:(Tn - 1))
  idx <- seq_along(thresholds)
  seeds <- base_seed + seq_along(idx)
  one_thr <- function(i) {
    withr::with_seed(seeds[i],{
      thr <- thresholds[i]
      done <- 0L
      sG <- 0.0; sH <- 0.0
      while (done < NREP_oracle) {
        nb <- min(batch_oracle , NREP_oracle - done)
        X <- simu_X(n = nb, Tn = Tn, mu = mu, sigma = sigma, 
                                  rho = rho, tau = tau, thresh = thr,
                                  mu_trend = mu_trend)
        Z <- X[, 1:Tn, drop = FALSE]
        Y <- -X[, 2:(Tn + 1), drop = FALSE]
        sG <- sG + sum(Y %*% w)
        sH <- sH + sum((Z >= thr) %*% w)
        done <- done + nb
      }
    })
    c(thresh = thr, ValG = sG / NREP_oracle, ValH = sH / NREP_oracle)
  }
  
  out <- if (mc.cores > 1) {
    parallel::mclapply(idx, one_thr, mc.cores = min(mc.cores, 2 * length(idx)), mc.set.seed = FALSE)
  } else {
    lapply(idx, one_thr)
  }
  as.data.frame(do.call(rbind, out))
}

# --------------------------------------------------
#   tau_RD oracle computation
#   tau_RD(c0) = dValG/dc / dValH/dc  (density cancels)
# --------------------------------------------------
compute_tauRD <- function(c0, grid, Tn, mu, rho, tau, sigma, disc_factor,
                          NREP_oracle = 1e7, batch_oracle = 2e5, cores = 8,
                          base_seed = 42, mu_trend = 1) {
  oo <- oracle_values(thresholds = grid, NREP_oracle = NREP_oracle, disc_factor = disc_factor,
                      Tn = Tn, mu = mu, rho = rho, tau = tau, sigma = sigma,
                      batch_oracle = batch_oracle, mc.cores = cores, base_seed = base_seed,
                      mu_trend = mu_trend)
  sG <- smooth.spline(oo$thresh, oo$ValG, df = 4)
  sH <- smooth.spline(oo$thresh, oo$ValH, df = 4)
  dG <- predict(sG, c0, deriv = 1)$y
  dH <- predict(sH, c0, deriv = 1)$y
  dG / dH
}

# --------------------------------------------------
#   Run experiments in parallel
# --------------------------------------------------
bench_run <- function(n, Tn, reps, thresh, h, mu, sigma, rho, tau,
                      disc_factor = 1, weights_by_time = TRUE,
                      K = K_uniform, mc.cores = 8,
                      seed_bench = 42, time_fe = FALSE, 
                      time_fe_Z = FALSE, time_fe_ZxW = FALSE,
                      mu_trend = 1,
                      save_errors = TRUE) {
  
  one_rep <- function(i) {
    tryCatch({
      withr::with_seed(seed_bench + i, {
        Xmat <- simu_X(n, Tn, mu, sigma, rho, tau, thresh, mu_trend)
        X_t  <- Xmat[, 1:Tn,        drop = FALSE]
        Ymat <- -Xmat[, 2:(Tn + 1), drop = FALSE]
        
        h_stat  <- h
        h_dyn   <- h
        h_naive <- h
        
        if (is.null(h)) {
          kernel_name <- if (identical(K, K_triangular)) "triangular" else "uniform"
          
          h_stat <- IK_bandwidth(
            Y         = as.vector(Ymat),
            X         = as.vector(X_t),
            threshold = thresh,
            kernel    = kernel_name
          )$bandwidth
          
          Gmat <- Ymat
          if (Tn > 1L && disc_factor != 0) {
            for (t in (Tn - 1L):1L) {
              Gmat[, t] <- Ymat[, t] + disc_factor * Gmat[, t + 1L]
            }
          }
          
          h_dyn <- IK_bandwidth(
            Y         = as.vector(Gmat),
            X         = as.vector(X_t),
            threshold = thresh,
            kernel    = kernel_name
          )$bandwidth
          
          h_naive <- IK_bandwidth(
            Y         = Gmat[, 1L],
            X         = X_t[, 1L],
            threshold = thresh,
            kernel    = kernel_name
          )$bandwidth
        }
        
        s  <- rdd_static(Xmat, thresh, h_stat, K = K, time_fe = time_fe)
        nv <- rdd_naive(Xmat, thresh, h_naive, disc_factor, K = K)
        d  <- rdd_dynamic(Xmat, thresh, h_dyn, disc_factor, weights_by_time,
                          K = K, time_fe = time_fe,
                          time_fe_Z = time_fe_Z,
                          time_fe_ZxW = time_fe_ZxW)
        
        c(
          static_est  = s$est,
          static_se   = s$se,
          naive_est   = nv$est,
          naive_se    = nv$se,
          dynamic_est = d$est,
          dynamic_se  = d$se,
          failed      = 0,
          rep_id      = i
        )
      })
    }, error = function(e) {
      msg <- conditionMessage(e)
      
      if (save_errors) {
        cat(sprintf(
          "[bench_run failure] rep=%d seed=%d n=%d gamma=%s: %s\n",
          i, seed_bench + i, n, as.character(disc_factor), msg
        ))
        flush.console()
      }
      
      c(
        static_est  = NA_real_,
        static_se   = NA_real_,
        naive_est   = NA_real_,
        naive_se    = NA_real_,
        dynamic_est = NA_real_,
        dynamic_se  = NA_real_,
        failed      = 1,
        rep_id      = i
      )
    })
  }
  
  pbmcapply::pbmclapply(
    X = seq_len(reps),
    FUN = one_rep,
    mc.cores = mc.cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

# ---------- formatting ----------
.format_row <- function(n, cov_dyn, wid_dyn, cov_naive, wid_naive,
                        cov_stat_tauRD, wid_stat, cov_stat_tauPE) {
  sprintf("%6d  %5.1f%%  %7.3f   %5.1f%%  %7.3f   %5.1f%%  %7.3f   %5.1f%%  %7.3f",
          n,
          100 * cov_dyn,       wid_dyn,
          100 * cov_naive,     wid_naive,
          100 * cov_stat_tauRD, wid_stat,
          100 * cov_stat_tauPE, wid_stat)
}

# ---------- precompute tau_RD(gamma) ----------
.precompute_tauRDs <- function(gammas, c0, grid, Tn, mu, rho, tau, sigma,
                               NREP_oracle = 1e7, batch_oracle = 2e5,
                               cores = n.cores, mu_trend) {
  out <- vapply(
    gammas,
    function(g) compute_tauRD(
      c0 = c0, grid = grid, Tn = Tn, mu = mu, 
      rho = rho, tau = tau, sigma = sigma,
      disc_factor = g, 
      NREP_oracle = NREP_oracle, batch_oracle = batch_oracle, cores = cores,
      mu_trend = mu_trend
    ),
    numeric(1)
  )
  names(out) <- as.character(gammas)
  out
}

# ---------- one gamma block ----------
.run_one_gamma <- function(gamma, n_grid, reps, c0, h, Tn, 
                           mu, sigma, rho, tau,
                           tauRD_map, 
                           cores = n.cores,
                           K = K_uniform,
                           time_fe = FALSE,
                           time_fe_Z = FALSE,
                           time_fe_ZxW = FALSE,
                           mu_trend) {
  tauRD <- unname(tauRD_map[as.character(gamma)])
  tau_partial <- tau * (c0 - mu)
  
  cat("\n=====================================================\n")
  cat(sprintf(" gamma = %.1f   tau_RD ≈ %.3f   tau_partial = %.3f\n",
              gamma, tauRD, tau_partial))
  cat("-----------------------------------------------------\n")
  cat("    n   | proposed τ_RD  | naive τ_RD     | std LLR τ_RD   | std LLR τ_partial\n")
  cat("        | cov     width  | cov    width   | cov    width   | cov      width\n")
  cat("--------+----------------+----------------+----------------+------------------\n")
  
  rows <- vector("list", length(n_grid))
  for (i in seq_along(n_grid)) {
    n <- n_grid[i]
    Rlist <- bench_run(
      n = n, Tn = Tn, reps = reps, thresh = c0, h = h,
      mu = mu, sigma = sigma, rho = rho, tau = tau,
      disc_factor = gamma, weights_by_time = TRUE, K = K,
      mc.cores = cores, time_fe = time_fe,
      time_fe_Z = time_fe_Z, time_fe_ZxW = time_fe_ZxW,
      mu_trend = mu_trend
    )
    
    R <- do.call(rbind, Rlist)
    
    failed <- R[, "failed"] == 1
    n_failed <- sum(failed)
    
    # Keep only successfully completed replications for coverage/width summaries.
    R_ok <- R[!failed, , drop = FALSE]
    
    se_stat  <- R_ok[, "static_se"]
    se_naive <- R_ok[, "naive_se"]
    se_dyn   <- R_ok[, "dynamic_se"]
    
    ci_stat  <- 1.96 * se_stat
    ci_naive <- 1.96 * se_naive
    ci_dyn   <- 1.96 * se_dyn
    
    cov_dyn         <- mean(abs(R_ok[, "dynamic_est"] - tauRD)       <= ci_dyn,   na.rm = TRUE)
    cov_naive_tauRD <- mean(abs(R_ok[, "naive_est"]   - tauRD)       <= ci_naive, na.rm = TRUE)
    cov_stat_tauRD  <- mean(abs(R_ok[, "static_est"]  - tauRD)       <= ci_stat,  na.rm = TRUE)
    cov_stat_tauPE  <- mean(abs(R_ok[, "static_est"]  - tau_partial) <= ci_stat,  na.rm = TRUE)
    
    wid_dyn   <- mean(2 * 1.96 * se_dyn,   na.rm = TRUE)
    wid_naive <- mean(2 * 1.96 * se_naive, na.rm = TRUE)
    wid_stat  <- mean(2 * 1.96 * se_stat,  na.rm = TRUE)
    
    med_dyn   <- median(2 * 1.96 * se_dyn,   na.rm = TRUE)
    med_naive <- median(2 * 1.96 * se_naive, na.rm = TRUE)
    med_stat  <- median(2 * 1.96 * se_stat,  na.rm = TRUE)
    
    n_na <- sum(!is.finite(se_dyn) | !is.finite(se_naive)) + n_failed
    
    cat(.format_row(n, cov_dyn, wid_dyn, cov_naive_tauRD, wid_naive,
                    cov_stat_tauRD, wid_stat, cov_stat_tauPE), "\n")
    cat(sprintf("        | med=%6.3f | med=%6.3f | med=%6.3f |  (NA reps: %d)\n",
                med_dyn, med_naive, med_stat, n_na))
    flush.console()
    
    rows[[i]] <- data.frame(
      n = n,
      proposed_cov = cov_dyn,
      proposed_width = wid_dyn,
      proposed_med_width = med_dyn,
      naive_cov = cov_naive_tauRD,
      naive_width = wid_naive,
      naive_med_width = med_naive,
      stdLLR_cov_tauRD = cov_stat_tauRD,
      stdLLR_width = wid_stat,
      stdLLR_med_width = med_stat,
      stdLLR_cov_tauPE = cov_stat_tauPE,
      n_na = n_na,
      n_failed = n_failed,
      n_success = nrow(R_ok),
      tauRD = tauRD,
      tau_partial = tau_partial,
      gamma = gamma,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---------- master wrapper ----------
run_all_tables <- function(
    gammas = c(0.5, 0.8, 1.0),
    n_grid = 1000*c(1, 2, 4, 8, 16, 32, 64, 128),
    reps = 2000,
    c0 = 110, h = NULL,
    Tn = 12, mu = 100, sigma = 4, rho = 0.9, tau = 0.1,
    K = K_uniform, 
    time_fe = TRUE,
    time_fe_Z = FALSE, 
    time_fe_ZxW = FALSE,
    NREP_oracle = 5e6, 
    batch_oracle = 2e5, grid_span = 5,
    cores = n.cores, mu_trend = 1,
    save_dir = "./tables",
    save_tables = TRUE
) {
  if (save_tables) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  grid <- seq(c0 - grid_span, c0 + grid_span, by = 1)
  
  tauRD_map <- .precompute_tauRDs(
    gammas, c0, grid, Tn, mu, rho, tau, sigma,
    NREP_oracle = NREP_oracle, batch_oracle = batch_oracle, cores = cores,
    mu_trend = mu_trend
  )
  
  out <- vector("list", length(gammas))
  names(out) <- paste0("gamma_", gammas)
  
  for (g in seq_along(gammas)) {
    gamma_g <- gammas[g]
    
    out[[g]] <- .run_one_gamma(
      gamma = gamma_g, n_grid = n_grid, reps = reps,
      c0 = c0, h = h, Tn = Tn, mu = mu, sigma = sigma, 
      rho = rho, tau = tau,
      tauRD_map = tauRD_map, 
      cores = cores, 
      K = K, 
      time_fe = time_fe,
      time_fe_Z = time_fe_Z,
      time_fe_ZxW = time_fe_ZxW,
      mu_trend = mu_trend
    )
    
    if (save_tables) {
      gamma_tag <- gsub("\\.", "p", as.character(gamma_g))
      
      file_name <- sprintf(
        "table_gamma_%s_mt%d_reps%d.csv",
        gamma_tag, mu_trend, reps
      )
      
      file_path <- file.path(save_dir, file_name)
      
      write.csv(out[[g]], file = file_path, row.names = FALSE)
      
      cat(sprintf("\nSaved table for gamma = %s to: %s\n", gamma_g, file_path))
      flush.console()
    }
  }
  
  if (save_tables) {
    all_tables <- do.call(rbind, out)
    all_file <- file.path(save_dir, sprintf("tables_all_gammas_mt%d_reps%d.csv", mu_trend, reps))
    write.csv(all_tables, file = all_file, row.names = FALSE)
    
    cat(sprintf("\nSaved combined table to: %s\n", all_file))
    flush.console()
  }
  
  return(out)
}

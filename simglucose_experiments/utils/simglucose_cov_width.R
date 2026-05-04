# simglucose_cov_width.R finds coverage and width of confidence intervals for 
# the simglucose experiments in the paper https://arxiv.org/abs/2512.15244
#
# Example usage:
#   source("./utils/simglucose_cov_width.R")
#   out <- run_sim(tau_RDs = taus, dir = "data", c0 = 150, mc.cores = parallel::detectCores())

suppressPackageStartupMessages({
  if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow", quiet = TRUE)
  if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply", quiet = TRUE)
  # If running on a cluster, comment the above and install the packages once
  # If having trouble installing arrow, try nanoparquet and 
  #   replace arrow::read_parquet with nanoparquet::read_parquet
  library(arrow)
  library(pbmcapply)
  library(parallel)
})


# ── Kernels ───────────────────────────────────────────────────────────────────
K_triangular   <- function(u) pmax(1 - abs(u), 0)
K_uniform      <- function(u) 0.5 * (abs(u) <= 1)
K_epanechnikov <- function(u) 0.75 * pmax(1 - u^2, 0)

# ── IK bandwidth ─────────────────────────────────────────────────────────────
IK_bandwidth <- function(Y, X, threshold, kernel=c("triangular","uniform")) {
  if (threshold >= max(X) || threshold <= min(X))
    stop("RD threshold is outside the running variable range.")
  kernel <- match.arg(kernel); x <- X - threshold; n <- length(Y)
  h.sil <- 1.84 * stats::sd(x) * n^(-1/5)
  ip1 <- x >= 0 & x <=  h.sil;  im1 <- x < 0 & x >= -h.sil
  s2  <- (sum((Y[ip1]-mean(Y[ip1]))^2) + sum((Y[im1]-mean(Y[im1]))^2)) /
    (sum(ip1)+sum(im1))
  fc  <- (sum(ip1)+sum(im1)) / (2*n*h.sil)
  mid <- x >= stats::median(x[x<0]) & x <= stats::median(x[x>=0])
  xm  <- x[mid]; ym <- Y[mid]
  tt  <- cbind(1, xm>=0, xm, xm^2, xm^3)
  gv  <- solve(t(tt)%*%tt, t(tt)%*%ym)
  td  <- 6*gv[5]
  hp2 <- 3.56*(s2/(fc*max(td^2,0.01)))^(1/7)*sum(x>=0)^(-1/7)
  hm2 <- 3.56*(s2/(fc*max(td^2,0.01)))^(1/7)*sum(x<0)^(-1/7)
  ip3 <- x>=0 & x<=hp2;  im3 <- x<0 & x>=-hm2
  if (sum(im1)<=2||sum(ip3)<=2) stop("Insufficient obs near discontinuity.")
  fit2 <- function(xi,yi) {
    tt2 <- cbind(1,xi,xi^2); solve(t(tt2)%*%tt2,t(tt2)%*%yi)
  }
  sd2p <- 2*fit2(x[ip3],Y[ip3])[3]; sd2m <- 2*fit2(x[im3],Y[im3])[3]
  rp   <- 2160*s2/(sum(ip3)*hp2^4); rm <- 2160*s2/(sum(im3)*hm2^4)
  h.opt <- 3.4375*((2*s2)/(fc*((sd2p-sd2m)^2+rp+rm)))^(1/5)*n^(-1/5)
  list(bandwidth=as.numeric(h.opt))
}

# ── Design matrix ─────────────────────────────────────────────────────────────
.build_X <- function(Zk, Wk, tk, Tn,
                     groupk = NULL,
                     time_fe=TRUE, time_fe_Z=FALSE, time_fe_ZxW=FALSE,
                     group_fe=TRUE) {
  X <- cbind(Intercept=1.0, Z=Zk, W=Wk, ZxW=Zk*Wk)

  if (!is.null(groupk) && group_fe) {
    G <- model.matrix(~factor(groupk, levels=sort(unique(groupk))))[, -1, drop=FALSE]
    colnames(G) <- paste0("gFE_", sort(unique(groupk))[-1])
    X <- cbind(X, G)
  }

  if (Tn>1 && (time_fe||time_fe_Z||time_fe_ZxW)) {
    D <- model.matrix(~factor(tk,levels=seq_len(Tn)))[,-1,drop=FALSE]
    colnames(D) <- paste0("tFE_",seq_len(Tn-1))
    if (time_fe)     X <- cbind(X, D)
    if (time_fe_Z)   {
      DZ  <- sweep(D,1,Zk,"*")
      colnames(DZ) <- paste0("Z:tFE_",seq_len(Tn-1))
      X <- cbind(X,DZ)
    }
    if (time_fe_ZxW) {
      DZW <- sweep(D,1,Zk*Wk,"*")
      colnames(DZW) <- paste0("ZxW:tFE_",seq_len(Tn-1))
      X <- cbind(X,DZW)
    }
  }
  X
}

# ── Cholesky solver ───────────────────────────────────────────────────────────
chol_solve <- function(A, B=diag(nrow(A))) {
  R <- chol(A); backsolve(R, forwardsolve(t(R), B))
}

# ── Fast WLS + cluster-robust VCV ─────────────────────────────────────────────
.fast_wls_cr <- function(Xw, yw, ck, vcov_type = c("CR0", "CR1", "CR3")) {
  vcov_type <- match.arg(vcov_type)

  A    <- crossprod(Xw)
  Ainv <- chol_solve(A)
  beta <- Ainv %*% crossprod(Xw, yw)

  e <- as.numeric(yw - Xw %*% beta)

  U0 <- rowsum(Xw * e, group = ck, reorder = FALSE)

  if (vcov_type == "CR0") {
    meat <- crossprod(U0)

  } else if (vcov_type == "CR1") {
    meat <- crossprod(U0)

  } else if (vcov_type == "CR3") {
    idx <- split(seq_len(nrow(Xw)), ck)
    G   <- length(idx)
    p   <- ncol(Xw)
    U3  <- matrix(0, nrow = G, ncol = p)

    for (g in seq_along(idx)) {
      ii <- idx[[g]]
      Xg <- Xw[ii, , drop = FALSE]
      eg <- e[ii]

      Mg <- diag(length(ii)) - Xg %*% Ainv %*% t(Xg)

      rg <- tryCatch(
        qr.solve(Mg, eg),
        error = function(err) solve(Mg + 1e-10 * diag(nrow(Mg)), eg)
      )

      U3[g, ] <- drop(crossprod(Xg, rg))
    }

    meat <- crossprod(U3)
  }

  V <- Ainv %*% meat %*% Ainv

  if (vcov_type == "CR1") {
    G <- nrow(U0)
    N <- nrow(Xw)
    p <- ncol(Xw)
    if (G > 1 && N > p) {
      V <- V * (G / (G - 1)) * ((N - 1) / (N - p))
    }
  }

  list(beta = drop(beta), vcov = V, Ainv = Ainv, U = U0)
}

# ── gamma -> column label and names ──────────────────────────────────────────
# Label convention: gamma = 1 - 1/d  ->  label = as.character(d)
#                   gamma = 1         ->  label = "Inf"
# Column names: G_{label}, H_{label}
.gamma_label <- function(g) {
  if (abs(g - 1.0) < 1e-12) "Inf" else as.character(round(1 / (1 - g)))
}

.gamma_cols <- function(g) {
  lab <- .gamma_label(g)
  list(G = paste0("G_", lab), H = paste0("H_", lab), label = lab)
}

# ── RD estimators ─────────────────────────────────────────────────────────────
.delta_method_se <- function(fitG, fitH, iW) {
  bG <- fitG$beta[iW]; bH <- fitH$beta[iW]
  if (!is.finite(bH)||abs(bH)<.Machine$double.eps) stop("Denominator near zero.")
  VG  <- fitG$vcov[iW,iW]; VH <- fitH$vcov[iW,iW]
  Cov <- sum((fitG$U%*%t(fitG$Ainv))[,iW]*(fitH$U%*%t(fitH$Ainv))[,iW])
  tau <- bG/bH
  se  <- sqrt(max((VG+tau^2*VH-2*tau*Cov)/bH^2, 0))
  list(est=tau, se=se, jump_G=bG, jump_H=bH, VG=VG, VH=VH, Cov_GH=Cov)
}

rdd_static_from_df <- function(df, thresh, h, K=K_uniform,
                               time_fe=TRUE, time_fe_Z=FALSE, time_fe_ZxW=FALSE,
                               group_fe=TRUE, vcov_type = c("CR0", "CR1", "CR3")) {
  vcov_type <- match.arg(vcov_type)
  stopifnot(h>0); Tn <- max(df$t)+1
  Z <- df$Z_it-thresh; W <- as.numeric(Z>=0); wK <- K(Z/h); keep <- wK>0
  if (!any(keep)) stop("No obs in bandwidth (static).")
  Zk<-Z[keep]; Wk<-W[keep]; Yk<-df$Y_it[keep]
  ck<-as.integer(factor(df$i[keep])); tk<-df$t[keep]+1
  gk <- df$group[keep]
  X <- .build_X(Zk,Wk,tk,Tn,
                groupk=gk,
                time_fe=time_fe,time_fe_Z=time_fe_Z,time_fe_ZxW=time_fe_ZxW,
                group_fe=group_fe)
  sw <- sqrt(wK[keep]); Xw <- X*sw
  kc <- colSums(abs(Xw))>0
  if (!kc[which(colnames(X)=="W")]) stop("No variation in W (static).")
  Xw<-Xw[,kc,drop=FALSE]; X<-X[,kc,drop=FALSE]; iW<-which(colnames(X)=="W")
  fit <- .fast_wls_cr(Xw, Yk*sw, ck, vcov_type)
  list(est=fit$beta[iW], se=sqrt(max(fit$vcov[iW,iW],0)))
}

rdd_dynamic_from_df <- function(df, gamma, thresh, h,
                                weights_by_time=TRUE,
                                K=K_uniform,
                                time_fe=TRUE,
                                time_fe_Z=FALSE,
                                time_fe_ZxW=FALSE,
                                group_fe=TRUE,
                                vcov_type = c("CR0", "CR1", "CR3")) {
  vcov_type <- match.arg(vcov_type)
  stopifnot(h>0); gc<-.gamma_cols(gamma); Tn<-max(df$t)+1
  if(!(gc$G%in%names(df))) stop("Missing column: ",gc$G)
  if(!(gc$H%in%names(df))) stop("Missing column: ",gc$H)
  Z<-df$Z_it-thresh; W<-as.numeric(Z>=0); wK<-K(Z/h); keep<-wK>0
  if(!any(keep)) stop("No obs in bandwidth (dynamic).")
  Zk<-Z[keep]; Wk<-W[keep]; Gk<-df[[gc$G]][keep]; Hk<-df[[gc$H]][keep]
  ck<-as.integer(factor(df$i[keep])); tk<-df$t[keep]+1
  sw <- sqrt(wK[keep]*(if(weights_by_time) gamma^(tk-1) else 1.0))
  gk <- df$group[keep]
  X <- .build_X(Zk,Wk,tk,Tn,
                groupk=gk,
                time_fe=time_fe,time_fe_Z=time_fe_Z,time_fe_ZxW=time_fe_ZxW,
                group_fe=group_fe)
  Xw <- X*sw; kc <- colSums(abs(Xw))>0
  if(!kc[which(colnames(X)=="W")]) stop("No variation in W (dynamic).")
  Xw<-Xw[,kc,drop=FALSE]; X<-X[,kc,drop=FALSE]; iW<-which(colnames(X)=="W")
  fitG<-.fast_wls_cr(Xw,Gk*sw,ck,vcov_type)
  fitH<-.fast_wls_cr(Xw,Hk*sw,ck,vcov_type)
  .delta_method_se(fitG, fitH, iW)
}

rdd_naive_from_df <- function(df, gamma, thresh, h, K=K_uniform,
                              group_fe=TRUE,
                              vcov_type = c("CR0", "CR1", "CR3")) {
  vcov_type <- match.arg(vcov_type)
  stopifnot(h>0); gc<-.gamma_cols(gamma)
  if(!(gc$G%in%names(df))) stop("Missing column: ",gc$G)
  if(!(gc$H%in%names(df))) stop("Missing column: ",gc$H)
  # First meal step per unit: lowest t with finite Z_it.
  meal_df <- df[is.finite(df$Z_it), , drop=FALSE]
  d0      <- meal_df[!duplicated(meal_df$i), , drop=FALSE]
  Z<-d0$Z_it-thresh; W<-as.numeric(Z>=0); wK<-K(Z/h); keep<-wK>0
  if(!any(keep)) stop("No obs in bandwidth (naive).")
  Zk<-Z[keep]; Wk<-W[keep]; Gk<-d0[[gc$G]][keep]; Hk<-d0[[gc$H]][keep]
  ck<-as.integer(factor(d0$i[keep]))
  gk <- d0$group[keep]
  X <- .build_X(Zk,Wk,rep(1,length(Zk)),Tn=1,
                groupk=gk,
                group_fe=group_fe)
  sw <- sqrt(wK[keep]); Xw <- X*sw; kc <- colSums(abs(Xw))>0
  if(!kc[which(colnames(X)=="W")]) stop("No variation in W (naive).")
  Xw<-Xw[,kc,drop=FALSE]; X<-X[,kc,drop=FALSE]; iW<-which(colnames(X)=="W")
  fitG<-.fast_wls_cr(Xw,Gk*sw,ck,vcov_type)
  fitH<-.fast_wls_cr(Xw,Hk*sw,ck,vcov_type)
  .delta_method_se(fitG, fitH, iW)
}

# ── compute_cis_for_df ────────────────────────────────────────────────────────
compute_cis_for_df <- function(df, 
                               gamma, 
                               c0=150, 
                               alpha=0.05, 
                               K=K_uniform,
                               time_fe=TRUE, time_fe_Z=FALSE, time_fe_ZxW=FALSE,
                               group_fe=TRUE,
                               weights_by_time_dynamic=TRUE,
                               h_static=NULL, h_dynamic=NULL, h_naive=NULL,
                               ik_kernel=c("uniform","triangular"),
                               vcov_type=c("CR0","CR1","CR3")) {
  vcov_type <- match.arg(vcov_type)
  ik_kernel <- match.arg(ik_kernel); zcrit <- qnorm(1-alpha/2)
  # Non-meal rows have Z_it = -Inf; restrict to meal steps for bandwidth selection.
  meal_rows <- is.finite(df$Z_it)
  if (is.null(h_static)) {
    h_static <- IK_bandwidth(Y=df$Y_it[meal_rows], X=df$Z_it[meal_rows],
                             threshold=c0, kernel=ik_kernel)$bandwidth
  }
  if (is.null(h_dynamic)) {
    gc        <- .gamma_cols(gamma)
    h_dynamic <- IK_bandwidth(Y=df[[gc$G]][meal_rows], X=df$Z_it[meal_rows],
                              threshold=c0, kernel=ik_kernel)$bandwidth
  }
  if (is.null(h_naive)) {
    # In the new setting t is the outer-step index (0..191); t=0 is almost never
    # a meal step, so filter to the first finite-Z row per unit instead.
    meal_df <- df[is.finite(df$Z_it), , drop=FALSE]
    d0      <- meal_df[!duplicated(meal_df$i), , drop=FALSE]
    h_naive <- IK_bandwidth(Y=d0[[.gamma_cols(gamma)$G]], X=d0$Z_it,
                            threshold=c0, kernel=ik_kernel)$bandwidth
  }
  stat <- rdd_static_from_df(df,thresh=c0,h=h_static,K=K,
                             time_fe=time_fe,time_fe_Z=time_fe_Z,time_fe_ZxW=time_fe_ZxW,
                             group_fe=group_fe, vcov_type=vcov_type)
  dyn  <- rdd_dynamic_from_df(df,gamma=gamma,thresh=c0,h=h_dynamic,
                              weights_by_time=weights_by_time_dynamic,K=K,
                              time_fe=time_fe,time_fe_Z=time_fe_Z,time_fe_ZxW=time_fe_ZxW,
                              group_fe=group_fe, vcov_type=vcov_type)
  nai  <- rdd_naive_from_df(df,gamma=gamma,thresh=c0,h=h_naive,K=K,
                            group_fe=group_fe, vcov_type=vcov_type)
  data.frame(
    gamma=gamma, c0=c0,
    h_static=h_static,
    h_dynamic=h_dynamic,
    h_naive=h_naive,
    static_est=stat$est,   static_se=stat$se,
    static_ci_l=stat$est   - zcrit*stat$se, static_ci_u=stat$est   + zcrit*stat$se,
    dynamic_est=dyn$est,   dynamic_se=dyn$se,
    dynamic_ci_l=dyn$est   - zcrit*dyn$se,  dynamic_ci_u=dyn$est   + zcrit*dyn$se,
    naive_est=nai$est,     naive_se=nai$se,
    naive_ci_l=nai$est     - zcrit*nai$se,  naive_ci_u=nai$est     + zcrit*nai$se,
    row.names=NULL
  )
}

# ── compute_value_columns ─────────────────────────────────────────────────────
# For each gamma, computes backward discounted sums:
#   G_t = Y_t + gamma * Y_{t+1} + ...  (df$Y_it[t] = instantaneous reward at outer step t)
#   H_t = A_t + gamma * A_{t+1} + ...  (A_t = 0 at non-meal steps)
# Column names: G_{label} and H_{label} where label = round(1/(1-gamma))
# or "Inf" for gamma=1.
compute_value_columns <- function(df, gammas = 1 - 1/c(15, 30, 60, 120, Inf)) {
  df <- df[order(df$i, df$t),,drop=FALSE]
  idx_list <- split(seq_len(nrow(df)), df$i)
  for (gam in gammas) {
    lab <- .gamma_label(gam)
    G_all <- numeric(nrow(df)); H_all <- numeric(nrow(df))
    for (rows in idx_list) {
      Yi <- df$Y_it[rows]; Ai <- df$A_it[rows]; n <- length(rows)
      Gi <- numeric(n); Hi <- numeric(n)
      Gi[n] <- Yi[n]; Hi[n] <- Ai[n]
      if (n>1) for (k in (n-1):1) {
        Gi[k] <- Yi[k] + gam * Gi[k+1]
        Hi[k] <- Ai[k] + gam * Hi[k+1]
      }
      G_all[rows] <- Gi; H_all[rows] <- Hi
      # scale <- if (abs(gam - 1.0) < 1e-12) 1.0 / n else (1.0 - gam)
      # G_all[rows] <- Gi * scale; H_all[rows] <- Hi * scale
    }
    df[[paste0("G_", lab)]] <- G_all
    df[[paste0("H_", lab)]] <- H_all
  }
  df
}

# ── .add_group_from_time_reset ────────────────────────────────────────────────
.add_group_from_time_reset <- function(raw, time_col = "t") {
  tt <- as.integer(raw[[time_col]])
  new_traj <- c(TRUE, diff(tt) < 0)
  traj_id  <- cumsum(new_traj)
  group    <- ((traj_id - 1) %% 10) + 1
  raw$traj_id_reconstructed <- traj_id
  raw$group <- group
  raw
}

# ── .prep_df ──────────────────────────────────────────────────────────────────
.prep_df <- function(raw, gammas = 1 - 1/c(24, 48, 96, 192, Inf)) {
  names(raw)[names(raw)=="Z"] <- "Z_it"
  names(raw)[names(raw)=="A"] <- "A_it"
  names(raw)[names(raw)=="Y"] <- "Y_it"
  raw$i <- as.integer(raw$i); raw$t <- as.integer(raw$t)
  raw <- .add_group_from_time_reset(raw, time_col = "t")
  # overwrite i with globally unique trajectory ID so that compute_value_columns
  # and all rdd_* clustering are correct whether processing one rep or pooled reps
  raw$i <- raw$traj_id_reconstructed
  compute_value_columns(raw, gammas=gammas)
}

# ── Process one rep ───────────────────────────────────────────────────────────
.process_one_rep <- function(path, 
                             tau_RDs, 
                             gammas, 
                             c0, 
                             K,
                             time_fe, time_fe_Z, time_fe_ZxW, group_fe,
                             weights_by_time_dynamic,
                             vcov_type = c("CR0","CR1","CR3"), ...) {
  vcov_type <- match.arg(vcov_type)
  tryCatch({
    df  <- .prep_df(as.data.frame(arrow::read_parquet(path)), gammas=gammas)
    cov_mat <- matrix(NA_real_, nrow=length(gammas), ncol=3,
                      dimnames=list(paste0("gamma=",gammas), c("static","naive","dynamic")))
    wid_mat <- cov_mat
    bw_mat  <- matrix(NA_real_, nrow=length(gammas), ncol=3,
                      dimnames=list(paste0("gamma=",gammas), c("h_static","h_dynamic","h_naive")))
    for (gid in seq_along(gammas)) {
      out <- compute_cis_for_df(df, gamma=gammas[gid], c0=c0, K=K,
                                time_fe=time_fe,
                                group_fe=group_fe,
                                time_fe_Z=time_fe_Z,
                                time_fe_ZxW=time_fe_ZxW,
                                weights_by_time_dynamic=weights_by_time_dynamic,
                                vcov_type=vcov_type, ...)
      tr <- tau_RDs[gid]
      cov_mat[gid,"static"]  <- as.numeric(out$static_ci_l  <= tr & tr <= out$static_ci_u)
      cov_mat[gid,"naive"]   <- as.numeric(out$naive_ci_l   <= tr & tr <= out$naive_ci_u)
      cov_mat[gid,"dynamic"] <- as.numeric(out$dynamic_ci_l <= tr & tr <= out$dynamic_ci_u)
      wid_mat[gid,"static"]  <- out$static_ci_u  - out$static_ci_l
      wid_mat[gid,"naive"]   <- out$naive_ci_u   - out$naive_ci_l
      wid_mat[gid,"dynamic"] <- out$dynamic_ci_u - out$dynamic_ci_l
      bw_mat[gid, "h_static"]  <- out$h_static
      bw_mat[gid, "h_dynamic"] <- out$h_dynamic
      bw_mat[gid, "h_naive"]   <- out$h_naive
    }
    list(ok=TRUE, coverage=cov_mat, width=wid_mat, bw=bw_mat, err=NA_character_)
  }, error=function(e) list(ok=FALSE, coverage=NULL, width=NULL, bw=bw_mat, err=conditionMessage(e)))
}

# ── main computations ─────────────────────────────────────────────────────────
run_sim <- function(tau_RDs,
                    dir      = "data",
                    mc.cores = parallel::detectCores(),
                    gammas   = 1 - 1/c(24, 48, 96, 192, Inf),
                    c0       = 150,
                    K        = K_uniform,
                    time_fe  = TRUE,
                    group_fe = TRUE,
                    time_fe_Z = FALSE,
                    time_fe_ZxW = FALSE,
                    weights_by_time_dynamic = TRUE,
                    vcov_type = "CR3",
                    ...) {
  stopifnot(length(tau_RDs) == length(gammas))
  files <- sort(list.files(dir, pattern="^rep[0-9]{4}\\.parquet$", full.names=TRUE))
  if (length(files)==0) stop("No rep*.parquet files found in: ", dir)
  cat(sprintf("Found %d rep files in '%s'\n", length(files), dir))

  worker <- function(path)
    .process_one_rep(path, tau_RDs=tau_RDs, gammas=gammas, c0=c0, K=K,
                     time_fe=time_fe,
                     time_fe_Z=time_fe_Z,
                     time_fe_ZxW=time_fe_ZxW,
                     group_fe=group_fe,
                     weights_by_time_dynamic=weights_by_time_dynamic,
                     vcov_type=vcov_type, ...)

  # results <- pbmcapply::pbmclapply(
  #   files,
  #   worker,
  #   mc.cores = mc.cores,
  #   mc.preschedule = TRUE
  # )
  # cat("\n")
  
  results <- if (mc.cores > 1L) {
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(arrow))
      NULL
    })
    
    parallel::clusterExport(
      cl,
      varlist = c(
        ".add_group_from_time_reset",
        ".process_one_rep", ".prep_df", "compute_value_columns",
        "compute_cis_for_df", "rdd_static_from_df", "rdd_dynamic_from_df",
        "rdd_naive_from_df", ".gamma_cols", ".gamma_label", ".delta_method_se",
        ".fast_wls_cr", "chol_solve", ".build_X", "IK_bandwidth",
        "K_uniform", "K_triangular", "K_epanechnikov",
        "tau_RDs", "gammas", "c0", "K",
        "time_fe", "time_fe_Z","time_fe_ZxW", "group_fe",
        "weights_by_time_dynamic", "vcov_type"
      ),
      envir = environment()
    )
    
    parallel::parLapply(cl, files, worker)
  } else {
    lapply(seq_along(files), function(i) {
      cat(sprintf("  rep %d/%d\r", i, length(files)))
      worker(files[i])
    })
  }
  cat("\n")

  ok <- vapply(results, function(x) is.list(x) && length(x$ok) == 1 && isTRUE(x$ok), logical(1))
  good <- results[ok]; bad <- results[!ok]
  if (length(good)==0) stop("All reps failed. Example: ", bad[[1]]$err)

  cov_arr <- simplify2array(lapply(good, `[[`, "coverage"))
  wid_arr <- simplify2array(lapply(good, `[[`, "width"))

  avg_cov <- apply(cov_arr, c(1,2), mean,   na.rm=TRUE)
  avg_wid <- apply(wid_arr, c(1,2), mean,   na.rm=TRUE)
  med_wid <- apply(wid_arr, c(1,2), median, na.rm=TRUE)

  cat(sprintf("\nReps used: %d  |  Failed: %d\n", length(good), length(bad)))
  if (length(bad)>0) for (i in seq_along(bad))
    cat("  FAILED:", basename(files[which(!ok)[i]]), "->", bad[[i]]$err, "\n")

  cat("\nEmpirical coverage:\n");  print(round(avg_cov,3))
  cat("\nMean CI width:\n");       print(round(avg_wid,4))
  cat("\nMedian CI width:\n");     print(round(med_wid,4))
  
  bw_arr <- simplify2array(lapply(good, `[[`, "bw"))   # [gammas x 2 x reps]
  
  avg_bw <- apply(bw_arr, c(1,2), mean,   na.rm=TRUE)
  min_bw <- apply(bw_arr, c(1,2), min,    na.rm=TRUE)
  max_bw <- apply(bw_arr, c(1,2), max,    na.rm=TRUE)
  
  cat("\nMean bandwidth (h_static, h_dynamic, h_naive):\n")
  print(round(avg_bw, 4))
  # cat("\nBandwidth range [min, max]:\n")
  # bw_rng <- array(NA_real_, dim=c(nrow(avg_bw), 3, 2),
  #                 dimnames=list(rownames(avg_bw), c("h_static","h_dynamic","h_naive"), c("min","max")))
  # bw_rng[,,"min"] <- min_bw
  # bw_rng[,,"max"] <- max_bw
  # for (col in c("h_static","h_dynamic","h_naive")) {
  #   cat(sprintf("  %s:\n", col))
  #   m <- cbind(min=min_bw[,col], max=max_bw[,col])
  #   rownames(m) <- rownames(avg_bw)
  #   print(round(m, 4))
  # }
  
  invisible(list(avg_coverage=avg_cov, avg_width=avg_wid, med_width=med_wid,
                 avg_bw=avg_bw,
                 n_reps_used=length(good), n_reps_failed=length(bad),
                 failures=if(length(bad)>0) files[which(!ok)] else NULL))
}

# ── CLI ───────────────────────────────────────────────────────────────────────
.parse_args <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  out  <- list(tau="", dir="data", cores=parallel::detectCores(), c0=150, save="")
  for (a in args) {
    kv <- strsplit(a,"=",fixed=TRUE)[[1]]; if (length(kv)!=2) next
    k<-kv[1]; v<-kv[2]
    if(k=="--tau")   out$tau   <- v
    if(k=="--dir")   out$dir   <- v
    if(k=="--cores") out$cores <- as.integer(v)
    if(k=="--c0")    out$c0    <- as.numeric(v)
    if(k=="--save")  out$save  <- v
  }
  out
}

main <- function() {
  opt <- .parse_args()
  if (!nzchar(opt$tau)) stop("Provide --tau=v1,v2,v3,v4,v5")
  tau_vec <- as.numeric(strsplit(opt$tau,",")[[1]])
  gammas  <- 1 - 1/c(24, 48, 96, 192, Inf)
  if (length(tau_vec) != length(gammas) || any(!is.finite(tau_vec)))
    stop("--tau needs 5 finite numbers matching gammas 1-1/24, 1-1/48, 1-1/96, 1-1/192, 1")
  names(tau_vec) <- paste0("gamma=", round(gammas, 6))
  cat("Oracle tau_RDs:\n"); print(tau_vec)
  out <- run_sim(tau_RDs=tau_vec, dir=opt$dir, mc.cores=opt$cores, c0=opt$c0)
  if (nzchar(opt$save)) { saveRDS(out,file=opt$save); cat("Saved to:",opt$save,"\n") }
  invisible(out)
}
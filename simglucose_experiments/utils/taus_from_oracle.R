pkgs <- c("dplyr", "tidyr", "ggplot2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

suppressPackageStartupMessages(library(grid))

compute_tau <- function(dir = "./oracle", 
                        n_seeds = 10000,
                        thresholds = c(145, 150, 155),
                        eval_at = 150) {
  # Read all threshold files and stack
  all_data <- do.call(rbind, lapply(thresholds, function(cv) {
    f <- file.path(dir, paste0("oracle_c", cv, "_", n_seeds, "seeds.csv"))
    if (!file.exists(f)) stop("Missing: ", f)
    df <- read.csv(f); df$c <- cv; df
  }))
  
  # For each (gamma, days): fit natural cubic spline, evaluate derivative at eval_at
  keys <- unique(all_data[, c("gamma", "days")])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    g <- keys$gamma[i]; d <- keys$days[i]
    sub <- all_data[abs(all_data$gamma - g) < 1e-9 & all_data$days == d, ]
    sub <- sub[order(sub$c), ]
    dVY <- splinefun(sub$c, sub$VY, method = "natural")(eval_at, deriv = 1)
    dVA <- splinefun(sub$c, sub$VA, method = "natural")(eval_at, deriv = 1)
    data.frame(gamma = g, days = d, dVY = dVY, dVA = dVA, tau = dVY / dVA)
  }))
  
  out <- out[order(out$gamma, out$days), ]
  print(out, digits = 6, row.names = FALSE)
  invisible(list(tau = out, all_data = all_data))
}

# -------------------------------------------------------------------
# Helpers for plots
# -------------------------------------------------------------------


plot_oracle <- function(result, 
                        eval_at = 150, 
                        all.plots = F) {
  all_data  <- result$all_data
  gammas    <- sort(unique(all_data$gamma))
  days_vals <- sort(unique(all_data$days))
  
  base_cols <- c("orange", "#E41A1C", "#377EB8", "#4DAF4A")
  cols <- setNames(base_cols[seq_along(days_vals)], as.character(days_vals))
  
  tau_all <- list()
  k <- 1
  
  par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
  
  for (g in gammas) {
    g_lbl <- if (abs(g) < 1e-12) {
      "γ = 0"
    } else if (abs(g - 1) < 1e-12) {
      "γ = 1"
    } else {
      paste0("γ = 1-1/", round(1 / (1 - g)))
    }
    
    sub_g <- all_data[abs(all_data$gamma - g) < 1e-9, ]
    
    tau_store <- vector("list", length(days_vals))
    names(tau_store) <- as.character(days_vals)
    tau_range <- c(Inf, -Inf)
    
    for (d in days_vals) {
      sub_d <- sub_g[sub_g$days == d, ]
      sub_d <- sub_d[order(sub_d$c), ]
      
      c_mid <- sub_d$c[2:(nrow(sub_d) - 1)]
      tau   <- (sub_d$VY[3:nrow(sub_d)] - sub_d$VY[1:(nrow(sub_d) - 2)]) /
        (sub_d$VA[3:nrow(sub_d)] - sub_d$VA[1:(nrow(sub_d) - 2)])
      
      keep <- is.finite(tau)
      if (any(keep)) {
        tau_range[1] <- min(tau_range[1], min(tau[keep]))
        tau_range[2] <- max(tau_range[2], max(tau[keep]))
      }
      
      tau_store[[as.character(d)]] <- list(c_mid = c_mid, tau = tau)
      tau_all[[k]] <- data.frame(gamma = g, days = d, c_mid = c_mid, tau = tau)
      k <- k + 1
    }
    
    if(all.plots){
      for (yvar in c("VY", "VA")) {
        plot(NA,
             xlim = range(all_data$c, na.rm = TRUE),
             ylim = range(sub_g[[yvar]], na.rm = TRUE),
             xlab = "threshold c  (mg/dL)",
             ylab = yvar,
             main = paste(g_lbl, "—", yvar))
        abline(v = eval_at, lty = 3, col = "grey60")
        
        for (d in days_vals) {
          sub_d <- sub_g[sub_g$days == d, ]
          sub_d <- sub_d[order(sub_d$c), ]
          
          xs <- seq(min(sub_d$c), max(sub_d$c), length.out = 200)
          ys <- splinefun(sub_d$c, sub_d[[yvar]], method = "natural")(xs)
          
          lines(xs, ys, col = cols[as.character(d)], lwd = 2)
          points(sub_d$c, sub_d[[yvar]],
                 col = cols[as.character(d)], pch = 19, cex = 1.2)
        }
        
        legend("right",
               legend = paste(days_vals, "days"),
               col = cols[as.character(days_vals)],
               lwd = 2, pch = 19, bty = "n", cex = 0.85)
      }
      
      if (!all(is.finite(tau_range))) tau_range <- c(-1, 1)
      
      plot(NA,
           xlim = range(all_data$c, na.rm = TRUE),
           ylim = tau_range,
           xlab = "threshold c  (mg/dL)",
           ylab = expression(tau[RD](c)),
           main = paste(g_lbl, "—", expression(tau[RD](c))))
      abline(v = eval_at, lty = 3, col = "grey60")
      abline(h = 0, lty = 2, col = "grey70")
      
      for (d in days_vals) {
        obj <- tau_store[[as.character(d)]]
        keep <- is.finite(obj$tau)
        lines(obj$c_mid[keep], obj$tau[keep],
              col = cols[as.character(d)], lwd = 2)
        points(obj$c_mid[keep], obj$tau[keep],
               col = cols[as.character(d)], pch = 19, cex = 1.2)
      }
      
      legend("right",
             legend = paste(days_vals, "days"),
             col = cols[as.character(days_vals)],
             lwd = 2, pch = 19, bty = "n", cex = 0.85)
    }
  }
  
  # final combined plot: only days = max(days)
  tau_all <- do.call(rbind, tau_all)
  tau_full <- tau_all[tau_all$days == max(tau_all$days), ]
  
  par(mfrow = c(1, 1), mar = c(4.1, 4.1, 3.1, 2.1))
  
  # soft red for gamma = 0; green/blue gradient for positive gammas
  gamma0_col <- "tomato"
  pos_gammas <- gammas[gammas > 0]
  pos_cols <- if (length(pos_gammas) > 0) {
    grDevices::colorRampPalette(c("#B8E186", "#7BCCC4", "#43A2CA", "#0868AC"))(length(pos_gammas))
  } else {
    character(0)
  }
  
  gamma_cols <- setNames(character(length(gammas)), as.character(gammas))
  for (g in gammas) {
    if (abs(g) < 1e-12) {
      gamma_cols[as.character(g)] <- gamma0_col
    } else {
      gamma_cols[as.character(g)] <- pos_cols[which(pos_gammas == g)]
    }
  }
  
  ylim_tau <- range(tau_full$tau[is.finite(tau_full$tau)], na.rm = TRUE)
  
  plot(NA,
       xlim = range(tau_full$c_mid, na.rm = TRUE),
       ylim = ylim_tau,
       xlab = "Threshold c (mg/dL)",
       ylab = expression(tau[RD](c)),
       main = expression("Dynamic vs Static RD parameters at different thresholds (oracle)"))
  abline(v = eval_at, lty = 3, col = "grey60")
  abline(h = 0, lty = 2, col = "grey70")
  
  for (g in gammas) {
    sub <- tau_full[abs(tau_full$gamma - g) < 1e-9, ]
    sub <- sub[order(sub$c_mid), ]
    keep <- is.finite(sub$tau)
    
    lines(sub$c_mid[keep], sub$tau[keep],
          col = gamma_cols[as.character(g)],
          lwd = 2)
    points(sub$c_mid[keep], sub$tau[keep],
           col = gamma_cols[as.character(g)],
           pch = 19, cex = 1.2)
  }
  
  gamma_labels <- as.expression(lapply(gammas, function(g) {
    if (abs(g) < 1e-12) {
      bquote(gamma == 0)
    } else if (abs(g - 1) < 1e-12) {
      bquote(gamma == 1)
    } else {
      m <- round(1 / (1 - g), 2)
      days_eff <- round(m / 48, 2)   # 48 outer steps per day
      bquote(gamma == 1 - 1/.(m) ~ "(eff. " * .(days_eff) ~ "days" * ")")
    }
  }))
  
  legend("topleft",
         legend = gamma_labels,
         col = gamma_cols[as.character(gammas)],
         lwd = 2,
         pch = 19,
         bty = "n",
         cex = 0.95)
}

.oracle_quantity_levels <- c("minus_dVY", "minus_dVA", "tau")

# plotmath labels for facet strips
.oracle_quantity_labels <- c(
  minus_dVY = "-partialdiff*V*'('*pi[c]*')'/partialdiff*c",
  minus_dVA = "-partialdiff*V^A*'('*pi[c]*')'/partialdiff*c",
  tau       = "tau[plain(RD)]"
)

.oracle_base_theme <- function(base_size = 15) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.position = "bottom",
      axis.title.x = element_text(size = 18),
      axis.text = element_text(size = 13)
    )
}

.gamma_expr <- function(x) {
  x0 <- gsub("\\s+", "", x)
  
  if (x0 %in% c("0", "0.0")) return("gamma == 0")
  if (x0 %in% c("1", "1.0")) return("gamma == 1")
  
  if (grepl("^1-1/[0-9.]+$", x0)) {
    denom <- sub("^1-1/", "", x0)
    return(paste0("gamma == 1 - frac(1, ", denom, ")"))
  }
  
  paste0("gamma == ", x0)
}

.make_gamma_legend_labels_inline <- function(gamma_text, eff_days_text) {
  stopifnot(length(gamma_text) == length(eff_days_text))
  
  out <- character(length(gamma_text))
  
  for (i in seq_along(gamma_text)) {
    g_part <- .gamma_expr(gamma_text[i])
    
    d_raw <- eff_days_text[i]
    d0 <- gsub("\\s+", "", d_raw)
    
    if (tolower(d0) %in% c("inf", "infty", "infinity")) {
      day_part <- "infinity*' days'"
    } else {
      d_num <- suppressWarnings(as.numeric(d0))
      
      if (!is.na(d_num)) {
        unit <- if (d_num < 2) " day" else " days"
        day_part <- paste0("'", format(d_num, trim = TRUE, scientific = FALSE), unit, "'")
      } else {
        day_part <- paste0("'", d_raw, "'")
      }
    }
    
    out[i] <- paste0(g_part, "*' ('*", day_part, "*')'")
  }
  
  parse(text = out)
}

.eff_expr <- function(x) {
  x0 <- gsub("\\s+", "", x)
  
  if (x0 == "") return("phantom(0)")
  if (tolower(x0) %in% c("inf", "infty", "infinity")) return("infinity")
  
  x0
}

.make_gamma_legend_labels <- function(gamma_text, eff_days_text) {
  gamma_part <- vapply(gamma_text, .gamma_expr, character(1))
  eff_part   <- vapply(eff_days_text, .eff_expr, character(1))
  
  if (all(eff_days_text == "")) {
    parse(text = gamma_part)
  } else {
    parse(text = paste0("atop(", gamma_part, ", ", eff_part, ")"))
  }
}

# pick a moderate number of marker locations so the plot is clearer
.make_point_df <- function(df, group_vars, target_n = 14) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      .idx = row_number(),
      .n = n(),
      .step = pmax(floor(.n / target_n), 1L)
    ) %>%
    filter(.idx %% .step == 1L | .idx == .n) %>%
    ungroup() %>%
    select(-.idx, -.n, -.step)
}

# colors for plot 2:
# first is tomato-ish for gamma = 0,
# all remaining are shades of blue from light to dark
.make_plot2_colors <- function(n) {
  stopifnot(n >= 1)
  
  tomato_col <- "#E76F51"
  
  if (n == 1) {
    return(tomato_col)
  }
  
  blue_fun <- grDevices::colorRampPalette(c("#BDD7EE", "#08519C"))
  c(tomato_col, blue_fun(n - 1))
}

# -------------------------------------------------------------------
# Plot 1: fixed gamma, varying days
# -------------------------------------------------------------------

plot_oracle_one_gamma <- function(result,
                                  gamma_target = 1 - 1/96,
                                  spar = 0.55,
                                  n_grid = 400,
                                  colors = c(
                                    "#E69F00",  # orange
                                    "#56B4E9",  # sky blue
                                    "#009E73",  # bluish green
                                    "#CC79A7"   # reddish purple
                                  ),
                                  linetypes = c(
                                    "solid",
                                    "longdash",
                                    "dashed",
                                    "dotdash"
                                  ),
                                  shapes = c(16, 17, 15, 18),
                                  line_width = 0.90,
                                  point_size = 1.7,
                                  point_n = 14) {
  tol <- 1e-8
  
  deriv_df <- result$all_data %>%
    filter(abs(gamma - gamma_target) < tol) %>%
    arrange(days, c) %>%
    group_by(days) %>%
    group_modify(~ {
      x <- .x$c
      grid <- seq(min(x), max(x), length.out = n_grid)
      
      fit_y <- smooth.spline(x, .x$VY, spar = spar)
      fit_a <- smooth.spline(x, .x$VA, spar = spar)
      
      dvy <- predict(fit_y, x = grid, deriv = 1)$y
      dva <- predict(fit_a, x = grid, deriv = 1)$y
      
      tibble(
        c = grid,
        minus_dVY = -dvy,
        minus_dVA = -dva,
        tau = dvy / dva
      )
    }) %>%
    ungroup()
  
  days_levels <- sort(unique(deriv_df$days))
  stopifnot(length(colors) >= length(days_levels))
  stopifnot(length(linetypes) >= length(days_levels))
  stopifnot(length(shapes) >= length(days_levels))
  
  days_chr <- as.character(days_levels)
  
  k <- round(1 / (1 - gamma_target))
  
  # Dummy legend entry placed at the end.
  discount_key <- ".discount_factor"
  legend_breaks <- c(days_chr, discount_key)
  
  legend_labels <- parse(text = c(
    days_chr,
    paste0(
      "paste('Discount factor (', gamma, ') = ', 1 - frac(1, ",
      k,
      "))"
    )
  ))
  
  color_map <- setNames(
    c(colors[seq_along(days_levels)], "#00000000"),
    legend_breaks
  )
  
  linetype_map <- setNames(
    c(linetypes[seq_along(days_levels)], "blank"),
    legend_breaks
  )
  
  shape_map <- setNames(
    c(shapes[seq_along(days_levels)], 32),  # 32 = invisible point
    legend_breaks
  )
  
  plot_df <- deriv_df %>%
    mutate(days = factor(as.character(days), levels = legend_breaks)) %>%
    pivot_longer(
      cols = all_of(.oracle_quantity_levels),
      names_to = "quantity",
      values_to = "value"
    ) %>%
    mutate(
      quantity = factor(
        quantity,
        levels = .oracle_quantity_levels,
        labels = .oracle_quantity_labels
      )
    )
  
  point_df <- .make_point_df(
    plot_df,
    group_vars = c("quantity", "days"),
    target_n = point_n
  )
  
  ggplot(
    plot_df,
    aes(x = c, y = value, color = days, linetype = days, shape = days, group = days)
  ) +
    geom_line(linewidth = line_width, lineend = "round") +
    geom_point(
      data = point_df,
      size = point_size,
      stroke = 0.2
    ) +
    facet_wrap(
      ~ quantity,
      scales = "free_y",
      nrow = 1,
      labeller = label_parsed
    ) +
    scale_color_manual(
      values = color_map,
      breaks = legend_breaks,
      labels = legend_labels,
      name = "Days observed",
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = linetype_map,
      breaks = legend_breaks,
      labels = legend_labels,
      name = "Days observed",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = shape_map,
      breaks = legend_breaks,
      labels = legend_labels,
      name = "Days observed",
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    .oracle_base_theme(base_size = 15) +
    theme(
      legend.title = element_text(size = 14, face = "plain"),
      legend.text  = element_text(size = 14, face = "plain"),
      legend.key.width = unit(1.9, "lines"),
      legend.spacing.x = unit(0.65, "lines")
    )
}

# -------------------------------------------------------------------
# Plot 2: fixed max days, varying gamma
# -------------------------------------------------------------------

plot_oracle_days_max <- function(result,
                                 spar = 0.5,
                                 n_grid = 50,
                                 gamma_values = sort(unique(result$all_data$gamma[result$all_data$days == days_max])),
                                 gamma_text = NULL,
                                 eff_days_text = NULL,
                                 colors = NULL,
                                 line_width = 0.95,
                                 tol = 1e-8) {
  days_max = max(result$all_data$days)
  if (is.null(gamma_text)) {
    gamma_text <- format(gamma_values, digits = 4, trim = TRUE)
  }
  if (is.null(eff_days_text)) {
    eff_days_text <- rep("", length(gamma_values))
  }
  
  stopifnot(length(gamma_values) == length(gamma_text))
  stopifnot(length(gamma_values) == length(eff_days_text))
  
  if (is.null(colors)) {
    colors <- .make_plot2_colors(length(gamma_values))
  }
  
  stopifnot(length(colors) >= length(gamma_values))
  
  gamma_keys <- paste0("g", seq_along(gamma_values))
  color_map <- setNames(colors[seq_along(gamma_values)], gamma_keys)
  
  legend_labels <- .make_gamma_legend_labels_inline(
    gamma_text = gamma_text,
    eff_days_text = eff_days_text
  )
  
  deriv_df <- result$all_data %>%
    filter(days == days_max) %>%
    mutate(
      gamma_idx = vapply(
        gamma,
        function(g) {
          j <- which.min(abs(gamma_values - g))
          if (length(j) == 0 || abs(gamma_values[j] - g) > tol) {
            NA_integer_
          } else {
            j
          }
        },
        integer(1)
      )
    ) %>%
    filter(!is.na(gamma_idx)) %>%
    mutate(
      gamma_f = factor(gamma_keys[gamma_idx], levels = gamma_keys)
    ) %>%
    arrange(gamma_idx, c) %>%
    group_by(gamma_idx, gamma_f) %>%
    group_modify(~ {
      x <- .x$c
      grid <- seq(min(x), max(x), length.out = n_grid)
      
      fit_y <- smooth.spline(x, .x$VY, spar = spar)
      fit_a <- smooth.spline(x, .x$VA, spar = spar)
      
      dvy <- predict(fit_y, x = grid, deriv = 1)$y
      dva <- predict(fit_a, x = grid, deriv = 1)$y
      
      tibble(
        c = grid,
        minus_dVY = -dvy,
        minus_dVA = -dva,
        tau = dvy / dva
      )
    }) %>%
    ungroup()
  
  plot_df <- deriv_df %>%
    pivot_longer(
      cols = all_of(.oracle_quantity_levels),
      names_to = "quantity",
      values_to = "value"
    ) %>%
    mutate(
      quantity = factor(
        quantity,
        levels = .oracle_quantity_levels,
        labels = .oracle_quantity_labels
      )
    )
  
  ggplot(plot_df, aes(x = c, y = value, color = gamma_f, group = gamma_f)) +
    geom_line(linewidth = line_width, lineend = "round") +
    facet_wrap(
      ~ quantity,
      scales = "free_y",
      nrow = 1,
      labeller = label_parsed
    ) +
    scale_color_manual(
      values = color_map,
      breaks = gamma_keys,
      labels = legend_labels,
      name = NULL
    ) +
    labs(
      x = NULL,#"Threshold c"
      y = NULL#,
      #title = paste0("Oracle derivative curves for days = ", days_max)
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    .oracle_base_theme(base_size = 15) +
    theme(
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12, lineheight = 0.95),
      legend.key.width = unit(1.8, "lines"),
      legend.spacing.x = unit(0.7, "lines")
    )
}

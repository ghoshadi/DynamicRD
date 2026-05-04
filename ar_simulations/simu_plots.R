library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)

# ================================================================
# Load table values from ./tables
# ================================================================

TABLE_DIR <- "./tables"

.parse_dynamic_rd_file_info <- function(path) {
  file <- basename(path)
  mt_match <- regmatches(file, regexpr("mt[0-9]+", file))
  gamma_match <- regmatches(file, regexpr("gamma_[^_]+", file))
  
  if (length(mt_match) == 0L || identical(mt_match, character(0))) {
    stop("Could not parse mt from file name: ", file, call. = FALSE)
  }
  if (length(gamma_match) == 0L || identical(gamma_match, character(0))) {
    stop("Could not parse gamma from file name: ", file, call. = FALSE)
  }
  
  mt <- as.integer(sub("^mt", "", mt_match))
  gamma_from_file <- as.numeric(gsub("p", ".", sub("^gamma_", "", gamma_match), fixed = TRUE))
  
  if (is.na(mt)) {
    stop("Parsed mt is NA for file: ", file, call. = FALSE)
  }
  if (is.na(gamma_from_file)) {
    stop("Parsed gamma is NA for file: ", file, call. = FALSE)
  }
  
  list(
    mt = mt,
    gamma = gamma_from_file,
    setting = paste0("Setting ", mt + 1L),
    delta = mt
  )
}

.as_coverage_probability <- function(x) {
  # The current CSV files store coverage as probabilities, e.g. 0.95.
  x <- as.numeric(x)
  ifelse(!is.na(x) & x > 1, x / 100, x)
}

.gamma_plot_label <- function(gamma) {
  paste0("gamma == ", format(gamma, trim = TRUE, scientific = FALSE))
}

load_dynamic_rd_tables <- function(table_dir = TABLE_DIR,
                                   pattern = "^table_gamma_.*_mt[0-9]+_.*\\.csv$") {
  files <- list.files(table_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0L) {
    stop(
      "No table CSV files found in '", table_dir, "'. Expected files like ",
      "table_gamma_0p5_mt0_reps2000.csv.",
      call. = FALSE
    )
  }
  
  required_cols <- c(
    "n",
    "proposed_cov", "proposed_med_width",
    "naive_cov", "naive_med_width",
    "stdLLR_cov_tauRD", "stdLLR_med_width",
    "tauRD", "gamma"
  )
  
  dynamic_rd_table_raw <- lapply(files, function(path) {
    info <- .parse_dynamic_rd_file_info(path)
    df <- read.csv(path, check.names = FALSE)
    
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0L) {
      stop(
        "File '", basename(path), "' is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
    
    df %>%
      transmute(
        setting = info$setting,
        delta = info$delta,
        gamma = as.numeric(gamma),
        gamma_from_file = info$gamma,
        tau_rd = as.numeric(tauRD),
        n = as.integer(n),
        std_rd_coverage = .as_coverage_probability(stdLLR_cov_tauRD),
        std_rd_length = as.numeric(stdLLR_med_width),
        naive_coverage = .as_coverage_probability(naive_cov),
        naive_length = as.numeric(naive_med_width),
        proposed_coverage = .as_coverage_probability(proposed_cov),
        proposed_length = as.numeric(proposed_med_width)
      )
  }) %>%
    bind_rows()
  
  gamma_mismatch <- dynamic_rd_table_raw %>%
    filter(!is.na(gamma), abs(gamma - gamma_from_file) > 1e-10)
  if (nrow(gamma_mismatch) > 0L) {
    warning(
      "Some files have a gamma value in the file name that differs from the gamma column. ",
      "Using the gamma column for plotting.",
      call. = FALSE
    )
  }
  
  dynamic_rd_table_raw %>%
    select(-gamma_from_file) %>%
    arrange(delta, gamma, n)
}

dynamic_rd_table_raw <- load_dynamic_rd_tables(TABLE_DIR)

# ================================================================
# Tidy data for the plot
# Three methods shown in the figure:
# Standard LLR (targeting tau_RD), Naive long-run LLR, Proposed method
# ================================================================

dynamic_rd_plot_data <- dynamic_rd_table_raw %>%
  select(
    setting, delta, gamma, tau_rd, n,
    std_rd_coverage, std_rd_length,
    naive_coverage, naive_length,
    proposed_coverage, proposed_length
  ) %>%
  pivot_longer(
    cols = matches("_(coverage|length)$"),
    names_to = c("method_key", "metric"),
    names_pattern = "(.+)_(coverage|length)",
    values_to = "value"
  ) %>%
  mutate(
    # CSV coverage columns are already probabilities, e.g. 0.95.
    value = as.numeric(value),
    metric = recode(
      metric,
      coverage = "Coverage",
      length   = "Length"
    ),
    metric = factor(metric, levels = c("Coverage", "Length")),
    method = recode(
      method_key,
      std_rd   = "Standard LLR",
      naive    = "Naive long-run LLR",
      proposed = "Proposed method"
    ),
    method = factor(
      method,
      levels = c("Standard LLR", "Naive long-run LLR", "Proposed method")
    ),
    gamma_lab = .gamma_plot_label(gamma),
    gamma_lab = factor(
      gamma_lab,
      levels = .gamma_plot_label(sort(unique(dynamic_rd_table_raw$gamma)))
    )
  )

.dynamic_rd_theme <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
      panel.border = element_blank(),
      
      strip.background = element_blank(),
      strip.text.x = element_text(size = 13, face = "plain"),
      strip.text.y = element_text(size = 13, face = "plain"),
      
      axis.title.x = element_text(size = 14, margin = margin(t = 8)),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 12, color = "grey35"),
      
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.key.width = unit(2.2, "lines"),
      legend.spacing.x = unit(0.55, "lines"),
      
      plot.margin = margin(6, 10, 6, 6)
    )
}

.method_colors <- c(
  "Standard LLR"       = "#F8766D",
  "Naive long-run LLR" = "#4C9BF5",
  "Proposed method"    = "#38A365"
)

# ================================================================
# Main plotting function
# ================================================================

plot_dynamic_rd_table <- function(setting_to_plot = "Setting 1",
                                  point_size = 2.25,
                                  line_width = 0.8) {
  
  plot_df <- dynamic_rd_plot_data %>%
    filter(setting == setting_to_plot)
  
  ggplot(
    plot_df,
    aes(
      x = n,
      y = value,
      color = method,
      group = method
    )
  ) +
    geom_line(linewidth = line_width, lineend = "round") +
    geom_point(size = point_size) +
    facet_grid(
      metric ~ gamma_lab,
      scales = "free_y",
      labeller = labeller(
        metric = label_value,
        gamma_lab = label_parsed
      )
    ) +
    scale_x_log10(
      breaks = c(1e3, 1e4, 1e5),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_color_manual(values = .method_colors, drop = FALSE) +
    labs(
      x = expression("Sample size " * italic(n) * " (log scale)")
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.0,
          size = 2.5
        )
      )
    ) +
    .dynamic_rd_theme(base_size = 14)
}

# ================================================================
# Usage
# ================================================================

p_setting1 <- plot_dynamic_rd_table("Setting 1")
p_setting2 <- plot_dynamic_rd_table("Setting 2")

p_setting1
p_setting2

ggsave(
  "ar_setting1.pdf",
  p_setting1,
  width = 10.5,
  height = 5,
  device = grDevices::pdf,
  useDingbats = FALSE
)

ggsave(
  "ar_setting2.pdf",
  p_setting2,
  width = 10.5,
  height = 5,
  device = grDevices::pdf,
  useDingbats = FALSE
)
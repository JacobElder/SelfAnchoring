# make_param_recovery_fig.R
# Parameter Recovery Heatmap: true vs. recovered for Sym_Lambda and Asym_Lambda.
#
# Reads from Results/parameter_recovery_symlambda_results.csv and
# Results/parameter_recovery_asymlambda_results.csv when available.
# Falls back to simulated placeholder data if not yet run.
#
# Run `Parameter Recovery/run_parameter_recovery.R` first to generate real data.
# Output: Figures/fig_param_recovery.tiff

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})

TIFF_DPI   <- 300
TIFF_UNITS <- "in"

# ── Load or simulate data ─────────────────────────────────────────────────────
sym_file  <- here("Results", "parameter_recovery_symlambda_results.csv")
asym_file <- here("Results", "parameter_recovery_asymlambda_results.csv")

using_placeholder <- FALSE

load_csv <- function(path) {
  if (file.exists(path)) { read.csv(path) } else { NULL }
}

sym_data  <- load_csv(sym_file)
asym_data <- load_csv(asym_file)

# Placeholder: simulates realistic recovery (r ~ .80–.90) if CSVs absent
gen_placeholder <- function(model_name, params, n = 61, seed = 42) {
  set.seed(seed)
  param_ranges <- list(
    m      = c(0.3, 9.7, 2.80),   # min, max, mean
    m_in   = c(0.3, 9.7, 3.50),
    m_out  = c(0.3, 9.7, 2.50),
    bias   = c(0.05, 0.95, 0.36),
    lambda = c(0.3, 4.9,  3.56),
    w      = c(0.10, 0.90, 0.56)
  )
  do.call(rbind, lapply(params, function(p) {
    rng  <- param_ranges[[p]]
    rng_width <- rng[2] - rng[1]
    # Draw true params from a beta-ish distribution
    true_pr <- rnorm(n, qnorm((rng[3] - rng[1]) / rng_width), 0.5)
    true_vals <- pmin(rng[2], pmax(rng[1], pnorm(true_pr) * rng_width + rng[1]))
    # Recovered adds noise (recovery r ~ .83)
    noise     <- rnorm(n, 0, rng_width * 0.22)
    rec_vals  <- pmin(rng[2], pmax(rng[1], true_vals + noise))
    data.frame(
      model     = model_name,
      param     = p,
      subj_idx  = seq_len(n),
      true      = true_vals,
      recovered = rec_vals
    )
  }))
}

if (is.null(sym_data) || is.null(asym_data)) {
  using_placeholder <- TRUE
  message("  Recovery CSVs not found — using placeholder data.")
  message("  Run: caffeinate -i Rscript 'Parameter Recovery/run_parameter_recovery.R'")
  if (is.null(sym_data))
    sym_data  <- gen_placeholder("Sym_Lambda",  c("m", "bias", "lambda", "w"),          seed = 42)
  if (is.null(asym_data))
    asym_data <- gen_placeholder("Asym_Lambda", c("m_in", "m_out", "bias", "lambda", "w"), seed = 43)
} else {
  message("  Loaded real recovery data from Results/")
}

# ── Pretty parameter labels ───────────────────────────────────────────────────
param_label_map <- c(
  m      = "\u03b1 (projection rate)",
  m_in   = "\u03b1[in]",
  m_out  = "\u03b1[out]",
  bias   = "\u03b3 (ingroup bias)",
  lambda = "\u03bb (generalization)",
  w      = "w (lapse rate)"
)

# Desired panel order per model
sym_order  <- c("m", "bias", "lambda", "w")
asym_order <- c("m_in", "m_out", "bias", "lambda", "w")

label_data <- function(df, param_order) {
  df |>
    filter(param %in% param_order) |>
    mutate(
      param_label = factor(
        param_label_map[param],
        levels = unname(param_label_map[param_order])
      )
    )
}

sym_df  <- label_data(sym_data,  sym_order)
asym_df <- label_data(asym_data, asym_order)

# ── Recovery correlation per panel ────────────────────────────────────────────
r_labels <- function(df) {
  df |>
    group_by(param_label) |>
    summarise(
      r    = cor(true, recovered, use = "complete.obs"),
      xpos = quantile(true,      0.05, na.rm = TRUE),
      ypos = quantile(recovered, 0.97, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(lab = sprintf("r = %.2f", r))
}

r_sym  <- r_labels(sym_df)
r_asym <- r_labels(asym_df)

# ── Build one row (one model) ─────────────────────────────────────────────────
make_row <- function(df, r_df, row_title) {
  ggplot(df, aes(x = true, y = recovered)) +
    stat_bin2d(aes(fill = after_stat(count)), bins = 14, drop = FALSE) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "firebrick", linewidth = 0.65) +
    geom_text(data = r_df,
              aes(x = xpos, y = ypos, label = lab),
              inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 2.8, colour = "grey15", fontface = "bold") +
    scale_fill_distiller(palette = "Blues", direction = 1,
                         name = "Count", limits = c(1, NA), na.value = "grey96") +
    facet_wrap(~ param_label, nrow = 1, scales = "free") +
    labs(
      title = row_title,
      x     = "True parameter value",
      y     = "Recovered (posterior median)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.border      = element_rect(colour = "grey70", fill = NA, linewidth = 0.45),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey93"),
      strip.text        = element_text(face = "bold", size = 9),
      plot.title        = element_text(face = "bold", size = 10.5),
      plot.background   = element_rect(fill = "white", colour = NA),
      legend.position   = "right",
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 7.5)
    )
}

p_sym  <- make_row(sym_df,  r_sym,  "Symmetric + \u03bb  [m, \u03b3, \u03bb, w]")
p_asym <- make_row(asym_df, r_asym, "Asymmetric + \u03bb  [m_in, m_out, \u03b3, \u03bb, w]")

placeholder_note <- if (using_placeholder) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, hjust = 0.5, size = 3, colour = "grey50",
             label = "[Placeholder data — run Parameter Recovery/run_parameter_recovery.R to populate]") +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0, 0, 0, 0))
} else {
  NULL
}

if (!is.null(placeholder_note)) {
  fig_pr <- p_sym / p_asym / placeholder_note +
    plot_layout(heights = c(1, 1, 0.12))
} else {
  fig_pr <- p_sym / p_asym
}

fig_pr <- fig_pr +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

ggsave(here("Figures", "fig_param_recovery.tiff"),
       fig_pr, width = 11, height = 7,
       dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")

message(sprintf("Saved: Figures/fig_param_recovery.tiff  [%s]",
                if (using_placeholder) "PLACEHOLDER" else "real data"))

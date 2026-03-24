# make_param_recovery_fig.R
# Generates two parameter-recovery figures:
#
#   fig_param_recovery.tiff         — 2D-bin scatter (true vs. recovered), one row per model
#   fig_param_recovery_heatmap.tiff — confusion matrix heatmap, side-by-side models
#
# Reads:
#   Results/parameter_recovery_symlambda_results.csv
#   Results/parameter_recovery_asymlambda_results.csv
#
# Falls back to simulated placeholder data if CSVs are absent.
# Run `Parameter Recovery/run_parameter_recovery.R` first to populate real data.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})

TIFF_DPI   <- 300
TIFF_UNITS <- "in"

# ── Load or simulate data ─────────────────────────────────────────────────────
sym_file  <- here("Results", "parameter_recovery_symlambdanow_results.csv")
asym_file <- here("Results", "parameter_recovery_asymlambdanow_results.csv")

using_placeholder <- FALSE

load_csv <- function(path) {
  if (file.exists(path)) read.csv(path) else NULL
}

sym_data  <- load_csv(sym_file)
asym_data <- load_csv(asym_file)

gen_placeholder <- function(model_name, params, n = 61, seed = 42) {
  set.seed(seed)
  param_ranges <- list(
    m      = c(0.3, 9.7, 2.80),
    m_in   = c(0.3, 9.7, 3.50),
    m_out  = c(0.3, 9.7, 2.50),
    bias   = c(0.05, 0.95, 0.36),
    lambda = c(0.3, 4.9,  3.56),
    w      = c(0.10, 0.90, 0.56)
  )
  do.call(rbind, lapply(params, function(p) {
    rng       <- param_ranges[[p]]
    rng_width <- rng[2] - rng[1]
    true_pr   <- rnorm(n, qnorm((rng[3] - rng[1]) / rng_width), 0.5)
    true_vals <- pmin(rng[2], pmax(rng[1], pnorm(true_pr) * rng_width + rng[1]))
    noise     <- rnorm(n, 0, rng_width * 0.22)
    rec_vals  <- pmin(rng[2], pmax(rng[1], true_vals + noise))
    data.frame(model = model_name, param = p, subj_idx = seq_len(n),
               true = true_vals, recovered = rec_vals)
  }))
}

if (is.null(sym_data) || is.null(asym_data)) {
  using_placeholder <- TRUE
  message("  Recovery CSVs not found — using placeholder data.")
  message("  Run: caffeinate -i Rscript 'Parameter Recovery/run_parameter_recovery.R'")
  if (is.null(sym_data))
    sym_data  <- gen_placeholder("Sym_Lambda_NoW",  c("m", "bias", "lambda"),             seed = 42)
  if (is.null(asym_data))
    asym_data <- gen_placeholder("Asym_Lambda_NoW", c("m_in", "m_out", "bias", "lambda"), seed = 43)
} else {
  message("  Loaded real recovery data from Results/")
}

# ── Labels ────────────────────────────────────────────────────────────────────
param_label_map <- c(
  m      = "\u03b1 (projection rate)",
  m_in   = "\u03b1\u1d35\u207f",         # α_in
  m_out  = "\u03b1\u1d52\u1d58\u1d57",   # α_out
  bias   = "\u03b3 (ingroup bias)",
  lambda = "\u03bb (generalization)",
  w      = "w (lapse rate)"
)

# Short labels for heatmap axes
param_short_map <- c(
  m      = "\u03b1",
  m_in   = "\u03b1[in]",
  m_out  = "\u03b1[out]",
  bias   = "\u03b3",
  lambda = "\u03bb",
  w      = "w"
)

sym_order  <- c("m", "bias", "lambda")
asym_order <- c("m_in", "m_out", "bias", "lambda")

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

# ── Figure 1: 2D-bin scatter (true vs. recovered) ─────────────────────────────
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

make_scatter_row <- function(df, r_df, row_title) {
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
    labs(title = row_title,
         x     = "True parameter value",
         y     = "Recovered (posterior median)") +
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

p_sym  <- make_scatter_row(sym_df,  r_sym,  "Symmetric + \u03bb  [\u03b1, \u03b3, \u03bb]")
p_asym <- make_scatter_row(asym_df, r_asym, "Asymmetric + \u03bb  [\u03b1\u1d35\u207f, \u03b1\u1d52\u1d58\u1d57, \u03b3, \u03bb]")

placeholder_note <- if (using_placeholder) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, hjust = 0.5, size = 3, colour = "grey50",
             label = "[Placeholder — run Parameter Recovery/run_parameter_recovery.R to populate]") +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0, 0, 0, 0))
} else { NULL }

fig_scatter <- if (!is.null(placeholder_note)) {
  p_sym / p_asym / placeholder_note + plot_layout(heights = c(1, 1, 0.12))
} else {
  p_sym / p_asym
}

fig_scatter <- fig_scatter +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(here("Figures", "fig_param_recovery.tiff"),
       fig_scatter, width = 11, height = 7,
       dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
message(sprintf("Saved: Figures/fig_param_recovery.tiff  [%s]",
                if (using_placeholder) "PLACEHOLDER" else "real data"))

# ── Figure 2: Confusion-matrix heatmap ───────────────────────────────────────
# Each cell (i, j) = Pearson r between true_param_i and recovered_param_j.
# Diagonal = self-recovery (bold text, outlined tile).
# Colour scale: white at 0 (no correlation) → deep blue at +1,
#               red at −1 (parameter confusability in opposite direction).

make_heatmap_panel <- function(raw_df, param_order, panel_title) {

  short <- param_short_map[param_order]  # short unicode labels in order

  # ── Compute confusion matrix ─────────────────────────────────────────────
  true_wide <- raw_df |>
    filter(param %in% param_order) |>
    select(subj_idx, param, true) |>
    pivot_wider(names_from = param, values_from = true, names_prefix = "T_")

  rec_wide <- raw_df |>
    filter(param %in% param_order) |>
    select(subj_idx, param, recovered) |>
    pivot_wider(names_from = param, values_from = recovered, names_prefix = "R_")

  wide <- left_join(true_wide, rec_wide, by = "subj_idx") |> select(-subj_idx)

  # r(true_i, rec_j) for all i,j
  cor_df <- expand.grid(true_p = param_order, rec_p = param_order,
                        stringsAsFactors = FALSE) |>
    rowwise() |>
    mutate(r = cor(wide[[paste0("T_", true_p)]],
                   wide[[paste0("R_", rec_p)]],
                   use = "complete.obs")) |>
    ungroup() |>
    mutate(
      true_label = factor(short[true_p], levels = rev(short)),  # y-axis: top=first
      rec_label  = factor(short[rec_p],  levels = short),       # x-axis: left=first
      on_diag    = true_p == rec_p,
      r_text     = sprintf("%.2f", r)
    )

  ggplot(cor_df, aes(x = rec_label, y = true_label, fill = r)) +
    # Background tiles
    geom_tile(colour = "white", linewidth = 0.8) +
    # Diagonal outline
    geom_tile(data = filter(cor_df, on_diag),
              colour = "#1a1a1a", linewidth = 1.1, fill = NA) +
    # Correlation text — bold on diagonal
    geom_text(aes(label = r_text,
                  fontface = ifelse(on_diag, "bold", "plain")),
              size = 3.6, colour = ifelse(cor_df$r > 0.6 | cor_df$r < -0.6,
                                          "white", "grey20")) +
    scale_fill_gradient2(
      low      = "#c0392b",
      mid      = "white",
      high     = "#1a5276",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "r",
      breaks   = c(-1, -0.5, 0, 0.5, 1)
    ) +
    labs(
      title = panel_title,
      x     = "Recovered parameter",
      y     = "True parameter"
    ) +
    scale_x_discrete(position = "bottom") +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid     = element_blank(),
      axis.text      = element_text(size = 10, face = "bold"),
      axis.title     = element_text(size = 9),
      plot.title     = element_text(face = "bold", size = 10.5, hjust = 0.5),
      legend.title   = element_text(size = 8),
      legend.text    = element_text(size = 8),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

p_heat_sym  <- make_heatmap_panel(
  sym_data,  sym_order,
  "Symmetric + \u03bb"
)
p_heat_asym <- make_heatmap_panel(
  asym_data, asym_order,
  "Asymmetric + \u03bb"
)

fig_heatmap <- (p_heat_sym | p_heat_asym) +
  plot_annotation(
    title    = "Parameter Recovery: Confusion Matrix",
    subtitle = paste("Each cell = r(true\u1d62, recovered\u2C7C). Diagonal (outlined) = self-recovery.",
                     if (using_placeholder) " [PLACEHOLDER DATA]" else ""),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, colour = "grey40", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave(here("Figures", "fig_param_recovery_heatmap.tiff"),
       fig_heatmap, width = 9, height = 5,
       dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
message(sprintf("Saved: Figures/fig_param_recovery_heatmap.tiff  [%s]",
                if (using_placeholder) "PLACEHOLDER" else "real data"))

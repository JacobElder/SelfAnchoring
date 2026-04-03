# make_figures.R
# Generates all TIFF figures for ElderJacob_JPSP_Submission.qmd
#
# Prerequisites (run before this script):
#   - marginal_effects_s1.R, marginal_effects_s2.R, marginal_effects_s3.R
#     (creates Results/marginal_effects_s*_predictions.csv)
#   - run_model_comparison_s*.R (creates Fits/loo_s*_*.rds)
#   - LOO files in Fits/ (git-ignored; regenerate with run_model_comparison_*.R)
#
# Output: Figures/*.tiff (300 DPI, git-ignored; tracked via this R script)
#
# Figure inventory:
#   fig04_parameter_illustration.tiff — from/self-ratings to ingroup predictions
#   fig05_lambda_identifiability.tiff  — Shepard's law / lambda identifiability
#   fig07_marginal_effects.tiff        — 3×3 behavioral marginal effects grid
#   fig_elpd_raincloud_s1.tiff         — per-subject ELPD by model, Study 1
#   fig_elpd_raincloud_s2.tiff         — per-subject ELPD by model, Study 2
#   fig_elpd_raincloud_s3.tiff         — per-subject ELPD by model, Study 3
#   fig_elpd_raincloud_pooled.tiff     — per-subject ELPD by model, Pooled

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})

# ── Optional packages (graceful fallback) ─────────────────────────────────────
has_ggdist <- requireNamespace("ggdist", quietly = TRUE)
has_loo    <- requireNamespace("loo",    quietly = TRUE)
if (!has_ggdist) message("ggdist not installed — ELPD rainclouds will be skipped. Install with: install.packages('ggdist')")
if (!has_loo)    message("loo not installed — ELPD rainclouds will be skipped. Install with: install.packages('loo')")

dir.create(here("Figures"), showWarnings = FALSE)

# ── Shared visual theme ────────────────────────────────────────────────────────
theme_apa <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92"),
      legend.position   = "bottom",
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      strip.text        = element_text(face = "bold", size = base_size),
      axis.title        = element_text(size = base_size - 0.5),
      axis.text         = element_text(size = base_size - 1),
      plot.title        = element_text(face = "bold", size = base_size + 0.5)
    )
}

TIFF_DPI    <- 300
TIFF_UNITS  <- "in"

# ── Color palettes ─────────────────────────────────────────────────────────────
# Model colors (ELPD raincloud)
model_cols <- c(
  "Bias"    = "#6C757D",
  "Sym"     = "#2B5C8A",
  "Sym+λ"   = "#4E9A9A",
  "Asym+λ"  = "#C7522A"
)

# ── Colorblind-safe + greyscale-compatible palettes ───────────────────────────
# Okabe-Ito palette — safe for deuteranopia/protanopia, distinct in greyscale
# via luminance differences. Dual-coded with linetype for redundancy.

# Study 1 (no condition): single black line
s1_col      <- "black"
s1_lty      <- "solid"

# Study 2: 3 outgroup conditions (Negation=ref, High-Status=UCLA, Low-Status=CSU LA)
s2_cols     <- c("Not UCR" = "#E69F00",   # orange
                 "UCLA"    = "#0072B2",   # dark blue
                 "CSU LA"  = "#009E73")   # green
s2_ltys     <- c("Not UCR" = "solid", "UCLA" = "dashed", "CSU LA" = "dotdash")
s2_labels   <- c("Not UCR" = "Negation", "UCLA" = "High-Status", "CSU LA" = "Low-Status")

# Study 3: 2 conditions
s3_cols     <- c("Minority" = "#CC79A7",  # pink/mauve
                 "Majority" = "#0072B2")  # dark blue
s3_ltys     <- c("Minority" = "solid", "Majority" = "dashed")
s3_labels   <- c("Minority" = "Racial Minority Outgroup",
                 "Majority" = "Racial Majority Outgroup")

# ── Dissertation-style theme (white bg, black border, no grid, bold axes) ────
theme_dissert <- function(base_size = 11) {
  theme(
    # Panel
    panel.background  = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    # Axes
    axis.line         = element_blank(),   # border covers it
    axis.ticks        = element_line(colour = "black", linewidth = 0.4),
    axis.text         = element_text(size = base_size - 1, colour = "black"),
    axis.title        = element_text(size = base_size,     face = "bold", colour = "black"),
    # Legend
    legend.background = element_rect(fill = NA),
    legend.key        = element_rect(fill = NA),
    legend.text       = element_text(size = base_size - 1),
    legend.title      = element_blank(),
    # Strip / facet
    strip.background  = element_blank(),
    strip.text        = element_text(size = base_size, face = "bold"),
    # Titles
    plot.title        = element_text(size = base_size, face = "bold"),
    plot.background   = element_rect(fill = "white", colour = NA)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Marginal Effects: 3 Studies × 3 Predictors
# Smooth continuous prediction curves (dense grid) emulating ggpredict() style
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Figure 7: Marginal Effects ──")

read_preds <- function(path) {
  if (!file.exists(path)) {
    message("  Missing: ", path, " — re-run the marginal_effects script")
    return(NULL)
  }
  df <- read.csv(path)
  # Handle old format (only SS, 3 points, predicted.Z column)
  if (!"predictor" %in% names(df)) {
    if ("predicted.Z" %in% names(df) && !"outgroup" %in% names(df)) {
      df <- df |>
        dplyr::filter(is.na(novel)) |>
        dplyr::rename(predictor_val = predicted.Z) |>
        dplyr::mutate(predictor = "ss", model = "M3")
    } else if ("predicted.Z" %in% names(df) && "outgroup" %in% names(df)) {
      df <- df |>
        dplyr::rename(predictor_val = predicted.Z) |>
        dplyr::mutate(predictor = "ss", model = "M3")
    } else if ("predicted.Z" %in% names(df) && "condition" %in% names(df)) {
      df <- df |>
        dplyr::rename(predictor_val = predicted.Z) |>
        dplyr::mutate(predictor = "ss", model = "M3")
    }
  }
  df
}

# Make a single smooth prediction panel — dissertation style
make_smooth_panel <- function(df, predictor_filter,
                              condition_col  = NULL,
                              color_map      = NULL,
                              lty_map        = NULL,
                              label_map      = NULL,
                              y_lab          = "Probability of Ingroup Classification",
                              x_lab          = NULL,
                              legend_pos     = c(0.05, 0.75)) {

  if (is.null(df)) {
    return(ggplot() + theme_void() +
             annotate("text", x=0.5, y=0.5, label="Pending\n(re-run marginal_effects script)",
                      hjust=0.5, size=3, colour="grey60"))
  }

  pred_xlabs <- c(
    "desirability" = "Desirability (Z)",
    "selfResp"     = "Self-Evaluation (Z)",
    "ss"           = "Similarity-to-Self (Z)"
  )
  if (is.null(x_lab)) x_lab <- pred_xlabs[predictor_filter]

  d <- df |> dplyr::filter(predictor == predictor_filter)
  if (nrow(d) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x=0.5, y=0.5,
                      label = paste0(pred_xlabs[predictor_filter], "\nPending"),
                      hjust=0.5, size=3, colour="grey60"))
  }

  # Y limits from data (don't clip CI) + add small padding
  y_lo <- floor(min(d$conf.low,  na.rm=TRUE) * 20) / 20
  y_hi <- ceil2(max(d$conf.high, na.rm=TRUE) * 20) / 20
  y_lo <- max(0, y_lo)
  y_hi <- min(1, y_hi)

  if (!is.null(condition_col) && condition_col %in% names(d)) {
    d[[condition_col]] <- factor(d[[condition_col]])
    # Keep original factor levels so color/lty maps match; pass label_map to scale labels

    p <- ggplot(d, aes(x = predictor_val, y = estimate,
                       colour   = .data[[condition_col]],
                       fill     = .data[[condition_col]],
                       linetype = .data[[condition_col]])) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  colour = NA, alpha = 0.15) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = color_map, labels = label_map, name = NULL) +
      scale_fill_manual(values   = color_map, labels = label_map, name = NULL) +
      scale_linetype_manual(values = lty_map,  labels = label_map, name = NULL) +
      theme_dissert() +
      theme(legend.position = legend_pos,
            legend.justification = c("right","bottom"),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.9, "lines"),
            legend.spacing.y = unit(0.1, "cm"))
  } else {
    p <- ggplot(d, aes(x = predictor_val, y = estimate)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = "grey80", colour = NA, alpha = 0.6) +
      geom_line(colour = s1_col, linetype = s1_lty, linewidth = 0.9) +
      theme_dissert() +
      theme(legend.position = "none")
  }

  p +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    scale_y_continuous(labels = scales::percent_format(1),
                       limits = c(y_lo, y_hi),
                       expand = expansion(mult = 0.01)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    labs(x = x_lab, y = y_lab) +
    theme(plot.margin = margin(1, 2, 1, 0, "mm"))
}

ceil2 <- function(x) ceiling(x * 20) / 20   # round up to nearest 0.05

# Load and standardise prediction data
preds_s1 <- read_preds(here("Results", "marginal_effects_s1_predictions.csv"))
preds_s2 <- read_preds(here("Results", "marginal_effects_s2_predictions.csv"))
preds_s3 <- read_preds(here("Results", "marginal_effects_s3_predictions.csv"))

# ── Build all 9 panels ────────────────────────────────────────────────────────
y_lab_short <- "P(Ingroup)"

# Row 1 — Study 1 (no condition grouping)
p1_des  <- make_smooth_panel(preds_s1, "desirability", y_lab = y_lab_short)
p1_self <- make_smooth_panel(preds_s1, "selfResp",     y_lab = y_lab_short)
p1_ss   <- make_smooth_panel(preds_s1, "ss",           y_lab = y_lab_short)

# Row 2 — Study 2 (by outgroup condition)
p2_des  <- make_smooth_panel(preds_s2, "desirability", condition_col="outgroup",
                              color_map=s2_cols, lty_map=s2_ltys, label_map=s2_labels,
                              y_lab = y_lab_short, legend_pos = c(0.97, 0.08))
p2_self <- make_smooth_panel(preds_s2, "selfResp",     condition_col="outgroup",
                              color_map=s2_cols, lty_map=s2_ltys, label_map=s2_labels,
                              y_lab = y_lab_short, legend_pos = c(0.97, 0.08))
p2_ss   <- make_smooth_panel(preds_s2, "ss",           condition_col="outgroup",
                              color_map=s2_cols, lty_map=s2_ltys, label_map=s2_labels,
                              y_lab = y_lab_short, legend_pos = c(0.97, 0.08))

# Row 3 — Study 3 (by majority/minority condition)
p3_des  <- make_smooth_panel(preds_s3, "desirability", condition_col="condition",
                              color_map=s3_cols, lty_map=s3_ltys, label_map=s3_labels,
                              y_lab = y_lab_short, legend_pos = c(0.97, 0.08))
p3_self <- make_smooth_panel(preds_s3, "selfResp",     condition_col="condition",
                              color_map=s3_cols, lty_map=s3_ltys, label_map=s3_labels,
                              y_lab = y_lab_short, legend_pos = c(0.97, 0.08))
p3_ss   <- make_smooth_panel(preds_s3, "ss",           condition_col="condition",
                              color_map=s3_cols, lty_map=s3_ltys, label_map=s3_labels,
                              y_lab = y_lab_short, legend_pos = c(0.97, 0.08))

# ── Column-title rows (invisible spacer plots with bold label) ────────────────
col_label <- function(txt) {
  ggplot() + theme_void() +
    annotate("text", x=0.5, y=0.5, label=txt, fontface="bold", size=3.8) +
    theme(plot.margin = margin(0,0,0,0))
}

# Row labels via left-side spacer plots (no in-panel titles per APA 7)
row_label <- function(txt) {
  ggplot() + theme_void() +
    annotate("text", x=0.5, y=0.5, label=txt, fontface="bold", size=3.4,
             angle=90, hjust=0.5, vjust=0.5) +
    theme(plot.margin = margin(0,0,0,0))
}

# ── Embed study labels as in-panel titles (no spacer column) ──────────────────
embed_title <- function(p, title_txt) {
  p + labs(title = title_txt) +
    theme(plot.title = element_text(face = "bold", size = 8.5, hjust = 0,
                                    margin = margin(b = 1)))
}

p1_des_t <- embed_title(p1_des,  "Study 1 (Minimal Groups)")
p2_des_t <- embed_title(p2_des,  "Study 2 (University Groups)")
p3_des_t <- embed_title(p3_des,  "Study 3 (Racial Groups)")

# 3-column layout: col1 = Desirability, col2 = Self-Evaluation, col3 = Similarity-to-Self
# Column headers are carried by x-axis labels set in make_smooth_panel
fig7 <- (p1_des_t | p1_self | p1_ss) /
        (p2_des_t | p2_self | p2_ss) /
        (p3_des_t | p3_self | p3_ss) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

ggsave(here("Figures", "fig07_marginal_effects.tiff"),
       fig7, width = 11, height = 10, dpi = TIFF_DPI, units = TIFF_UNITS,
       compression = "lzw")
message("  Saved: Figures/fig07_marginal_effects.tiff")


# ══════════════════════════════════════════════════════════════════════════════
# ELPD SLOPEGRAPH — Per-Subject LOO-ELPD by Model, per Study
# Thin connecting lines show within-subject trajectories across model complexity;
# boxplots overlay the aggregate distribution per model.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── ELPD Slopegraphs ──")

if (!has_loo) {
  message("  Skipping ELPD slopegraphs (loo not available). Install with: install.packages('loo')")
} else {
  library(loo)

  # Helper: sum pointwise ELPD within subjects using pareto_k subject mapping
  get_subj_elpd <- function(loo_file, pk_df) {
    if (!file.exists(loo_file)) return(NULL)
    loo_obj <- readRDS(loo_file)
    pw      <- loo_obj$pointwise[, "elpd_loo"]

    pk <- pk_df |> arrange(subj_idx)
    n_trials    <- pk$n_trials
    trial_end   <- cumsum(n_trials)
    trial_start <- c(1L, head(trial_end, -1L) + 1L)

    subj_elpd <- sapply(seq_along(n_trials), function(i) {
      idx <- trial_start[i]:trial_end[i]
      if (any(idx > length(pw))) return(NA_real_)
      sum(pw[idx])
    })
    tibble(subID = pk$subID, subj_elpd = subj_elpd)
  }

  make_elpd_slopegraph <- function(elpd_long, study_label, available_models,
                                    condition_col = NULL, condition_labels = NULL,
                                    condition_colors = NULL) {
    n_models_possible <- 4
    if (length(available_models) < n_models_possible) {
      pending <- setdiff(c("Bias","Sym","Sym+λ","Asym+λ"), available_models)
      caption_txt <- paste0("Note. ", paste(pending, collapse = ", "),
                            " model(s) pending — not yet run for ", study_label, ".")
    } else {
      caption_txt <- paste0(
        "Note. Each line = one participant's summed LOO-ELPD across model architectures. ",
        "Higher (less negative) = better predictive fit. Boxplots show median and IQR; ",
        "diamonds = group mean."
      )
    }

    elpd_long$model <- factor(elpd_long$model, levels = names(model_cols))

    # Group means for diamond overlay (grouped by condition if faceting)
    group_vars <- if (!is.null(condition_col)) c("model", condition_col) else "model"
    mean_df <- elpd_long |>
      group_by(across(all_of(group_vars))) |>
      summarise(mean_elpd = mean(subj_elpd, na.rm = TRUE), .groups = "drop")

    # Condition color coding for connecting lines (same palette as fig07)
    use_cond_color <- !is.null(condition_col) && !is.null(condition_colors)

    if (use_cond_color) {
      line_geom <- geom_line(aes(group = subID, color = .data[[condition_col]]),
                             alpha = 0.45, linewidth = 0.3)
    } else {
      line_geom <- geom_line(aes(group = subID), color = "grey45",
                             alpha = 0.40, linewidth = 0.3)
    }

    p <- ggplot(elpd_long, aes(x = model, y = subj_elpd)) +
      line_geom +
      # Colored boxplots per model (no outlier points — trajectories show them)
      geom_boxplot(
        aes(fill = model),
        width = 0.38, outlier.shape = NA, alpha = 0.70,
        color = "grey25", linewidth = 0.45
      ) +
      # Group mean diamonds (filled by model, grey border)
      geom_point(data = mean_df,
                 aes(x = model, y = mean_elpd, fill = model),
                 size = 3.5, shape = 23, color = "grey25") +
      scale_fill_manual(values = model_cols, guide = "none") +
      scale_x_discrete(labels = c(
        "Bias"   = "Bias",
        "Sym"    = "Symmetric",
        "Sym+λ"  = "Sym+λ",
        "Asym+λ" = "Asym+λ"
      )) +
      labs(
        x = NULL,
        y = "Per-Subject LOO-ELPD"
      ) +
      theme_dissert() +
      theme(legend.position = "none")

    # Add condition color scale if using condition-coded lines
    if (use_cond_color) {
      p <- p + scale_color_manual(values = condition_colors, guide = "none")
    }

    # Add condition facets if requested
    if (!is.null(condition_col)) {
      labeller_fn <- if (!is.null(condition_labels)) {
        ggplot2::as_labeller(condition_labels)
      } else {
        ggplot2::label_value
      }
      p <- p + facet_wrap(as.formula(paste("~", condition_col)),
                           labeller = labeller_fn, nrow = 1)
    }

    p
  }

  # ── Study 1 ─────────────────────────────────────────────────────────────────
  pk_s1 <- tryCatch(
    read.csv(here("Results","pareto_k_subj_s1_sym_lambda.csv")) |> arrange(subj_idx),
    error = function(e) { message("  Cannot read S1 pareto_k file"); NULL }
  )

  if (!is.null(pk_s1)) {
    models_s1 <- list(
      "Bias"    = here("Fits","loo_s1_bias.rds"),
      "Sym"     = here("Fits","loo_s1_symmetric.rds"),
      "Sym+λ"   = here("Fits","loo_s1_sym_lambda.rds"),
      "Asym+λ"  = here("Fits","loo_s1_asym_lambda.rds")
    )

    elpd_s1 <- imap_dfr(models_s1, function(path, label) {
      res <- get_subj_elpd(path, pk_s1)
      if (is.null(res)) return(NULL)
      res |> mutate(model = label)
    })

    if (nrow(elpd_s1) > 0) {
      available_s1 <- unique(elpd_s1$model)
      p_s1 <- make_elpd_slopegraph(elpd_s1, "Study 1 (Minimal Groups)", available_s1)
      ggsave(here("Figures","fig_elpd_raincloud_s1.tiff"),
             p_s1, width = 7, height = 5.5, dpi = TIFF_DPI, units = TIFF_UNITS,
             compression = "lzw")
      message("  Saved: Figures/fig_elpd_raincloud_s1.tiff")
    }
  }

  # ── Study 2 ─────────────────────────────────────────────────────────────────
  pk_s2 <- tryCatch(
    read.csv(here("Results","pareto_k_subj_s2_bias.csv")) |> arrange(subj_idx),
    error = function(e) { message("  Cannot read S2 pareto_k file"); NULL }
  )

  if (!is.null(pk_s2)) {
    models_s2 <- list(
      "Bias"    = here("Fits","loo_s2_bias.rds"),
      "Sym"     = here("Fits","loo_s2_symmetric.rds"),
      "Sym+λ"   = here("Fits","loo_s2_sym_lambda.rds"),
      "Asym+λ"  = here("Fits","loo_s2_asym_lambda.rds")
    )

    elpd_s2 <- imap_dfr(models_s2, function(path, label) {
      res <- get_subj_elpd(path, pk_s2)
      if (is.null(res)) return(NULL)
      res |> mutate(model = label)
    })

    if (nrow(elpd_s2) > 0) {
      # Join outgroup condition for faceting
      cond_s2 <- tryCatch(
        read.csv(here("Study 2","Cleaning","output","fullTest.csv")) |>
          dplyr::distinct(subID, outgroup),
        error = function(e) NULL
      )
      if (!is.null(cond_s2)) elpd_s2 <- dplyr::left_join(elpd_s2, cond_s2, by = "subID")

      available_s2 <- unique(elpd_s2$model)
      s2_cond_labels <- c("Not UCR"="Negation", "UCLA"="High-Status", "CSU LA"="Low-Status")
      p_s2 <- make_elpd_slopegraph(elpd_s2, "Study 2 (University Groups)", available_s2,
                                    condition_col     = if ("outgroup" %in% names(elpd_s2)) "outgroup" else NULL,
                                    condition_labels  = s2_cond_labels,
                                    condition_colors  = s2_cols)
      ggsave(here("Figures","fig_elpd_raincloud_s2.tiff"),
             p_s2, width = 10, height = 5.5, dpi = TIFF_DPI, units = TIFF_UNITS,
             compression = "lzw")
      message("  Saved: Figures/fig_elpd_raincloud_s2.tiff")
    }
  }

  # ── Study 3 ─────────────────────────────────────────────────────────────────
  pk_s3 <- tryCatch(
    read.csv(here("Results","pareto_k_subj_s3_bias.csv")) |> arrange(subj_idx),
    error = function(e) { message("  Cannot read S3 pareto_k file"); NULL }
  )

  if (!is.null(pk_s3)) {
    models_s3 <- list(
      "Bias"    = here("Fits","loo_s3_bias.rds"),
      "Sym"     = here("Fits","loo_s3_symmetric.rds"),
      "Sym+λ"   = here("Fits","loo_s3_sym_lambda.rds"),
      "Asym+λ"  = here("Fits","loo_s3_asym_lambda.rds")
    )

    elpd_s3 <- imap_dfr(models_s3, function(path, label) {
      res <- get_subj_elpd(path, pk_s3)
      if (is.null(res)) return(NULL)
      res |> mutate(model = label)
    })

    if (nrow(elpd_s3) > 0) {
      # Join condition for faceting
      cond_s3 <- tryCatch(
        read.csv(here("Study 3","Cleaning","output","fullTest_fixed.csv")) |>
          dplyr::distinct(subID, condition),
        error = function(e) NULL
      )
      if (!is.null(cond_s3)) elpd_s3 <- dplyr::left_join(elpd_s3, cond_s3, by = "subID")

      available_s3 <- unique(elpd_s3$model)
      s3_cond_labels <- c("Minority"="Racial Minority Outgroup", "Majority"="Racial Majority Outgroup")
      p_s3 <- make_elpd_slopegraph(elpd_s3, "Study 3 (Racial Groups)", available_s3,
                                    condition_col     = if ("condition" %in% names(elpd_s3)) "condition" else NULL,
                                    condition_labels  = s3_cond_labels,
                                    condition_colors  = s3_cols)
      ggsave(here("Figures","fig_elpd_raincloud_s3.tiff"),
             p_s3, width = 8.5, height = 5.5, dpi = TIFF_DPI, units = TIFF_UNITS,
             compression = "lzw")
      message("  Saved: Figures/fig_elpd_raincloud_s3.tiff")
    }
  }
  # ── Pooled ──────────────────────────────────────────────────────────────────
  pk_pooled <- tryCatch(
    read.csv(here("Results","pareto_k_subj_pooled_bias.csv")) |> arrange(subj_idx),
    error = function(e) { message("  Cannot read pooled pareto_k file"); NULL }
  )

  if (!is.null(pk_pooled)) {
    models_pooled <- list(
      "Bias"    = here("Fits","loo_pooled_bias.rds"),
      "Sym"     = here("Fits","loo_pooled_symmetric.rds"),
      "Sym+λ"   = here("Fits","loo_pooled_sym_lambda.rds"),
      "Asym+λ"  = here("Fits","loo_pooled_asym_lambda.rds")
    )

    elpd_pooled <- imap_dfr(models_pooled, function(path, label) {
      res <- get_subj_elpd(path, pk_pooled)
      if (is.null(res)) return(NULL)
      res |> mutate(model = label)
    })

    if (nrow(elpd_pooled) > 0) {
      available_pooled <- unique(elpd_pooled$model)
      p_pooled <- make_elpd_slopegraph(elpd_pooled, "Pooled (N = 609)", available_pooled)
      ggsave(here("Figures","fig_elpd_raincloud_pooled.tiff"),
             p_pooled, width = 7, height = 5.5, dpi = TIFF_DPI, units = TIFF_UNITS,
             compression = "lzw")
      message("  Saved: Figures/fig_elpd_raincloud_pooled.tiff")
    }
  }

}  # end loo block


# ── Load S1 summary (used by gradient panels below) ──────────────────────────
summary_s1 <- tryCatch(
  read.csv(here("Results","summary_s1_sym_lambda.csv")),
  error = function(e) { message("  Missing summary_s1_sym_lambda.csv"); NULL }
)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Multi-Study Generalization Gradients
# S^λ decay curves by study and by condition within study.
# Layout: 3 panels (S1, S2, S3) sharing x/y axes.
# Individual participant gradients in background; condition means as solid lines.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Multi-Study Generalization Gradients ──")

S_seq <- seq(0, 1, by = 0.005)

# ── Load condition means from desc CSVs ──────────────────────────────────────
read_cond_lambda <- function(path, condition_map = NULL) {
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path)
  df <- df[df$param == "lambda", c("condition","M","SD","n")]
  if (!is.null(condition_map)) df$condition <- condition_map[df$condition]
  df
}

# S2 condition labels
s2_map <- c("Not UCR" = "Negation", "UCLA" = "High-Status", "CSU LA" = "Low-Status")
# (desc CSV already uses these labels directly)

cond_s2 <- read_cond_lambda(here("Results","param_by_condition_s2_desc.csv")) |>
  dplyr::mutate(condition = dplyr::recode(condition,
    "Negation" = "Not UCR", "High-Status" = "UCLA", "Low-Status" = "CSU LA"
  ))
cond_s3 <- read_cond_lambda(here("Results","param_by_condition_s3_desc.csv"))

# ── Load individual λ by condition ───────────────────────────────────────────
ind_s2_grad <- tryCatch(
  read.csv(here("Results","ind_diffs_s2_full.csv")) |>
    dplyr::select(subj_idx, outgroup, lambda),
  error = function(e) NULL
)
ind_s3_grad <- tryCatch(
  read.csv(here("Results","ind_diffs_s3_full.csv")) |>
    dplyr::select(subj_idx, condition, lambda),
  error = function(e) NULL
)

# ── Helper: build gradient data frame from a tibble of condition × lambda ─────
make_grad_df <- function(cond_df, cond_col = "condition") {
  do.call(rbind, lapply(seq_len(nrow(cond_df)), function(i) {
    lam    <- cond_df$M[i]
    lam_lo <- max(0.01, lam - cond_df$SD[i])
    lam_hi <- lam + cond_df$SD[i]
    data.frame(
      S       = S_seq,
      g_med   = S_seq ^ lam,
      g_lo    = S_seq ^ lam_hi,  # higher λ = steeper decay = lower g
      g_hi    = S_seq ^ lam_lo,
      condition = cond_df[[cond_col]][i]
    )
  }))
}

# ── Panel S1 ─────────────────────────────────────────────────────────────────
make_grad_panel_s1 <- function() {
  if (is.null(summary_s1)) return(ggplot() + theme_void())

  lambda_row <- summary_s1 |> dplyr::filter(grepl("^\"?mu_pr\\[3\\]", variable))
  if (nrow(lambda_row) == 0) return(ggplot() + theme_void())

  lam_med <- pnorm(as.numeric(lambda_row$median)) * 5
  lam_lo  <- pnorm(as.numeric(lambda_row$q5))     * 5
  lam_hi  <- pnorm(as.numeric(lambda_row$q95))    * 5

  grad_s1 <- data.frame(
    S    = S_seq,
    g_med = S_seq ^ lam_med,
    g_lo  = S_seq ^ lam_hi,
    g_hi  = S_seq ^ lam_lo
  )

  # Individual gradients
  params_s1 <- tryCatch(read.csv(here("Results","params_ind_s1_sym_lambda.csv")), error = function(e) NULL)
  ind_lam_s1 <- if (!is.null(params_s1)) {
    params_s1 |>
      dplyr::filter(grepl("^\"?lambda\\[", variable)) |>
      dplyr::pull(median)
  } else { numeric(0) }

  p <- ggplot()
  if (length(ind_lam_s1) > 0) {
    ind_df <- do.call(rbind, lapply(seq_along(ind_lam_s1), function(i)
      data.frame(S = S_seq, g = S_seq ^ ind_lam_s1[i], subj_i = i)))
    p <- p + geom_line(data = ind_df, aes(x = S, y = g, group = subj_i),
                       colour = s1_col, alpha = 0.06, linewidth = 0.25)
  }
  p +
    geom_ribbon(data = grad_s1, aes(x = S, ymin = g_lo, ymax = g_hi),
                fill = s1_col, alpha = 0.20) +
    geom_line(data = grad_s1, aes(x = S, y = g_med),
              colour = s1_col, linewidth = 1.1) +
    annotate("text", x = 0.62, y = 0.87,
             label = sprintf("\u03bb = %.2f [%.2f, %.2f]", lam_med, lam_lo, lam_hi),
             hjust = 0, size = 3, colour = s1_col) +
    scale_x_continuous(name = "Semantic Similarity (Dice)", breaks = seq(0,1,0.25)) +
    scale_y_continuous(name = expression("Generalization Weight " ~ (S^lambda)), limits = c(0,1), breaks = seq(0,1,0.2)) +
    theme_dissert(base_size = 10) +
    theme(legend.position = "none")
}

# ── Panel S2 ─────────────────────────────────────────────────────────────────
make_grad_panel_s2 <- function() {
  if (is.null(cond_s2)) return(ggplot() + theme_void() +
                                annotate("text", x=0.5, y=0.5, label="S2 pending", hjust=0.5, size=3, colour="grey50"))

  grad_s2 <- make_grad_df(cond_s2, "condition")
  grad_s2$condition <- factor(grad_s2$condition, levels = names(s2_cols))

  ind_ribbons <- if (!is.null(ind_s2_grad)) {
    ind_s2_grad |>
      dplyr::filter(!is.na(outgroup)) |>
      dplyr::group_by(outgroup) |>
      dplyr::group_modify(~ {
        do.call(rbind, lapply(seq_len(nrow(.x)), function(i)
          data.frame(S = S_seq, g = S_seq ^ .x$lambda[i], subj_i = i)))
      }) |>
      dplyr::ungroup()
  } else { NULL }

  p <- ggplot()
  if (!is.null(ind_ribbons)) {
    p <- p + geom_line(data = ind_ribbons, aes(x = S, y = g, group = interaction(outgroup, subj_i),
                                                colour = outgroup),
                       alpha = 0.04, linewidth = 0.2)
  }
  p +
    geom_ribbon(data = grad_s2,
                aes(x = S, ymin = g_lo, ymax = g_hi, fill = condition),
                alpha = 0.18) +
    geom_line(data = grad_s2,
              aes(x = S, y = g_med, colour = condition, linetype = condition),
              linewidth = 1.0) +
    scale_colour_manual(values = s2_cols,  labels = s2_labels, name = NULL) +
    scale_fill_manual(values   = s2_cols,  labels = s2_labels, name = NULL) +
    scale_linetype_manual(values = s2_ltys, labels = s2_labels, name = NULL) +
    scale_x_continuous(name = "Semantic Similarity (Dice)", breaks = seq(0,1,0.25)) +
    scale_y_continuous(name = expression("Generalization Weight " ~ (S^lambda)), limits = c(0,1), breaks = seq(0,1,0.2)) +
    theme_dissert(base_size = 10) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right","top"),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.9, "lines"))
}

# ── Panel S3 ─────────────────────────────────────────────────────────────────
make_grad_panel_s3 <- function() {
  if (is.null(cond_s3)) return(ggplot() + theme_void() +
                                annotate("text", x=0.5, y=0.5, label="S3 pending", hjust=0.5, size=3, colour="grey50"))

  grad_s3 <- make_grad_df(cond_s3, "condition")
  grad_s3$condition <- factor(grad_s3$condition, levels = names(s3_cols))

  ind_ribbons3 <- if (!is.null(ind_s3_grad)) {
    ind_s3_grad |>
      dplyr::filter(!is.na(condition)) |>
      dplyr::group_by(condition) |>
      dplyr::group_modify(~ {
        do.call(rbind, lapply(seq_len(nrow(.x)), function(i)
          data.frame(S = S_seq, g = S_seq ^ .x$lambda[i], subj_i = i)))
      }) |>
      dplyr::ungroup()
  } else { NULL }

  p <- ggplot()
  if (!is.null(ind_ribbons3)) {
    p <- p + geom_line(data = ind_ribbons3,
                       aes(x = S, y = g, group = interaction(condition, subj_i),
                           colour = condition),
                       alpha = 0.04, linewidth = 0.2)
  }
  p +
    geom_ribbon(data = grad_s3,
                aes(x = S, ymin = g_lo, ymax = g_hi, fill = condition),
                alpha = 0.18) +
    geom_line(data = grad_s3,
              aes(x = S, y = g_med, colour = condition, linetype = condition),
              linewidth = 1.0) +
    scale_colour_manual(values = s3_cols,  labels = s3_labels, name = NULL) +
    scale_fill_manual(values   = s3_cols,  labels = s3_labels, name = NULL) +
    scale_linetype_manual(values = s3_ltys, labels = s3_labels, name = NULL) +
    scale_x_continuous(name = "Semantic Similarity (Dice)", breaks = seq(0,1,0.25)) +
    scale_y_continuous(name = expression("Generalization Weight " ~ (S^lambda)), limits = c(0,1), breaks = seq(0,1,0.2)) +
    theme_dissert(base_size = 10) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right","top"),
          legend.text = element_text(size = 7.5),
          legend.key.size = unit(0.9, "lines"))
}

pg1 <- make_grad_panel_s1()
pg2 <- make_grad_panel_s2()
pg3 <- make_grad_panel_s3()

# Embed study label as a compact title — no spacer column needed
grad_title_theme <- theme(
  plot.title = element_text(face = "bold", size = 8.5, hjust = 0,
                             margin = margin(b = 1, t = 1))
)
pg1_t <- pg1 + labs(title = "Study 1 (Minimal Groups)")   + grad_title_theme
pg2_t <- pg2 + labs(title = "Study 2 (University Groups)") + grad_title_theme
pg3_t <- pg3 + labs(title = "Study 3 (Racial Groups)")    + grad_title_theme

fig_grad_all <- (pg1_t / pg2_t / pg3_t) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

ggsave(here("Figures","fig_generalization_gradients.tiff"),
       fig_grad_all, width = 6, height = 8,
       dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
message("  Saved: Figures/fig_generalization_gradients.tiff")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — MCR × Individual Differences (Study 1 populated; S2/S3 placeholders)
# subject_mcr from asym_lambda model (tied with sym_lambda winner; MCR GQ not
# present in sym_lambda Stan model).  Key correlate: SING.Ind r = .53 p < .001.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── MCR × Individual Differences ──")

ind_diffs_s1 <- tryCatch(
  read.csv(here("Results","ind_diffs_s1_full.csv")),
  error = function(e) { message("  Missing ind_diffs_s1_full.csv"); NULL }
)

make_placeholder_panel <- function(label) {
  ggplot() + theme_void() +
    annotate("rect", xmin=0, xmax=1, ymin=0, ymax=1,
             fill="grey95", color="grey70", linewidth=0.5) +
    annotate("text", x=0.5, y=0.5, label=label,
             hjust=0.5, vjust=0.5, size=3.2, color="grey50") +
    theme(plot.background = element_rect(fill="white", colour=NA))
}

if (!is.null(ind_diffs_s1)) {
  # ── Study 1: MCR vs. SING.Ind (strongest correlate, r = .53) ────────────────
  r_val <- cor.test(ind_diffs_s1$subject_mcr, ind_diffs_s1$SING.Ind)
  r_lab <- sprintf("r = %.2f, p < .001", r_val$estimate)

  p_mcr_s1 <- ggplot(ind_diffs_s1, aes(x = SING.Ind, y = subject_mcr)) +
    geom_point(color = s1_col, size = 2, alpha = 0.65) +
    geom_smooth(method = "lm", se = TRUE, color = s1_col,
                fill = "grey75", linewidth = 0.9) +
    annotate("label", x = Inf, y = Inf,
             label = r_lab, parse = FALSE,
             hjust = 1.05, vjust = 1.3, size = 3.2, color = "grey20",
             fill = "white", linewidth = 0.3) +
    scale_x_continuous(name = "Social Identity Importance\n(SING Independence)") +
    scale_y_continuous(name = "Metacontrast Ratio") +
    theme_dissert()

  p_mcr_s2 <- make_placeholder_panel(
    "Study 2\n(University Groups)\nMCR pending\nS2 model completion"
  )
  p_mcr_s3 <- make_placeholder_panel(
    "Study 3\n(Racial Groups)\nMCR pending\nS3 model completion"
  )

  fig_mcr <- p_mcr_s1 | p_mcr_s2 | p_mcr_s3

  message("  Skipped: fig_mcr_indiff.tiff (removed from manuscript)")
} else {
  message("  Skipping MCR figure (ind_diffs_s1_full.csv not found)")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Task Schematic (convert JPG to TIFF)
# Source: Figures/Old Figures/SA_TaskSchematic/TaskSchematic.jpg
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Figure 2: Task Schematic ──")

convert_jpg_to_tiff <- function(jpg_path, tiff_path) {
  if (!file.exists(jpg_path)) {
    message("  Missing: ", jpg_path)
    return(invisible(NULL))
  }
  has_magick <- requireNamespace("magick", quietly = TRUE)
  if (has_magick) {
    img <- magick::image_read(jpg_path)
    magick::image_write(img, path = tiff_path, format = "tiff",
                        quality = NULL, density = TIFF_DPI)
    message("  Saved: ", tiff_path)
  } else {
    tryCatch({
      img_data <- jpeg::readJPEG(jpg_path)
      tiff(tiff_path, width = dim(img_data)[2], height = dim(img_data)[1],
           units = "px", res = TIFF_DPI, compression = "lzw")
      grid::grid.raster(img_data)
      dev.off()
      message("  Saved via grDevices: ", tiff_path)
    }, error = function(e) message("  Could not convert (install magick): ", e$message))
  }
}

convert_jpg_to_tiff(
  here("Figures", "Old Figures", "SA_TaskSchematic", "TaskSchematic.jpg"),
  here("Figures", "fig02_task_schematic.tiff")
)

# ── Notes on Figures 3, 4, 5 ─────────────────────────────────────────────────
# These are generated by their own scripts which save TIFF directly via ggsave:
#   fig03: Scripts/plot_parameter_space.R
#   fig04: Figures/make_parameter_illustration.R
#   fig05: Scripts/plot_lambda_explanation.R
# Run those scripts independently to regenerate.


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Cross-Study Parameter Forest Plot
# Group-level posterior medians + 90% CIs for α, γ, λ across S1, S2, S3.
# Back-transforms mu_pr from probit scale using pnorm() * scale.
# Shows parameter stability (or divergence) across intergroup contexts.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Cross-Study Parameter Forest Plot ──")

# Back-transform helper: pnorm(x) * scale
bt <- function(x, scale = 1) pnorm(x) * scale

extract_group_params <- function(csv_path, study_label) {
  if (!file.exists(csv_path)) {
    message("  Missing: ", csv_path); return(NULL)
  }
  sum_df <- read.csv(csv_path)
  # mu_pr order for sym_lambda NoW: [m, bias, lambda] → indices 1-3
  # variables named "mu_pr[1]" through "mu_pr[3]"
  params <- list(
    list(name = "\u03b1",       idx = 1, scale = 10),
    list(name = "\u03b3",       idx = 2, scale = 1),
    list(name = "\u03bb",       idx = 3, scale = 5)
  )
  do.call(rbind, lapply(params, function(p) {
    pat  <- paste0("mu_pr\\[", p$idx, "\\]")
    row  <- sum_df[grepl(pat, sum_df$variable), ]
    if (nrow(row) == 0) return(NULL)
    data.frame(
      study  = study_label,
      param  = p$name,
      median = bt(as.numeric(row$median), p$scale),
      lo     = bt(as.numeric(row$q5),     p$scale),
      hi     = bt(as.numeric(row$q95),    p$scale)
    )
  }))
}

gp_s1 <- extract_group_params(here("Results","summary_s1_sym_lambda.csv"), "Study 1\n(Minimal Groups)")
gp_s2 <- extract_group_params(here("Results","summary_s2_sym_lambda.csv"), "Study 2\n(University Groups)")
gp_s3 <- extract_group_params(here("Results","summary_s3_sym_lambda.csv"), "Study 3\n(Racial Groups)")

gp_all <- do.call(rbind, Filter(Negate(is.null), list(gp_s1, gp_s2, gp_s3)))

if (!is.null(gp_all) && nrow(gp_all) > 0) {

  study_cols <- c(
    "Study 1\n(Minimal Groups)"    = "#2B5C8A",
    "Study 2\n(University Groups)" = "#E69F00",
    "Study 3\n(Racial Groups)"     = "#CC79A7"
  )
  study_shapes <- c(
    "Study 1\n(Minimal Groups)"    = 16,
    "Study 2\n(University Groups)" = 17,
    "Study 3\n(Racial Groups)"     = 15
  )

  param_levels <- c("\u03b1", "\u03b3", "\u03bb")
  param_xlabs  <- c(
    "\u03b1" = "Projection Rate (\u03b1)\n[0 – 10]",
    "\u03b3" = "Ingroup Bias (\u03b3)\n[0 – 1]",
    "\u03bb" = "Generalization Sensitivity (\u03bb)\n[0 – 5]"
  )

  gp_all$param  <- factor(gp_all$param,  levels = param_levels)
  gp_all$study  <- factor(gp_all$study,  levels = names(study_cols))

  forest_panels <- lapply(param_levels, function(p) {
    d    <- gp_all[gp_all$param == p, ]
    xlab <- param_xlabs[p]
    ggplot(d, aes(x = median, y = study, colour = study, shape = study)) +
      geom_errorbarh(aes(xmin = lo, xmax = hi),
                     height = 0.18, linewidth = 0.8) +
      geom_point(size = 3.5) +
      scale_colour_manual(values = study_cols, guide = "none") +
      scale_shape_manual(values = study_shapes, guide = "none") +
      labs(x = xlab, y = NULL) +
      theme_dissert(base_size = 10) +
      theme(axis.text.y = element_text(size = 8.5),
            plot.margin = margin(4, 8, 4, 4, "mm"))
  })

  fig_forest <- wrap_plots(forest_panels, nrow = 1) +
    plot_annotation(
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )

  ggsave(here("Figures","fig_cross_study_params.tiff"),
         fig_forest, width = 12, height = 3.8,
         dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
  message("  Saved: Figures/fig_cross_study_params.tiff")
} else {
  message("  Skipping forest plot — summary CSVs not found")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Study 2 Condition Effects: γ vs. α (individual points + boxplot)
# Core finding: α invariant across conditions; γ drops dramatically in
# High-Status condition. Side-by-side panels for α (projection rate) and
# γ (ingroup bias).
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Study 2 Condition Effects (γ vs α) ──")

ind_s2 <- tryCatch(
  read.csv(here("Results","ind_diffs_s2_full.csv")),
  error = function(e) { message("  Missing ind_diffs_s2_full.csv"); NULL }
)

if (!is.null(ind_s2) && "outgroup" %in% names(ind_s2)) {

  cond_order  <- c("Not UCR", "UCLA", "CSU LA")
  cond_labels <- c("Not UCR" = "Negation", "UCLA" = "High-Status", "CSU LA" = "Low-Status")
  cond_cols   <- c("Not UCR" = "#E69F00", "UCLA" = "#0072B2", "CSU LA" = "#009E73")

  ind_s2$outgroup <- factor(ind_s2$outgroup, levels = cond_order)

  make_cond_panel <- function(df, y_var, y_label, y_lim = NULL) {
    # Compute condition means for crossbar
    means <- df |>
      dplyr::group_by(outgroup) |>
      dplyr::summarise(m = mean(.data[[y_var]], na.rm = TRUE), .groups = "drop")

    p <- ggplot(df, aes(x = outgroup, y = .data[[y_var]], colour = outgroup)) +
      geom_jitter(width = 0.18, size = 1.4, alpha = 0.55) +
      geom_boxplot(aes(fill = outgroup), alpha = 0.25, colour = "grey30",
                   width = 0.45, outlier.shape = NA, linewidth = 0.6) +
      scale_colour_manual(values = cond_cols, guide = "none") +
      scale_fill_manual(values = cond_cols, guide = "none") +
      scale_x_discrete(labels = cond_labels) +
      labs(x = NULL, y = y_label) +
      theme_dissert(base_size = 10) +
      theme(axis.text.x = element_text(size = 9))

    if (!is.null(y_lim)) p <- p + coord_cartesian(ylim = y_lim)
    p
  }

  p_alpha  <- make_cond_panel(ind_s2, "m",      "\u03b1 (Projection Rate)", c(0, 10))
  p_gamma  <- make_cond_panel(ind_s2, "bias",   "\u03b3 (Ingroup Bias)",    c(0, 1))
  p_lambda <- make_cond_panel(ind_s2, "lambda", "\u03bb (Generalization Sensitivity)", c(0, 5))

  fig_cond_s2 <- (p_alpha | p_gamma | p_lambda) +
    plot_annotation(
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )

  ggsave(here("Figures","fig_condition_effects_s2.tiff"),
         fig_cond_s2, width = 12, height = 4.5,
         dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
  message("  Saved: Figures/fig_condition_effects_s2.tiff")
} else {
  message("  Skipping condition effects figure (ind_diffs_s2_full.csv not found or missing outgroup column)")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Study 3 Condition Effects: α, γ, λ (Minority vs. Majority)
# Analogous to fig_condition_effects_s2 for racial group contrast.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Study 3 Condition Effects (\u03b1, \u03b3, \u03bb) ──")

ind_s3_cond <- tryCatch(
  read.csv(here("Results","ind_diffs_s3_full.csv")),
  error = function(e) { message("  Missing ind_diffs_s3_full.csv"); NULL }
)

if (!is.null(ind_s3_cond) && "condition" %in% names(ind_s3_cond)) {

  cond_order_s3  <- c("Minority", "Majority")
  cond_labels_s3 <- c("Minority" = "Racial Minority", "Majority" = "Racial Majority")
  cond_cols_s3   <- c("Minority" = "#CC79A7", "Majority" = "#0072B2")

  ind_s3_cond$condition <- factor(ind_s3_cond$condition, levels = cond_order_s3)

  make_cond_panel_s3 <- function(df, y_var, y_label, y_lim = NULL) {
    p <- ggplot(df, aes(x = condition, y = .data[[y_var]], colour = condition)) +
      geom_jitter(width = 0.18, size = 1.4, alpha = 0.55) +
      geom_boxplot(aes(fill = condition), alpha = 0.25, colour = "grey30",
                   width = 0.45, outlier.shape = NA, linewidth = 0.6) +
      scale_colour_manual(values = cond_cols_s3, guide = "none") +
      scale_fill_manual(values = cond_cols_s3, guide = "none") +
      scale_x_discrete(labels = cond_labels_s3) +
      labs(x = NULL, y = y_label) +
      theme_dissert(base_size = 10) +
      theme(axis.text.x = element_text(size = 9))
    if (!is.null(y_lim)) p <- p + coord_cartesian(ylim = y_lim)
    p
  }

  p_alpha_s3  <- make_cond_panel_s3(ind_s3_cond, "m",      "\u03b1 (Projection Rate)",         c(0, 10))
  p_gamma_s3  <- make_cond_panel_s3(ind_s3_cond, "bias",   "\u03b3 (Ingroup Bias)",             c(0, 1))
  p_lambda_s3 <- make_cond_panel_s3(ind_s3_cond, "lambda", "\u03bb (Generalization Sensitivity)", c(0, 5))

  fig_cond_s3 <- (p_alpha_s3 | p_gamma_s3 | p_lambda_s3) +
    plot_annotation(
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )

  ggsave(here("Figures","fig_condition_effects_s3.tiff"),
         fig_cond_s3, width = 12, height = 4.5,
         dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
  message("  Saved: Figures/fig_condition_effects_s3.tiff")
} else {
  message("  Skipping S3 condition effects figure (ind_diffs_s3_full.csv not found or missing condition column)")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Projection Rate (α) by Study and Condition (Pooled Estimates)
# Uses individual-level α from pooled model (params_ind_pooled_sym_lambda.csv)
# joined with condition info from ind_diffs files.
# Layout mirrors fig_generalization_gradients: 3 rows (one per study).
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Projection Rate (\u03b1) by Study/Condition (Pooled Estimates) ──")

pooled_ind <- tryCatch(
  read.csv(here("Results","params_ind_pooled_sym_lambda.csv")),
  error = function(e) { message("  Missing params_ind_pooled_sym_lambda.csv"); NULL }
)

ind_s1_pooled <- tryCatch(read.csv(here("Results","ind_diffs_s1_full.csv")), error = function(e) NULL)
ind_s2_pooled <- tryCatch(read.csv(here("Results","ind_diffs_s2_full.csv")), error = function(e) NULL)
ind_s3_pooled <- tryCatch(read.csv(here("Results","ind_diffs_s3_full.csv")), error = function(e) NULL)

if (!is.null(pooled_ind) && !is.null(ind_s1_pooled) && !is.null(ind_s2_pooled) && !is.null(ind_s3_pooled)) {

  n1 <- nrow(ind_s1_pooled)
  n2 <- nrow(ind_s2_pooled)
  n3 <- nrow(ind_s3_pooled)

  m_rows    <- pooled_ind[startsWith(pooled_ind$variable, "\"m[") |
                           startsWith(pooled_ind$variable, "m["), ]
  bias_rows <- pooled_ind[startsWith(pooled_ind$variable, "\"bias[") |
                           startsWith(pooled_ind$variable, "bias["), ]
  lam_rows  <- pooled_ind[startsWith(pooled_ind$variable, "\"lambda[") |
                           startsWith(pooled_ind$variable, "lambda["), ]

  if (nrow(m_rows) == 507) {
    pooled_alpha_df <- data.frame(
      alpha     = m_rows$median,
      gamma     = bias_rows$median,
      lambda    = lam_rows$median,
      study     = c(rep("Study 1", n1), rep("Study 2", n2), rep("Study 3", n3)),
      condition = c(
        rep("Minimal Groups", n1),
        dplyr::recode(ind_s2_pooled$outgroup,
                      "Not UCR" = "Negation", "UCLA" = "High-Status", "CSU LA" = "Low-Status"),
        ind_s3_pooled$condition
      )
    )

    study_cond_cols <- c(
      "Minimal Groups" = "#2B5C8A",
      "Negation"       = "#E69F00",
      "High-Status"    = "#0072B2",
      "Low-Status"     = "#009E73",
      "Minority"       = "#CC79A7",
      "Majority"       = "#0072B2"
    )

    make_alpha_panel <- function(df, study_label, cond_levels, cond_cols_sub) {
      d <- df[df$study == study_label, ]
      d$condition <- factor(d$condition, levels = cond_levels)
      means <- d |>
        dplyr::group_by(condition) |>
        dplyr::summarise(m = mean(alpha, na.rm = TRUE), .groups = "drop")
      ggplot(d, aes(x = condition, y = alpha, colour = condition)) +
        geom_jitter(width = 0.15, size = 1.2, alpha = 0.45) +
        geom_boxplot(aes(fill = condition), alpha = 0.22, colour = "grey30",
                     width = 0.40, outlier.shape = NA, linewidth = 0.55) +
        scale_colour_manual(values = cond_cols_sub, guide = "none") +
        scale_fill_manual(values = cond_cols_sub, guide = "none") +
        scale_y_continuous(name = "\u03b1 (Projection Rate)", limits = c(0, 10)) +
        labs(x = NULL, title = study_label) +
        theme_dissert(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, size = 9.5, face = "bold"),
              axis.text.x = element_text(size = 8.5))
    }

    pa1 <- make_alpha_panel(pooled_alpha_df, "Study 1",
                            c("Minimal Groups"),
                            c("Minimal Groups" = "#2B5C8A"))
    pa2 <- make_alpha_panel(pooled_alpha_df, "Study 2",
                            c("Negation","High-Status","Low-Status"),
                            c("Negation" = "#E69F00","High-Status" = "#0072B2","Low-Status" = "#009E73"))
    pa3 <- make_alpha_panel(pooled_alpha_df, "Study 3",
                            c("Minority","Majority"),
                            c("Minority" = "#CC79A7","Majority" = "#0072B2"))

    fig_alpha_cond <- (pa1 | pa2 | pa3) +
      plot_layout(widths = c(1, 3, 2)) +
      plot_annotation(
        theme = theme(plot.background = element_rect(fill = "white", colour = NA))
      )

    ggsave(here("Figures","fig_alpha_by_condition.tiff"),
           fig_alpha_cond, width = 12, height = 4.5,
           dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
    message("  Saved: Figures/fig_alpha_by_condition.tiff")

    # ── Combined α + Generalization Gradient Figure ────────────────────────────
    # Two-column (α slope | gradient), three-row (study) grid.
    # Left column: ingroup projection slope G_in(E) = plogis(α*(E-4)) per condition.
    #   Individual participant curves (faint) + condition mean (bold, colored).
    # Right column: S^λ generalization gradient per condition (no stats annotation).
    # Emulates fig07_marginal_effects. Consistent color coding throughout.
    message("\n── Combined \u03b1 + Generalization Gradient Figure ──")

    make_alpha_slope_panel <- function(df_study, study_label,
                                       cond_levels, cond_cols_sub,
                                       cond_labels_sub, cond_ltys_sub,
                                       legend_pos = c(0.03, 0.97)) {
      E_seq_loc <- seq(1, 7, length.out = 200)
      df_study$condition <- factor(df_study$condition, levels = cond_levels)

      # Individual participant curves
      ind_curves <- do.call(rbind, lapply(seq_len(nrow(df_study)), function(i) {
        data.frame(
          E         = E_seq_loc,
          G_in      = plogis(df_study$alpha[i] * (E_seq_loc - 4)),
          condition = df_study$condition[i],
          subj_i    = i
        )
      }))
      ind_curves$condition <- factor(ind_curves$condition, levels = cond_levels)

      # Condition mean curves
      mean_curves <- do.call(rbind, lapply(cond_levels, function(cond) {
        alpha_mean <- mean(df_study$alpha[df_study$condition == cond], na.rm = TRUE)
        data.frame(
          E         = E_seq_loc,
          G_in      = plogis(alpha_mean * (E_seq_loc - 4)),
          condition = factor(cond, levels = cond_levels)
        )
      }))

      ggplot() +
        geom_line(data = ind_curves,
                  aes(x = E, y = G_in, group = interaction(condition, subj_i),
                      colour = condition),
                  alpha = 0.06, linewidth = 0.2) +
        geom_vline(xintercept = 4, linetype = "dashed", color = "gray60", linewidth = 0.4) +
        geom_line(data = mean_curves,
                  aes(x = E, y = G_in, colour = condition, linetype = condition),
                  linewidth = 1.1) +
        scale_colour_manual(values = cond_cols_sub, labels = cond_labels_sub, name = NULL) +
        scale_linetype_manual(values = cond_ltys_sub, labels = cond_labels_sub, name = NULL) +
        scale_x_continuous(name = "Self-Rating on Trait",
                           breaks = c(1, 4, 7),
                           labels = c("1\n(Not at all)", "4\n(Neutral)", "7\n(Extremely)")) +
        scale_y_continuous(name = "Ingroup Projection Weight",
                           limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        labs(title = study_label) +
        theme_dissert(base_size = 9.5) +
        theme(legend.position        = legend_pos,
              legend.justification   = c("left", "top"),
              legend.text            = element_text(size = 7.5),
              legend.key.size        = unit(0.9, "lines"),
              plot.title             = element_text(face = "bold", size = 8.5, hjust = 0,
                                                    margin = margin(b = 1, t = 1)))
    }

    d_s1_c <- pooled_alpha_df[pooled_alpha_df$study == "Study 1", ]
    d_s2_c <- pooled_alpha_df[pooled_alpha_df$study == "Study 2", ]
    d_s3_c <- pooled_alpha_df[pooled_alpha_df$study == "Study 3", ]

    # S1: single condition
    pa_s1 <- make_alpha_slope_panel(
      d_s1_c, "Study 1 (Minimal Groups)",
      cond_levels    = "Minimal Groups",
      cond_cols_sub  = c("Minimal Groups" = "black"),
      cond_labels_sub = c("Minimal Groups" = "Minimal Groups"),
      cond_ltys_sub  = c("Minimal Groups" = "solid"),
      legend_pos     = "none"
    )
    # S2: three conditions — colors match s2_cols (original code keys recoded to display labels)
    pa_s2 <- make_alpha_slope_panel(
      d_s2_c, "Study 2 (University Groups)",
      cond_levels    = c("Negation", "High-Status", "Low-Status"),
      cond_cols_sub  = c("Negation" = "#E69F00", "High-Status" = "#0072B2", "Low-Status" = "#009E73"),
      cond_labels_sub = c("Negation" = "Negation", "High-Status" = "High-Status", "Low-Status" = "Low-Status"),
      cond_ltys_sub  = c("Negation" = "solid", "High-Status" = "dashed", "Low-Status" = "dotdash")
    )
    # S3: two conditions
    pa_s3 <- make_alpha_slope_panel(
      d_s3_c, "Study 3 (Racial Groups)",
      cond_levels    = c("Minority", "Majority"),
      cond_cols_sub  = c("Minority" = "#CC79A7", "Majority" = "#0072B2"),
      cond_labels_sub = c("Minority" = "Racial Minority", "Majority" = "Racial Majority"),
      cond_ltys_sub  = c("Minority" = "solid", "Majority" = "dashed")
    )

    # Gradient panels: rebuild S1 without the λ annotation (inconsistent in combined figure);
    # pg2/pg3 are now correctly built with fixed condition labels and subject grouping.
    make_grad_panel_s1_comb <- function() {
      if (is.null(summary_s1)) return(ggplot() + theme_void())
      lambda_row <- summary_s1 |> dplyr::filter(grepl("^\"?mu_pr\\[3\\]", variable))
      if (nrow(lambda_row) == 0) return(ggplot() + theme_void())
      lam_med <- pnorm(as.numeric(lambda_row$median)) * 5
      lam_lo  <- pnorm(as.numeric(lambda_row$q5))     * 5
      lam_hi  <- pnorm(as.numeric(lambda_row$q95))    * 5
      grad_s1_c <- data.frame(S = S_seq, g_med = S_seq ^ lam_med,
                               g_lo = S_seq ^ lam_hi, g_hi = S_seq ^ lam_lo)
      params_s1_c <- tryCatch(read.csv(here("Results","params_ind_s1_sym_lambda.csv")),
                               error = function(e) NULL)
      ind_lam_s1_c <- if (!is.null(params_s1_c)) {
        params_s1_c |> dplyr::filter(grepl("^\"?lambda\\[", variable)) |> dplyr::pull(median)
      } else { numeric(0) }
      p <- ggplot()
      if (length(ind_lam_s1_c) > 0) {
        ind_df_c <- do.call(rbind, lapply(seq_along(ind_lam_s1_c), function(i)
          data.frame(S = S_seq, g = S_seq ^ ind_lam_s1_c[i], subj_i = i)))
        p <- p + geom_line(data = ind_df_c, aes(x = S, y = g, group = subj_i),
                           colour = s1_col, alpha = 0.06, linewidth = 0.25)
      }
      p +
        geom_ribbon(data = grad_s1_c, aes(x = S, ymin = g_lo, ymax = g_hi),
                    fill = s1_col, alpha = 0.20) +
        geom_line(data = grad_s1_c, aes(x = S, y = g_med),
                  colour = s1_col, linewidth = 1.1) +
        scale_x_continuous(name = "Semantic Similarity (Dice)", breaks = seq(0, 1, 0.25)) +
        scale_y_continuous(name = expression("Generalization Weight " ~ (S^lambda)),
                           limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        labs(title = expression("Generalization Gradient " ~ (S^lambda))) +
        theme_dissert(base_size = 10) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", size = 8.5, hjust = 0,
                                        margin = margin(b = 1, t = 1)))
    }

    pg1_comb <- make_grad_panel_s1_comb()
    pg2_comb <- pg2 + labs(title = NULL)
    pg3_comb <- pg3 + labs(title = NULL)

    fig_combined_alpha_grad <- (pa_s1 | pg1_comb) /
                                (pa_s2 | pg2_comb) /
                                (pa_s3 | pg3_comb) +
      plot_annotation(
        theme = theme(plot.background = element_rect(fill = "white", colour = NA))
      )

    ggsave(here("Figures","fig_combined_alpha_gradient.tiff"),
           fig_combined_alpha_grad, width = 12, height = 9,
           dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
    message("  Saved: Figures/fig_combined_alpha_gradient.tiff")

  } else {
    message(sprintf("  Unexpected row count in pooled m params (%d vs 507 expected); skipping.", nrow(m_rows)))
  }
} else {
  message("  Skipping pooled alpha figure — required CSVs not found")
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — MCR × Individual Differences (all 3 studies, SING.Ind)
# MCR × SING.Ind replicates across all three studies (FDR-significant each).
# 3-panel scatter (one per study) with regression line and r annotation.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── MCR × Individual Differences (all 3 studies) ──")

load_ind <- function(path, study_label, condition_col = NULL) {
  if (!file.exists(path)) { message("  Missing: ", path); return(NULL) }
  df <- read.csv(path)
  df$study <- study_label
  df
}

id_s1 <- load_ind(here("Results","ind_diffs_s1_full.csv"), "Study 1\n(Minimal Groups)")
id_s2 <- load_ind(here("Results","ind_diffs_s2_full.csv"), "Study 2\n(University Groups)")
id_s3 <- load_ind(here("Results","ind_diffs_s3_full.csv"), "Study 3\n(Racial Groups)")

make_mcr_scatter <- function(df, x_var, x_label, study_col) {
  if (is.null(df) || !x_var %in% names(df) || !"subject_mcr" %in% names(df)) {
    return(ggplot() + theme_void() +
             annotate("rect", xmin=0,xmax=1,ymin=0,ymax=1,fill="grey95",color="grey70") +
             annotate("text", x=0.5, y=0.5,
                      label = paste0(unique(df$study), "\n[Pending]"),
                      hjust=0.5, size=3, color="grey50"))
  }
  d_clean <- df[!is.na(df[[x_var]]) & !is.na(df$subject_mcr), ]
  ct <- cor.test(d_clean[[x_var]], d_clean$subject_mcr)
  r  <- round(ct$estimate, 2)
  p  <- ct$p.value
  p_label <- if (p < .001) "p < .001" else sprintf("p = %.3f", p)
  r_label <- sprintf("r = %s\n%s", formatC(r, format="f", digits=2), p_label)

  ggplot(d_clean, aes_string(x = x_var, y = "subject_mcr")) +
    geom_point(colour = study_col, size = 1.8, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, colour = study_col,
                fill = study_col, alpha = 0.15, linewidth = 0.9) +
    annotate("label", x = Inf, y = Inf, hjust = 1.08, vjust = 1.35,
             label = r_label, size = 3, fill = "white", colour = "grey20",
             label.size = 0.3) +
    labs(x = x_label, y = "Metacontrast Ratio (MCR)",
         title = unique(df$study)) +
    theme_dissert(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 9.5))
}

study_cols_vec <- c(
  "Study 1\n(Minimal Groups)"    = "#2B5C8A",
  "Study 2\n(University Groups)" = "#E69F00",
  "Study 3\n(Racial Groups)"     = "#CC79A7"
)

x_var   <- "SING.Ind"
x_label <- "Independent Self-Construal (SING.Ind)"

p_mcr1 <- make_mcr_scatter(id_s1, x_var, x_label, study_cols_vec["Study 1\n(Minimal Groups)"])
p_mcr2 <- make_mcr_scatter(id_s2, x_var, x_label, study_cols_vec["Study 2\n(University Groups)"])
p_mcr3 <- make_mcr_scatter(id_s3, x_var, x_label, study_cols_vec["Study 3\n(Racial Groups)"])

fig_mcr_pers <- (p_mcr1 | p_mcr2 | p_mcr3) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

message("  Skipped: fig_mcr_personality.tiff (removed from manuscript)")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Trait Network (Supplementary Materials)
# Visualizes the semantic adjacency graph using igraph force-directed layout.
# Nodes = traits, edges = adjacency. Saved for Supplementary reference.
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Trait Network (Supplementary) ──")

has_igraph <- requireNamespace("igraph", quietly = TRUE)
has_ggraph <- requireNamespace("ggraph", quietly = TRUE)

if (has_igraph) {
  library(igraph)

  adj_path <- here("Pooled","input","adjacencyMatrix_p.csv")
  if (file.exists(adj_path)) {
    adj_mat   <- as.matrix(read.csv(adj_path, check.names = FALSE))
    # Use first 80 most-connected traits for legibility
    deg       <- rowSums(adj_mat)
    top_idx   <- order(deg, decreasing = TRUE)[1:min(80, nrow(adj_mat))]
    adj_sub   <- adj_mat[top_idx, top_idx]
    g         <- graph_from_adjacency_matrix(adj_sub, mode = "undirected", diag = FALSE)
    V(g)$degree <- degree(g)

    if (has_ggraph) {
      library(ggraph)
      set.seed(42)
      fig_net <- ggraph(g, layout = "fr") +
        geom_edge_link(colour = "grey75", alpha = 0.5, linewidth = 0.3) +
        geom_node_point(aes(size = degree, colour = degree)) +
        scale_colour_gradient(low = "#c7e9c0", high = "#006d2c",
                              name = "Degree") +
        scale_size(range = c(1.5, 5), guide = "none") +
        labs(
          title    = "Semantic Trait Network (Top 80 Traits by Connectivity)",
          subtitle = "Nodes = personality traits; edges = adjacency in SAGE database. Force-directed layout (Fruchterman-Reingold)."
        ) +
        theme_void(base_size = 9) +
        theme(
          plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.subtitle   = element_text(size = 8, colour = "grey40", hjust = 0.5),
          legend.position = "right",
          plot.background = element_rect(fill = "white", colour = NA)
        )
    } else {
      # Fallback: base R igraph plot via ggplot raster
      message("  ggraph not installed — using base igraph plot for trait network")
      tiff_out <- here("Figures","fig01_trait_network.tiff")
      tiff(tiff_out, width = 8, height = 8, units = TIFF_UNITS, res = TIFF_DPI,
           compression = "lzw")
      set.seed(42)
      plot(g, vertex.size = 4, vertex.label = NA, edge.color = "grey70",
           vertex.color = "#2B5C8A", main = "Semantic Trait Network (Top 80 Traits)")
      dev.off()
      message("  Saved (base igraph): Figures/fig01_trait_network.tiff")
      fig_net <- NULL
    }

    if (!is.null(fig_net)) {
      ggsave(here("Figures","fig01_trait_network.tiff"),
             fig_net, width = 8, height = 7,
             dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")
      message("  Saved: Figures/fig01_trait_network.tiff")
    }
  } else {
    message("  Missing: Pooled/input/adjacencyMatrix_p.csv")
  }
} else {
  message("  igraph not installed — skipping trait network figure")
}


message("\n── Done. Run marginal_effects_s*.R scripts to regenerate prediction CSVs,")
message("   then re-run this script to update fig07_marginal_effects.tiff.")
message("   Run run_model_comparison_s*.R to regenerate ELPD raincloud figures.\n")

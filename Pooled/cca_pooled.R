# cca_pooled.R
# Pooled Canonical Correlation Analysis: model parameters × personality
# All three studies combined (N ≈ 503), scales consistent across studies.
#
# Key design choice: Z-score X and Y within each study before pooling, so
# that study-level mean/variance differences do not drive the canonical
# solution. The CCA then captures individual-difference structure that
# generalises across studies.
#
# Set X (model-derived): m, bias, lambda, w, subject_mcr, delta_elpd
# Set Y (personality):   DS, Proto, SCC, SI, RSE, NTB, NFC, SING.Ind, SING.Inter
#
# Prerequisites:
#   Results/ind_diffs_s1_full.csv           (from loo_delta_indiff_s1.R)
#   Results/ind_diffs_s2_full_enriched.csv  (from param_comparison_conditions_s2.R)
#   Results/ind_diffs_s3_full_enriched.csv  (from param_comparison_conditions_s3.R)
#
# Output:
#   Results/cca_pooled_correlations.csv
#   Results/cca_pooled_loadings.csv
#   Results/cca_pooled_biplot.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── Load individual differences per study ─────────────────────────────────────
s1 <- read.csv(here("Results", "ind_diffs_s1_full.csv")) %>%
  mutate(study = "Study 1 (Minimal Group)", condition = "Minimal Group")

s2_file <- here("Results", "ind_diffs_s2_full_enriched.csv")
s2 <- read.csv(if (file.exists(s2_file)) s2_file else
               here("Results", "ind_diffs_s2_full.csv")) %>%
  mutate(
    study = "Study 2 (University Status)",
    condition = case_when(
      outgroup == "UCLA"    ~ "High-Status",
      outgroup == "Not UCR" ~ "Negation",
      outgroup == "CSU LA"  ~ "Low-Status",
      TRUE                  ~ as.character(outgroup)
    )
  )

s3_file <- here("Results", "ind_diffs_s3_full_enriched.csv")
s3 <- read.csv(if (file.exists(s3_file)) s3_file else
               here("Results", "ind_diffs_s3_full.csv")) %>%
  mutate(study = "Study 3 (Racial Groups)")

# ── Stack and label ───────────────────────────────────────────────────────────
possible_x <- c("m", "lambda", "bias", "w", "subject_mcr")
possible_y <- c("DS", "Proto", "SCC", "SI", "RSE", "NTB", "NFC", "SING.Ind", "SING.Inter")

keep_cols <- c("subID", "study", "condition",
               intersect(possible_x, c(names(s1), names(s2), names(s3))),
               intersect(possible_y, c(names(s1), names(s2), names(s3))))

df <- bind_rows(
  s1 %>% select(any_of(keep_cols)),
  s2 %>% select(any_of(keep_cols)),
  s3 %>% select(any_of(keep_cols))
)

x_vars <- intersect(possible_x, names(df))
y_vars <- intersect(possible_y, names(df))

cat(sprintf("Loaded: S1 N=%d, S2 N=%d, S3 N=%d\n",
            sum(df$study == "Study 1 (Minimal Group)"),
            sum(df$study == "Study 2 (University Status)"),
            sum(df$study == "Study 3 (Racial Groups)")))
cat(sprintf("X vars (%d): %s\n", length(x_vars), paste(x_vars, collapse = ", ")))
cat(sprintf("Y vars (%d): %s\n", length(y_vars), paste(y_vars, collapse = ", ")))

# ── Drop incomplete cases ─────────────────────────────────────────────────────
df_cc <- df %>%
  select(subID, study, condition, all_of(c(x_vars, y_vars))) %>%
  drop_na(all_of(c(x_vars, y_vars)))
message(sprintf("N complete cases: %d", nrow(df_cc)))

# ── Z-score within study (critical: removes study-level shifts) ───────────────
df_cc <- df_cc %>%
  group_by(study) %>%
  mutate(across(all_of(c(x_vars, y_vars)), ~ as.numeric(scale(.)))) %>%
  ungroup()

X <- as.matrix(df_cc[, x_vars])
Y <- as.matrix(df_cc[, y_vars])

# ── CCA ───────────────────────────────────────────────────────────────────────
cc    <- cancor(X, Y)
n_can <- min(ncol(X), ncol(Y))

cat("\n═══════════════════════════════════════\n")
cat("Pooled CCA — Canonical Correlations:\n")
cat("═══════════════════════════════════════\n")
for (i in seq_len(n_can)) cat(sprintf("  CV%d:  r = %.3f\n", i, cc$cor[i]))

# ── Permutation test ──────────────────────────────────────────────────────────
set.seed(2026)
n_perm    <- 5000
perm_cors <- replicate(n_perm, { Xp <- X[sample(nrow(X)), ]; cancor(Xp, Y)$cor })
p_perm    <- sapply(seq_len(n_can), function(i) mean(perm_cors[i, ] >= cc$cor[i]))

cat("\nPermutation p-values (B = 5000):\n")
for (i in seq_len(n_can)) {
  cat(sprintf("  CV%d:  r = %.3f  p_perm = %.3f\n", i, cc$cor[i], p_perm[i]))
}

# ── Loadings ──────────────────────────────────────────────────────────────────
scores_X <- X %*% cc$xcoef
scores_Y <- Y %*% cc$ycoef
load_X   <- cor(X, scores_X)
load_Y   <- cor(Y, scores_Y)

n_print <- min(3, n_can)
cat(sprintf("\n─── Loadings on CV1–CV%d (model parameters) ───\n", n_print))
print(round(load_X[, 1:n_print, drop = FALSE], 3))
cat(sprintf("\n─── Loadings on CV1–CV%d (personality) ───\n", n_print))
print(round(load_Y[, 1:n_print, drop = FALSE], 3))

# ── Save results ──────────────────────────────────────────────────────────────
cv_res <- data.frame(
  cv     = paste0("CV", seq_len(n_can)),
  r      = round(cc$cor, 4),
  p_perm = round(p_perm, 4)
)

load_X_df <- as.data.frame(load_X)
colnames(load_X_df) <- paste0("CV", seq_len(n_can))
load_X_df$variable <- x_vars; load_X_df$set <- "model"

load_Y_df <- as.data.frame(load_Y)
colnames(load_Y_df) <- paste0("CV", seq_len(n_can))
load_Y_df$variable <- y_vars; load_Y_df$set <- "personality"

write.csv(cv_res,
          here("Results", "cca_pooled_correlations.csv"), row.names = FALSE)
write.csv(bind_rows(load_X_df, load_Y_df),
          here("Results", "cca_pooled_loadings.csv"),     row.names = FALSE)

# ── Biplot ────────────────────────────────────────────────────────────────────
label_map <- c(
  m           = "\u03b1 (projection)",
  lambda      = "\u03bb (generalization)",
  bias        = "\u03b3 (ingroup bias)",
  w           = "w (lapse)",
  subject_mcr = "MCR",
  delta_elpd  = "\u0394ELPD",
  DS          = "Dial. Self",
  Proto       = "Prototypicality",
  SCC         = "SCC",
  SI          = "Social ID",
  RSE         = "Self-Esteem",
  NTB         = "Need to Belong",
  NFC         = "Need for Cog.",
  SING.Ind    = "SING-Ind",
  SING.Inter  = "SING-Inter"
)

ld <- bind_rows(
  data.frame(variable = x_vars, CV1 = load_X[, 1],
             CV2 = if (n_can >= 2) load_X[, 2] else 0, set = "Model (X)"),
  data.frame(variable = y_vars, CV1 = load_Y[, 1],
             CV2 = if (n_can >= 2) load_Y[, 2] else 0, set = "Personality (Y)")
)
ld$label <- ifelse(ld$variable %in% names(label_map), label_map[ld$variable], ld$variable)

r1_lab <- sprintf("CV1  (r = %.3f, p = %.3f)", cc$cor[1], p_perm[1])
r2_lab <- if (n_can >= 2) sprintf("CV2  (r = %.3f, p = %.3f)", cc$cor[2], p_perm[2]) else "CV2"

library(ggplot2)
p_biplot <- ggplot(ld, aes(x = CV1, y = CV2, colour = set, label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_segment(aes(x = 0, y = 0, xend = CV1, yend = CV2),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               linewidth = 0.6, alpha = 0.7) +
  scale_colour_manual(values = c("Model (X)" = "#2B5C8A", "Personality (Y)" = "#CC79A7"),
                      name = NULL) +
  labs(
    title    = sprintf("Pooled CCA Biplot — Studies 1–3 (N = %d)", nrow(df_cc)),
    subtitle = "Variable loadings on first two canonical variates (within-study Z-scored)",
    x = r1_lab, y = r2_lab,
    caption = paste0(
      "Model set: ", paste(x_vars, collapse = ", "), ".\n",
      "Personality set: ", paste(y_vars, collapse = ", "), ".\n",
      "X and Y Z-scored within study before pooling to remove study-level shifts."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.background  = element_rect(fill = "white", colour = NA))

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_biplot <- p_biplot +
    ggrepel::geom_text_repel(size = 3, fontface = "bold", max.overlaps = 20,
                             segment.size = 0.25, box.padding = 0.4)
} else {
  p_biplot <- p_biplot + geom_text(size = 2.8, vjust = -0.5)
}

ggsave(here("Results", "cca_pooled_biplot.png"), p_biplot, width = 8, height = 7, dpi = 200)

message("Saved: Results/cca_pooled_biplot.png")
message("Saved: Results/cca_pooled_correlations.csv, Results/cca_pooled_loadings.csv")
message(sprintf(
  "\nPooled CCA complete: N = %d subjects (%d + %d + %d), %d + %d variables.",
  nrow(df_cc),
  sum(df_cc$study == "Study 1 (Minimal Group)"),
  sum(df_cc$study == "Study 2 (University Status)"),
  sum(df_cc$study == "Study 3 (Racial Groups)"),
  length(x_vars), length(y_vars)
))

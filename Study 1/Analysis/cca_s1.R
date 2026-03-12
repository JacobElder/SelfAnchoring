# cca_s1.R
# Canonical Correlation Analysis: model parameters × personality (Study 1)
# Parallel structure to Study 2/Analysis/cca_s2.R and Study 3/Analysis/cca_s3.R
#
# Set X (model-derived): m, lambda, bias, w, subject_mcr, delta_elpd, groupHomoph
# Set Y (personality):   DS, Proto, SCC, SI, RSE, NTB, NFC, SING.Ind, SING.Inter
#
# N = 61; exploratory — not preregistered.
#
# Prerequisites:
#   - ind_diffs_s1_full.csv (from loo_delta_indiff_s1.R)
#
# Output:
#   Results/cca_s1_biplot.png
#   Results/cca_s1_correlations.csv
#   Results/cca_s1_loadings.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── Load and merge data ────────────────────────────────────────────────────────
ind <- read.csv(here("Results", "ind_diffs_s1_full.csv"))

# groupHomoph: one value per subject; grab from fullTest
ft <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) |>
  filter(!is.na(ingChoiceN)) |>
  group_by(subID) |>
  slice(1) |>
  ungroup() |>
  select(subID, groupHomoph)

df <- left_join(ind, ft, by = "subID")
message(sprintf("N subjects: %d", nrow(df)))

# ── Define variable sets ───────────────────────────────────────────────────────
possible_x <- c("m", "lambda", "bias", "w", "subject_mcr", "delta_elpd", "groupHomoph")
x_vars <- intersect(possible_x, names(df))

possible_y <- c("DS", "Proto", "SCC", "SI", "RSE", "NTB", "NFC", "SING.Ind", "SING.Inter")
y_vars <- intersect(possible_y, names(df))

cat(sprintf("X vars (%d): %s\n", length(x_vars), paste(x_vars, collapse = ", ")))
cat(sprintf("Y vars (%d): %s\n", length(y_vars), paste(y_vars, collapse = ", ")))

df_cc <- df |>
  select(all_of(c(x_vars, y_vars))) |>
  drop_na()

message(sprintf("N complete cases: %d", nrow(df_cc)))
if (nrow(df_cc) < 20) stop("Too few complete cases for CCA.")

X <- scale(df_cc[, x_vars])
Y <- scale(df_cc[, y_vars])

# ── CCA ───────────────────────────────────────────────────────────────────────
cc    <- cancor(X, Y)
n_can <- min(ncol(X), ncol(Y))

cat("\n═══════════════════════════════════════\n")
cat("Canonical Correlations:\n")
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

cat("\n─── Loadings on CV1–CV3 (model parameters) ───\n")
print(round(load_X[, 1:min(3, n_can), drop = FALSE], 3))
cat("\n─── Loadings on CV1–CV3 (personality) ───\n")
print(round(load_Y[, 1:min(3, n_can), drop = FALSE], 3))

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

write.csv(cv_res,                          here("Results", "cca_s1_correlations.csv"), row.names = FALSE)
write.csv(bind_rows(load_X_df, load_Y_df), here("Results", "cca_s1_loadings.csv"),     row.names = FALSE)

# ── Biplot ────────────────────────────────────────────────────────────────────
label_map <- c(
  m           = "\u03b1 (projection)",
  lambda      = "\u03bb (generalization)",
  bias        = "\u03b3 (ingroup bias)",
  w           = "w (lapse)",
  subject_mcr = "MCR",
  delta_elpd  = "\u0394ELPD",
  groupHomoph = "Trait Segregation",
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
  data.frame(variable = x_vars, CV1 = load_X[,1], CV2 = load_X[,2], set = "Model (X)"),
  data.frame(variable = y_vars, CV1 = load_Y[,1], CV2 = load_Y[,2], set = "Personality (Y)")
)
ld$label <- ifelse(ld$variable %in% names(label_map), label_map[ld$variable], ld$variable)

r1_lab <- sprintf("CV1  (r = %.3f, p = %.3f)", cc$cor[1], p_perm[1])
r2_lab <- sprintf("CV2  (r = %.3f, p = %.3f)", cc$cor[2], p_perm[2])

p_biplot <- ggplot(ld, aes(x = CV1, y = CV2, colour = set, label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_segment(aes(x = 0, y = 0, xend = CV1, yend = CV2),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               linewidth = 0.6, alpha = 0.7) +
  scale_colour_manual(values = c("Model (X)" = "#2B5C8A", "Personality (Y)" = "#CC79A7"),
                      name = NULL) +
  labs(
    title    = sprintf("CCA Biplot — Study 1, Minimal Groups (N = %d)", nrow(df_cc)),
    subtitle = "Variable loadings on first two canonical variates",
    x = r1_lab, y = r2_lab,
    caption = paste0(
      "Model set (blue): ", paste(x_vars, collapse = ", "), ".\n",
      "Personality set (pink): ", paste(y_vars, collapse = ", "), ".\n",
      "N = ", nrow(df_cc), " complete cases."
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

ggsave(here("Results", "cca_s1_biplot.png"), p_biplot, width = 8, height = 7, dpi = 200)

message("Saved: Results/cca_s1_biplot.png")
message("Saved: Results/cca_s1_correlations.csv, Results/cca_s1_loadings.csv")
message(sprintf("\nNote: N = %d; p = %d + %d = %d variables. CCA exploratory; first 1-2 CVs interpretable.",
                nrow(df_cc), length(x_vars), length(y_vars), length(x_vars) + length(y_vars)))

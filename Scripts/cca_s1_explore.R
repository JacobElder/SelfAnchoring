# cca_s1_explore.R
# Canonical Correlation Analysis: model parameters × personality (Study 1 only)
#
# Set X (model-derived): m, lambda, bias, w, subject_mcr, groupHomoph
# Set Y (personality):   DS, Proto, SCC, SI, RSE, NTB, NFC, SING.Ind, SING.Inter
#
# N = 61; exploratory only — not added to paper.
# Output: console loadings + Scripts/cca_s1_biplot.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── Load and merge data ────────────────────────────────────────────────────────
ind <- read.csv(here("Results", "ind_diffs_s1_full.csv"))

# groupHomoph: one value per subject; grab from fullTest
ft <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) |>
  group_by(subID) |>
  slice(1) |>
  ungroup() |>
  select(subID, groupHomoph)

df <- left_join(ind, ft, by = "subID")

message(sprintf("N subjects: %d", nrow(df)))

# ── Define variable sets ───────────────────────────────────────────────────────
x_vars <- c("m", "lambda", "bias", "w", "subject_mcr", "groupHomoph")
y_vars <- c("DS", "Proto", "SCC", "SI", "RSE", "NTB", "NFC", "SING.Ind", "SING.Inter")

df_cc <- df |>
  select(all_of(c(x_vars, y_vars))) |>
  drop_na()

message(sprintf("N complete cases: %d", nrow(df_cc)))

X <- scale(df_cc[, x_vars])
Y <- scale(df_cc[, y_vars])

# ── CCA via base cancor() ─────────────────────────────────────────────────────
cc <- cancor(X, Y)

n_can <- min(ncol(X), ncol(Y))  # = 7

cat("\n═══════════════════════════════════════\n")
cat("Canonical Correlations:\n")
cat("═══════════════════════════════════════\n")
for (i in seq_len(n_can)) {
  cat(sprintf("  CV%d:  r = %.3f\n", i, cc$cor[i]))
}

# ── Permutation test for each canonical correlation ───────────────────────────
set.seed(2026)
n_perm <- 5000

perm_cors <- replicate(n_perm, {
  Xp <- X[sample(nrow(X)), ]
  cancor(Xp, Y)$cor
})
# perm_cors is n_can × n_perm matrix

p_perm <- sapply(seq_len(n_can), function(i) {
  mean(perm_cors[i, ] >= cc$cor[i])
})

cat("\nPermutation p-values (B = 5000):\n")
for (i in seq_len(n_can)) {
  cat(sprintf("  CV%d:  r = %.3f  p_perm = %.3f\n", i, cc$cor[i], p_perm[i]))
}

# ── Loadings: correlation of original variables with canonical scores ──────────
# Canonical scores
scores_X <- X %*% cc$xcoef          # n × n_can
scores_Y <- Y %*% cc$ycoef          # n × n_can

load_X <- cor(X, scores_X)           # x_vars × n_can
load_Y <- cor(Y, scores_Y)           # y_vars × n_can

cat("\n─── Loadings on CV1 (model parameters) ───\n")
print(round(load_X[, 1:min(3, n_can), drop = FALSE], 3))

cat("\n─── Loadings on CV1 (personality) ───\n")
print(round(load_Y[, 1:min(3, n_can), drop = FALSE], 3))

# ── Biplot of variable loadings on CV1 × CV2 ─────────────────────────────────
ld <- bind_rows(
  data.frame(
    variable = rownames(load_X),
    CV1      = load_X[, 1],
    CV2      = load_X[, 2],
    set      = "Model (X)"
  ),
  data.frame(
    variable = rownames(load_Y),
    CV1      = load_Y[, 1],
    CV2      = load_Y[, 2],
    set      = "Personality (Y)"
  )
)

# Nicer labels
label_map <- c(
  m            = "α (projection)",
  lambda       = "λ (generalization)",
  bias         = "γ (ingroup bias)",
  w            = "w (lapse)",
  subject_mcr  = "MCR",
  delta_elpd   = "ΔELPD",
  groupHomoph  = "Trait Segregation",
  DS           = "Dial. Self",
  Proto        = "Prototypicality",
  SCC          = "SCC",
  SI           = "Social ID",
  RSE          = "Self-Esteem",
  NTB          = "Need to Belong",
  NFC          = "Need for Cog.",
  SING.Ind     = "SING-Ind",
  SING.Inter   = "SING-Inter"
)
ld$label <- label_map[ld$variable]

r1_lab <- sprintf("CV1  (r = %.3f, p = %.3f)", cc$cor[1], p_perm[1])
r2_lab <- sprintf("CV2  (r = %.3f, p = %.3f)", cc$cor[2], p_perm[2])

p_biplot <- ggplot(ld, aes(x = CV1, y = CV2, colour = set, label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_segment(aes(x = 0, y = 0, xend = CV1, yend = CV2),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               linewidth = 0.6, alpha = 0.7) +
  ggrepel::geom_text_repel(size = 3, fontface = "bold", max.overlaps = 20,
                            segment.size = 0.25, box.padding = 0.4) +
  scale_colour_manual(values = c("Model (X)" = "#2B5C8A", "Personality (Y)" = "#CC79A7"),
                      name = NULL) +
  labs(
    title = "CCA Biplot — Study 1 (N = 61)",
    subtitle = "Variable loadings on first two canonical variates",
    x = r1_lab,
    y = r2_lab,
    caption = "Arrows = correlation of original variable with canonical score.\nModel set (blue): α, λ, γ, w, MCR, ΔELPD, Trait Segregation.\nPersonality set (pink): DS, Proto, SCC, SI, RSE, NTB, NFC, SING."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# ggrepel check
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("ggrepel not installed; using geom_text instead")
  p_biplot <- p_biplot +
    geom_text(size = 2.8, vjust = -0.5)
}

ggsave(here("Scripts", "cca_s1_biplot.png"),
       p_biplot, width = 8, height = 7, dpi = 200)

message("\nSaved: Scripts/cca_s1_biplot.png")
message("\nNote: N = 61 with 6 + 9 = 15 variables. CCA is underpowered for later canonical variates.")
message("First 1-2 CVs are interpretable; treat all as exploratory.")

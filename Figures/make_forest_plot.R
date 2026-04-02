# Forest Plot: Sym+λ Parameters Across All Conditions
# Generates fig_forest_plot_pooled.tiff
# Conditions: Study 1 (Minimal Group), Study 2 (Negation/Low-Status/High-Status),
#             Study 3 (Minority/Majority), + Pooled Global diamond

library(tidyverse)
library(here)
library(ggplot2)

# ── Helper: parse subject index from variable name ────────────────────────────
parse_ind_params <- function(df, study_label) {
  df %>%
    filter(str_detect(variable, "^(m|bias|lambda)\\[")) %>%
    mutate(
      param    = str_extract(variable, "^[^\\[]+"),
      subj_idx = as.integer(str_extract(variable, "[0-9]+"))
    ) %>%
    select(param, subj_idx, median) %>%
    mutate(study = study_label)
}

# ── 1. Load individual-level Stan posteriors (already back-transformed) ───────
p1 <- parse_ind_params(read.csv(here("Results", "params_ind_s1_sym_lambda.csv")), "S1")
p2 <- parse_ind_params(read.csv(here("Results", "params_ind_s2_sym_lambda.csv")), "S2")
p3 <- parse_ind_params(read.csv(here("Results", "params_ind_s3_sym_lambda.csv")), "S3")

# ── 2. Build subject-to-condition maps (matching original model_comparison ordering) ─
# Study 1: all Minimal Group (no condition split)
s1_full  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))  %>% filter(!is.na(ingChoiceN))
s1_train <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) %>% filter(!is.na(selfResp))
s1_ids   <- sort(intersect(unique(s1_full$subID), unique(s1_train$subID)))
s1_cond  <- data.frame(
  subj_idx  = seq_along(s1_ids),
  condition = "Minimal Group\n(Study 1)",
  study     = "S1"
)

# Study 2: condition from outgroup column (UCLA=High-Status, Not UCR=Negation, CSU LA=Low-Status)
s2_full  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv"))       %>% filter(!is.na(ingChoiceN))
s2_train <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))
s2_ids   <- sort(intersect(unique(s2_full$subID), unique(s2_train$subID)))
s2_cond_map <- s2_full %>%
  filter(subID %in% s2_ids) %>%
  group_by(subID) %>%
  summarise(outgroup = first(outgroup), .groups = "drop") %>%
  mutate(condition = case_when(
    outgroup == "UCLA"    ~ "High-Status\n(Study 2)",
    outgroup == "Not UCR" ~ "Negation\n(Study 2)",
    outgroup == "CSU LA"  ~ "Low-Status\n(Study 2)"
  ))
s2_cond <- data.frame(subj_idx = seq_along(s2_ids), subID = s2_ids) %>%
  left_join(s2_cond_map %>% select(subID, condition), by = "subID") %>%
  mutate(study = "S2") %>%
  select(subj_idx, condition, study)

# Study 3: condition column (Minority/Majority)
s3_full  <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv"))  %>% filter(!is.na(ingChoiceN))
s3_train <- read.csv(here("Study 3/Cleaning/output/fullTrain_fixed.csv")) %>% filter(!is.na(selfResp))
s3_ids   <- sort(intersect(unique(s3_full$subID), unique(s3_train$subID)))
s3_cond_map <- s3_full %>%
  filter(subID %in% s3_ids) %>%
  group_by(subID) %>%
  summarise(condition_raw = first(condition), .groups = "drop")
s3_cond <- data.frame(subj_idx = seq_along(s3_ids), subID = s3_ids) %>%
  left_join(s3_cond_map %>% select(subID, condition_raw), by = "subID") %>%
  mutate(
    condition = paste0(condition_raw, "\n(Study 3)"),
    study     = "S3"
  ) %>%
  select(subj_idx, condition, study)

# ── 3. Join params with conditions ───────────────────────────────────────────
cond_all  <- bind_rows(s1_cond, s2_cond, s3_cond)
params_all <- bind_rows(p1, p2, p3) %>%
  left_join(cond_all, by = c("subj_idx", "study"))

# ── 4. Condition-level summary (mean ± 90% CI via 1.645*SE; n available) ─────
cond_summary <- params_all %>%
  group_by(condition, param) %>%
  summarise(
    n    = n(),
    mean = mean(median),
    se   = sd(median) / sqrt(n()),
    lo   = mean - 1.96 * se,
    hi   = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(is_pooled = FALSE, shape = "circle")

# ── 5. Pooled global from summary_pooled_sym_lambda.csv ──────────────────────
# global_mu_pr[1..3] in probit space; back-transform with Phi * scale
# Parameter order: [m, bias, lambda]; scales [10, 1, 5]
scales <- c(m = 10, bias = 1, lambda = 5)

pooled_sum <- read.csv(here("Results", "summary_pooled_sym_lambda.csv"))
global_rows <- pooled_sum %>%
  filter(str_detect(variable, "^global_mu_pr")) %>%
  distinct(variable, .keep_all = TRUE) %>%     # deduplicate
  arrange(variable)                              # [1],[2],[3],[4]

global_df <- data.frame(
  condition = "Pooled Global",
  param     = c("m", "bias", "lambda"),
  n         = NA_integer_,
  mean      = pnorm(global_rows$median) * scales,
  se        = NA_real_,
  lo        = pnorm(global_rows$q5)     * scales,
  hi        = pnorm(global_rows$q95)    * scales,
  is_pooled = TRUE,
  shape     = "diamond"
)

plot_df <- bind_rows(cond_summary, global_df)

# ── 6. Factor levels (top = Study 1, bottom = Pooled) ────────────────────────
cond_levels <- c(
  "Minimal Group\n(Study 1)",
  "Negation\n(Study 2)",
  "Low-Status\n(Study 2)",
  "High-Status\n(Study 2)",
  "Minority\n(Study 3)",
  "Majority\n(Study 3)",
  "Pooled Global"
)
param_labels <- c(
  m      = "\u03b1 (Projection Rate)",
  bias   = "\u03b3 (Ingroup Bias)",
  lambda = "\u03bb (Similarity Weight)"
)

plot_df <- plot_df %>%
  mutate(
    condition   = factor(condition, levels = rev(cond_levels)),  # rev so Study1 at top
    param_label = factor(param, levels = names(param_labels), labels = param_labels)
  )

# ── 7. Study band colors ──────────────────────────────────────────────────────
study_colors <- c(
  "Minimal Group\n(Study 1)"  = "#4E79A7",
  "Negation\n(Study 2)"       = "#E15759",
  "Low-Status\n(Study 2)"     = "#F28E2B",
  "High-Status\n(Study 2)"    = "#76B7B2",
  "Minority\n(Study 3)"       = "#59A14F",
  "Majority\n(Study 3)"       = "#B07AA1",
  "Pooled Global"             = "firebrick"
)

# ── 8. Plot ───────────────────────────────────────────────────────────────────
p <- ggplot(
  plot_df %>% filter(!is_pooled),
  aes(x = mean, y = condition, color = condition)
) +
  # Reference lines at pooled global mean
  geom_vline(
    data = global_df %>%
      mutate(param_label = factor(param, levels = names(param_labels), labels = param_labels)),
    aes(xintercept = mean),
    linetype = "dashed", color = "grey50", linewidth = 0.4, inherit.aes = FALSE
  ) +
  # Study condition CIs and points
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.7) +
  geom_point(aes(size = n), shape = 19) +
  # Pooled global diamond + CI
  geom_point(
    data = plot_df %>% filter(is_pooled),
    aes(x = mean, y = condition),
    shape = 18, size = 5, color = "firebrick", inherit.aes = FALSE
  ) +
  geom_errorbarh(
    data = plot_df %>% filter(is_pooled),
    aes(xmin = lo, xmax = hi, y = condition),
    height = 0.2, color = "firebrick", linewidth = 0.9, inherit.aes = FALSE
  ) +
  facet_wrap(~ param_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = study_colors, guide = "none") +
  scale_size_continuous(range = c(2.5, 5), guide = "none") +
  labs(
    x = "Posterior Median (individual-level)",
    y = NULL,
    caption = "Points = condition means; bars = 95% CI; diamond = pooled global (90% CI)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background  = element_rect(fill = "grey93", color = NA),
    strip.text        = element_text(face = "bold", size = 9),
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y       = element_text(size = 8.5),
    plot.caption      = element_text(size = 7, color = "grey50")
  )

# ── 9. Save ───────────────────────────────────────────────────────────────────
ggsave(
  here("Figures", "fig_forest_plot_pooled.tiff"),
  plot = p, width = 11, height = 4.5, dpi = 300,
  compression = "lzw"
)
message("Saved: Figures/fig_forest_plot_pooled.tiff")

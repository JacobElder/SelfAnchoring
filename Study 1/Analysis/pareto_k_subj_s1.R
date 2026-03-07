# Per-subject Pareto k diagnostics for Study 1
# Run after run_model_comparison_s1.R completes all 4 models.
# Produces: Results/pareto_k_subj_s1_<model>.csv

library(loo)
library(tidyverse)
library(here)

# --- Reconstruct subject info from data ---
fulldf  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))  %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) %>% filter(!is.na(selfResp))
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
fulldf  <- fulldf  %>% filter(subID %in% common_ids)
uIds    <- sort(common_ids)
maxTrials <- max(fulldf$trialTotal)

# Per-subject nTrials (must match what was passed to Stan)
nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf, subID == id)))

# --- Helper ---
compute_subj_pareto_k <- function(l, study_label, model_name, uIds, nTrials_vec, maxTrials) {
  pk <- l$diagnostics$pareto_k
  rows <- lapply(seq_along(uIds), function(i) {
    idx <- ((i - 1) * maxTrials + 1):((i - 1) * maxTrials + nTrials_vec[i])
    k_i <- pk[idx]
    data.frame(
      study        = study_label,
      model        = model_name,
      subID        = uIds[i],
      subj_idx     = i,
      n_trials     = nTrials_vec[i],
      k_mean       = round(mean(k_i), 4),
      k_max        = round(max(k_i), 4),
      n_good       = sum(k_i < 0.5),
      n_ok         = sum(k_i >= 0.5 & k_i < 0.7),
      n_bad        = sum(k_i >= 0.7 & k_i < 1.0),
      n_verybad    = sum(k_i >= 1.0),
      pct_reliable = round(100 * mean(k_i < 0.7), 1),
      concern      = ifelse(any(k_i >= 1.0), "HIGH",
                    ifelse(any(k_i >= 0.7), "MODERATE",
                    ifelse(any(k_i >= 0.5), "LOW", "NONE")))
    )
  })
  do.call(rbind, rows)
}

# --- Process each model ---
models <- c("bias", "symmetric", "sym_lambda", "asym_lambda")
all_rows <- list()

for (m in models) {
  loo_path <- here("Fits", paste0("loo_s1_", m, ".rds"))
  if (!file.exists(loo_path)) {
    message("Skipping ", m, " — LOO file not found yet.")
    next
  }
  l <- readRDS(loo_path)
  df <- compute_subj_pareto_k(l, "S1", m, uIds, nTrials_vec, maxTrials)
  write.csv(df, here("Results", paste0("pareto_k_subj_s1_", m, ".csv")), row.names = FALSE)
  message("Saved pareto_k_subj_s1_", m, ".csv")
  all_rows[[m]] <- df
}

# Combined across models
if (length(all_rows) > 0) {
  combined <- do.call(rbind, all_rows)
  write.csv(combined, here("Results", "pareto_k_subj_s1_all_models.csv"), row.names = FALSE)

  # Summary: subjects with HIGH concern in any model
  high_concern <- combined %>%
    filter(concern == "HIGH") %>%
    group_by(subID) %>%
    summarise(
      models_flagged = paste(model, collapse = ", "),
      max_k_any_model = round(max(k_max), 3),
      max_verybad = max(n_verybad),
      .groups = "drop"
    ) %>%
    arrange(desc(max_k_any_model))

  cat("\n=== Subjects with HIGH Pareto k concern (k > 1) ===\n")
  print(high_concern)
}

message("Study 1 per-subject Pareto k complete.")

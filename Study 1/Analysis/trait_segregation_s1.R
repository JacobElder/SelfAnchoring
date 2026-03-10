# Trait Segregation Correlations — Study 1
# Correlates groupHomoph (network nominal homophily) with:
#   - Computational parameters (m, bias, lambda, w)
#   - Model-derived MCR (subject_mcr)
#   - Per-subject ΔELPD
#   - Individual difference scales
# FDR correction applied across the full family of tests.

library(tidyverse)
library(here)

# ── 1. Load subject-level data ────────────────────────────────────────────────
fulldf <- read.csv(here("Study 1/Cleaning/output/fullTest.csv")) |>
  filter(!is.na(ingChoiceN))

uIds <- sort(unique(fulldf$subID))

scale_vars <- c("DS","Proto","SCC","SI","RSE","NTB","NFC","SING.Ind","SING.Inter")
avail_scales <- intersect(scale_vars, names(fulldf))

# One row per subject
id_df <- fulldf[!duplicated(fulldf$subID), ]
id_df <- id_df[order(match(id_df$subID, uIds)), ]
id_df$subj_idx <- seq_len(nrow(id_df))
id_df <- id_df[, c("subj_idx","subID", "groupHomoph", avail_scales)]

# ── 2. Load computational parameters and ΔELPD ────────────────────────────────
params_file <- here("Results","params_ind_s1_sym_lambda.csv")
if (!file.exists(params_file)) params_file <- here("Results","params_ind_s1_asym_lambda.csv")

params_sl <- read.csv(params_file)
params_sl <- params_sl[grepl("^(m|bias|lambda|w)\\[", params_sl$variable), ]
params_sl$subj_idx <- as.integer(regmatches(params_sl$variable,
                                             regexpr("[0-9]+", params_sl$variable)))
params_sl$param    <- sub("\\[.*", "", params_sl$variable)
params_wide <- reshape(params_sl[, c("subj_idx","param","median")],
  idvar="subj_idx", timevar="param", direction="wide")
names(params_wide) <- sub("median\\.", "", names(params_wide))

# MCR from asym_lambda
params_al <- read.csv(here("Results","params_ind_s1_asym_lambda.csv"))
mcr <- params_al[grepl("^subject_mcr\\[", params_al$variable), ]
mcr$subj_idx <- as.integer(regmatches(mcr$variable, regexpr("[0-9]+", mcr$variable)))
mcr <- mcr[, c("subj_idx","median")]; names(mcr)[2] <- "subject_mcr"

# ΔELPD
loo_path_sym  <- here("Fits","loo_s1_sym_lambda.rds")
loo_path_asym <- here("Fits","loo_s1_asym_lambda.rds")

delta_elpd_df <- NULL
if (file.exists(loo_path_sym) && file.exists(loo_path_asym)) {
  traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) |>
    filter(!is.na(selfResp))
  common_ids <- sort(intersect(uIds, unique(traindf$subID)))
  maxTrials <- max(fulldf$trialTotal, na.rm=TRUE)
  nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf, subID == id)))
  loo_sym  <- readRDS(loo_path_sym)
  loo_asym <- readRDS(loo_path_asym)
  elpd_sym  <- loo_sym$pointwise[, "elpd_loo"]
  elpd_asym <- loo_asym$pointwise[, "elpd_loo"]
  delta_elpd_vec <- sapply(seq_along(uIds), function(i) {
    idx <- ((i-1)*maxTrials + 1) : ((i-1)*maxTrials + nTrials_vec[i])
    sum(elpd_asym[idx]) - sum(elpd_sym[idx])
  })
  delta_elpd_df <- data.frame(subj_idx = seq_along(uIds), delta_elpd = delta_elpd_vec)
}

# ── 3. Merge ──────────────────────────────────────────────────────────────────
df <- merge(params_wide, mcr, by="subj_idx")
if (!is.null(delta_elpd_df)) df <- merge(df, delta_elpd_df, by="subj_idx")
df <- merge(df, id_df, by="subj_idx")
cat(sprintf("Merged N: %d\n", nrow(df)))

# ── 4. Correlate groupHomoph with params, MCR, ΔELPD, and scales ──────────────
# FDR correction over theoretically meaningful set:
#   alpha/m, lambda, subject_mcr, delta_elpd, and individual-difference scales.
# bias and w are reported descriptively but excluded from FDR family.
fdr_preds   <- intersect(c("m","lambda","subject_mcr","delta_elpd", avail_scales), names(df))
descr_preds <- intersect(c("bias","w"), names(df))

run_cors_1 <- function(pred_set, df) {
  rows <- list()
  for (p in pred_set) {
    ct <- cor.test(df$groupHomoph, df[[p]], use="complete.obs")
    rows[[length(rows)+1]] <- data.frame(
      predictor = p,
      r         = round(ct$estimate, 3),
      p_raw     = round(ct$p.value, 4),
      n         = sum(complete.cases(df[, c("groupHomoph", p)]))
    )
  }
  do.call(rbind, rows)
}

res_fdr   <- run_cors_1(fdr_preds,   df)
res_descr <- run_cors_1(descr_preds, df)
res_fdr$p_fdr   <- round(p.adjust(res_fdr$p_raw, method="BH"), 4)
res_descr$p_fdr <- NA_real_
res <- rbind(res_fdr, res_descr)
res <- res[order(res$p_raw), ]

cat("\n=== Trait Segregation correlations (FDR-corrected predictors) ===\n")
print(res_fdr[order(res_fdr$p_raw), ], row.names=FALSE)
cat("\n=== Descriptive: bias, w (not FDR-corrected) ===\n")
print(res_descr[order(res_descr$p_raw), ], row.names=FALSE)
cat(sprintf("\nDescriptives: M = %.3f, SD = %.3f, range [%.3f, %.3f]\n",
  mean(df$groupHomoph, na.rm=TRUE), sd(df$groupHomoph, na.rm=TRUE),
  min(df$groupHomoph, na.rm=TRUE), max(df$groupHomoph, na.rm=TRUE)))

write.csv(res, here("Results","trait_segregation_s1_correlations.csv"), row.names=FALSE)
message("Saved: Results/trait_segregation_s1_correlations.csv")

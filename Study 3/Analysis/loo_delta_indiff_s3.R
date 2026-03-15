# Per-Subject LOO ELPD Difference (asym_lambda - sym_lambda) — Study 3
# Run AFTER run_model_comparison_s3.R completes (requires Fits/loo_s3_*.rds)
# Parallel structure to loo_delta_indiff_s1.R / loo_delta_indiff_s2.R

library(tidyverse)
library(here)

# ── 1. Reconstruct subject trial counts ──────────────────────────────────────
fulldf  <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv"))  |> filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 3/Cleaning/output/fullTrain_fixed.csv")) |> filter(!is.na(selfResp))

common_ids <- sort(intersect(unique(fulldf$subID), unique(traindf$subID)))
fulldf <- filter(fulldf, subID %in% common_ids)

uIds      <- sort(common_ids)
nSubjects <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2)
nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf, subID == id)))

cat(sprintf("Subjects: %d, maxTrials: %d\n", nSubjects, maxTrials))

# ── 2. Load LOO objects ───────────────────────────────────────────────────────
loo_sym_path  <- here("Fits","loo_s3_sym_lambda.rds")
loo_asym_path <- here("Fits","loo_s3_asym_lambda.rds")

if (!file.exists(loo_sym_path) || !file.exists(loo_asym_path)) {
  stop("LOO files not found. Run run_model_comparison_s3.R first.")
}

loo_sym  <- readRDS(loo_sym_path)
loo_asym <- readRDS(loo_asym_path)
elpd_sym  <- loo_sym$pointwise[, "elpd_loo"]
elpd_asym <- loo_asym$pointwise[, "elpd_loo"]

# ── 3. Sum per subject (valid trials only) ────────────────────────────────────
subj_elpd <- data.frame(
  subj_idx  = seq_along(uIds),
  subID     = uIds,
  n_trials  = nTrials_vec,
  elpd_sym  = NA_real_,
  elpd_asym = NA_real_
)
for (i in seq_along(uIds)) {
  idx <- ((i-1)*maxTrials + 1) : ((i-1)*maxTrials + nTrials_vec[i])
  subj_elpd$elpd_sym[i]  <- sum(elpd_sym[idx])
  subj_elpd$elpd_asym[i] <- sum(elpd_asym[idx])
}
subj_elpd$delta_elpd <- subj_elpd$elpd_asym - subj_elpd$elpd_sym

cat("ΔELPD summary:\n"); print(summary(subj_elpd$delta_elpd))
cat(sprintf("Total ΔELPD = %.3f\n", sum(subj_elpd$delta_elpd)))

# ── 4. Merge with structural params, MCR, individual differences ──────────────
params_file <- here("Results","params_ind_s3_sym_lambda.csv")
if (!file.exists(params_file)) params_file <- here("Results","params_ind_s3_asym_lambda.csv")

params_sl <- read.csv(params_file)
params_sl <- params_sl[grepl("^(m|bias|lambda|w)\\[", params_sl$variable), ]
params_sl$subj_idx <- as.integer(regmatches(params_sl$variable, regexpr("[0-9]+", params_sl$variable)))
params_sl$param    <- sub("\\[.*", "", params_sl$variable)
params_wide <- reshape(params_sl[, c("subj_idx","param","median")],
  idvar="subj_idx", timevar="param", direction="wide")
names(params_wide) <- sub("median\\.", "", names(params_wide))

params_al <- read.csv(here("Results","params_ind_s3_asym_lambda.csv"))
mcr <- params_al[grepl("^subject_mcr\\[", params_al$variable), ]
mcr$subj_idx <- as.integer(regmatches(mcr$variable, regexpr("[0-9]+", mcr$variable)))
mcr <- mcr[, c("subj_idx","median")]; names(mcr)[2] <- "subject_mcr"

scale_vars <- c("DS","Proto","SCC","SI","RSE","NTB","NFC","SING.Ind","SING.Inter")
avail_scales <- intersect(scale_vars, names(fulldf))
id_df <- fulldf[!duplicated(fulldf$subID), ]
id_df <- id_df[order(match(id_df$subID, uIds)), ]
id_df$subj_idx <- seq_len(nrow(id_df))
id_df <- id_df[, c("subj_idx","subID","condition", avail_scales)]

df <- merge(params_wide, mcr, by="subj_idx")
df <- merge(df, subj_elpd[, c("subj_idx","delta_elpd")], by="subj_idx")
df <- merge(df, id_df, by="subj_idx")
cat(sprintf("Merged N: %d\n", nrow(df)))

# ── 5. Correlations ───────────────────────────────────────────────────────────
# Two separate FDR families:
#   Family A (mechanism-personality): m, bias/γ, lambda, subject_mcr × scales
#   Family B (ΔELPD-personality): delta_elpd × scales — separate question (model preference)
# w descriptive only.
mech_params  <- intersect(c("m","m_in","m_out","bias","lambda","subject_mcr"), names(df))
elpd_params  <- intersect(c("delta_elpd"), names(df))
descr_params <- intersect(c("w"), names(df))

run_cors <- function(param_set, df, scale_cols) {
  rows <- list()
  for (p in param_set) for (s in scale_cols) {
    ct <- cor.test(df[[p]], df[[s]], use="complete.obs")
    rows[[length(rows)+1]] <- data.frame(param=p, scale=s,
      r=round(ct$estimate,3), p_raw=round(ct$p.value,4),
      n=sum(complete.cases(df[,c(p,s)])))
  }
  do.call(rbind, rows)
}

res_mech  <- run_cors(mech_params,  df, avail_scales)
res_elpd  <- run_cors(elpd_params,  df, avail_scales)
res_descr <- run_cors(descr_params, df, avail_scales)
res_mech$p_fdr  <- round(p.adjust(res_mech$p_raw,  method="BH"), 4)
res_elpd$p_fdr  <- round(p.adjust(res_elpd$p_raw,  method="BH"), 4)
res_descr$p_fdr <- NA_real_
res <- rbind(res_mech, res_elpd, res_descr)
res <- res[order(res$p_raw), ]

cat("\n=== ΔELPD correlations ===\n")
print(res_elpd[order(res_elpd$p_raw), ], row.names=FALSE)
cat("\n=== Nominally significant (p_raw < .05) — mechanism params ===\n")
print(res_mech[res_mech$p_raw < .05, ], row.names=FALSE)
cat("\n=== FDR-significant mechanism params (p_fdr < .05) ===\n")
print(res_mech[!is.na(res_mech$p_fdr) & res_mech$p_fdr < .05, ], row.names=FALSE)
cat("\n=== FDR-significant ΔELPD (p_fdr < .05) ===\n")
print(res_elpd[!is.na(res_elpd$p_fdr) & res_elpd$p_fdr < .05, ], row.names=FALSE)
cat("\n=== Descriptive: w (not FDR-corrected) ===\n")
print(res_descr[order(res_descr$p_raw), ], row.names=FALSE)

write.csv(df,  here("Results","ind_diffs_s3_full.csv"),    row.names=FALSE)
write.csv(res, here("Results","correlations_s3_full.csv"), row.names=FALSE)
message("Saved: Results/ind_diffs_s3_full.csv, Results/correlations_s3_full.csv")

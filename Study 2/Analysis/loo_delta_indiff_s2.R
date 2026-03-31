# Per-Subject LOO ELPD Difference (asym_lambda - sym_lambda) — Study 2
# Run AFTER run_model_comparison_s2.R completes (requires Fits/loo_s2_*.rds)
# subject_mcr computed analytically from asym_lambda posterior medians (NoW).

library(tidyverse)
library(here)
library(igraph)

# ── 1. Load data ──────────────────────────────────────────────────────────────
fulldf  <- read.csv(here("Study 2/Cleaning/output/fullTest.csv"))       |> filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 2/Cleaning/output/fullTrain_fixed.csv")) |> filter(!is.na(selfResp))

common_ids <- sort(intersect(unique(fulldf$subID), unique(traindf$subID)))
fulldf  <- filter(fulldf,  subID %in% common_ids)
traindf <- filter(traindf, subID %in% common_ids)

posDf    <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
simMat   <- similarity(graph_from_adjacency_matrix(as.matrix(posDf), mode = "max"),
                       method = "dice")

uIds        <- sort(common_ids)
nSubjects   <- length(uIds)
maxTrials   <- max(fulldf$trialTotalT2)
nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf, subID == id)))

cat(sprintf("Subjects: %d, maxTrials: %d\n", nSubjects, maxTrials))

# ── 2. Load LOO objects ───────────────────────────────────────────────────────
loo_sym  <- readRDS(here("Fits", "loo_s2_sym_lambda.rds"))
loo_asym <- readRDS(here("Fits", "loo_s2_asym_lambda.rds"))
elpd_sym  <- loo_sym$pointwise[, "elpd_loo"]
elpd_asym <- loo_asym$pointwise[, "elpd_loo"]

# ── 3. Sum per subject (valid trials only) ────────────────────────────────────
subj_elpd <- data.frame(
  subj_idx  = seq_along(uIds), subID = uIds, n_trials = nTrials_vec,
  elpd_sym  = NA_real_, elpd_asym = NA_real_
)
for (i in seq_along(uIds)) {
  idx <- ((i-1)*maxTrials + 1) : ((i-1)*maxTrials + nTrials_vec[i])
  subj_elpd$elpd_sym[i]  <- sum(elpd_sym[idx])
  subj_elpd$elpd_asym[i] <- sum(elpd_asym[idx])
}
subj_elpd$delta_elpd <- subj_elpd$elpd_asym - subj_elpd$elpd_sym

cat("ΔELPD summary:\n"); print(summary(subj_elpd$delta_elpd))
cat(sprintf("Subjects favoring asym (delta > 0): %d / %d\n",
            sum(subj_elpd$delta_elpd > 0), nSubjects))
cat(sprintf("Total ΔELPD = %.3f\n", sum(subj_elpd$delta_elpd)))

# ── 4. Structural params from sym_lambda (winner) ────────────────────────────
params_sl <- read.csv(here("Results", "params_ind_s2_sym_lambda.csv"))
params_sl <- params_sl[grepl("^(m|bias|lambda)\\[", params_sl$variable), ]
params_sl$subj_idx <- as.integer(regmatches(params_sl$variable, regexpr("[0-9]+", params_sl$variable)))
params_sl$param    <- sub("\\[.*", "", params_sl$variable)
params_wide <- reshape(params_sl[, c("subj_idx","param","median")],
  idvar = "subj_idx", timevar = "param", direction = "wide")
names(params_wide) <- sub("median\\.", "", names(params_wide))

# ── 5. subject_mcr — analytical from sym_lambda posterior medians (winning model)
# Consistent with S1/S3. MCR derived from sym model params avoids relying on
# the non-winning asym architecture.
get_pv <- function(df, prefix) {
  rows <- df[grepl(paste0("^", prefix, "\\["), df$variable), ]
  rows$idx <- as.integer(regmatches(rows$variable, regexpr("[0-9]+", rows$variable)))
  rows$median[order(rows$idx)]
}
params_sl_raw <- read.csv(here("Results", "params_ind_s2_sym_lambda.csv"))
m_v   <- get_pv(params_sl_raw, "m")
lam_v <- get_pv(params_sl_raw, "lambda")

mcr_vals <- sapply(seq_along(uIds), function(s) {
  id      <- uIds[s]
  s_df    <- filter(fulldf,  subID == id)
  s_train <- filter(traindf, subID == id)
  GP       <- plogis(m_v[s] * (s_train$selfResp - 4))
  PS       <- simMat[s_df$Idx, s_train$Idx]^lam_v[s]
  simW_in  <- as.numeric(PS %*% GP) + 1e-9
  simW_out <- as.numeric(PS %*% (1 - GP)) + 1e-9
  mean(simW_in / simW_out)
})
mcr <- data.frame(subj_idx = seq_along(uIds), subject_mcr = mcr_vals)
cat("\nsubject_mcr summary:\n"); print(summary(mcr$subject_mcr))

# ── 6. Merge ──────────────────────────────────────────────────────────────────
scale_vars   <- c("DS","Proto","SCC","SI","RSE","NTB","NFC","SING.Ind","SING.Inter")
avail_scales <- intersect(scale_vars, names(fulldf))
id_df <- fulldf[!duplicated(fulldf$subID), ]
id_df <- id_df[order(match(id_df$subID, uIds)), ]
id_df$subj_idx <- seq_len(nrow(id_df))
id_df <- id_df[, c("subj_idx","subID","outgroup", avail_scales)]

df <- merge(params_wide,                               mcr,         by = "subj_idx")
df <- merge(df, subj_elpd[, c("subj_idx","delta_elpd")],            by = "subj_idx")
df <- merge(df,                                        id_df,        by = "subj_idx")
cat(sprintf("Merged N: %d\n", nrow(df)))

# ── 7. Correlations ───────────────────────────────────────────────────────────
# Family A: m, bias/γ, lambda, subject_mcr × 9 scales (k=36)
# Family B: delta_elpd × 9 scales (k=9, separate)
mech_params <- intersect(c("m","bias","lambda","subject_mcr"), names(df))
elpd_params <- "delta_elpd"

run_cors <- function(param_set, df, scale_cols) {
  rows <- list()
  for (p in param_set) for (s in scale_cols) {
    ct <- cor.test(df[[p]], df[[s]], use = "complete.obs")
    rows[[length(rows)+1]] <- data.frame(param=p, scale=s,
      r=round(ct$estimate,3), p_raw=round(ct$p.value,4),
      n=sum(complete.cases(df[,c(p,s)])))
  }
  do.call(rbind, rows)
}

res_mech <- run_cors(mech_params, df, avail_scales)
res_elpd <- run_cors(elpd_params, df, avail_scales)
res_mech$p_fdr <- round(p.adjust(res_mech$p_raw, method = "BH"), 4)
res_elpd$p_fdr <- round(p.adjust(res_elpd$p_raw, method = "BH"), 4)
res <- rbind(res_mech, res_elpd)
res <- res[order(res$p_raw), ]

cat("\n=== ΔELPD correlations ===\n")
print(res_elpd[order(res_elpd$p_raw), ], row.names = FALSE)
cat("\n=== Nominally significant (p_raw < .05) — mechanism params ===\n")
print(res_mech[res_mech$p_raw < .05, ], row.names = FALSE)
cat("\n=== FDR-significant mechanism params (p_fdr < .05) ===\n")
print(res_mech[!is.na(res_mech$p_fdr) & res_mech$p_fdr < .05, ], row.names = FALSE)
cat("\n=== FDR-significant ΔELPD (p_fdr < .05) ===\n")
print(res_elpd[!is.na(res_elpd$p_fdr) & res_elpd$p_fdr < .05, ], row.names = FALSE)

write.csv(df,  here("Results","ind_diffs_s2_full.csv"),    row.names = FALSE)
write.csv(res, here("Results","correlations_s2_full.csv"), row.names = FALSE)
message("Saved: Results/ind_diffs_s2_full.csv, Results/correlations_s2_full.csv")

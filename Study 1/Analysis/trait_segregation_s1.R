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
params_sl <- params_sl[grepl("^(m|bias|lambda)\\[", params_sl$variable), ]
params_sl$subj_idx <- as.integer(regmatches(params_sl$variable,
                                             regexpr("[0-9]+", params_sl$variable)))
params_sl$param    <- sub("\\[.*", "", params_sl$variable)
params_wide <- reshape(params_sl[, c("subj_idx","param","median")],
  idvar="subj_idx", timevar="param", direction="wide")
names(params_wide) <- sub("median\\.", "", names(params_wide))

# MCR computed analytically from asym_lambda posterior medians
# (subject_mcr not in Stan output for NoW models; computed here from params)
library(igraph)
traindf_ts <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) |>
  filter(!is.na(selfResp), subID %in% uIds)
posDf_ts   <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
simMat_ts  <- similarity(graph_from_adjacency_matrix(as.matrix(posDf_ts), mode="max"),
                         method="dice")
params_al  <- read.csv(here("Results","params_ind_s1_asym_lambda.csv"))
get_pv <- function(df, prefix) {
  rows <- df[grepl(paste0("^", prefix, "\\["), df$variable), ]
  rows$idx <- as.integer(regmatches(rows$variable, regexpr("[0-9]+", rows$variable)))
  rows$median[order(rows$idx)]
}
m_in_v  <- get_pv(params_al, "m_in")
m_out_v <- get_pv(params_al, "m_out")
lam_v   <- get_pv(params_al, "lambda")
mcr_vals <- sapply(seq_along(uIds), function(s) {
  id <- uIds[s]
  s_df    <- filter(fulldf,      subID == id)
  s_train <- filter(traindf_ts,  subID == id)
  GPin  <- plogis( m_in_v[s]  * (s_train$selfResp - 4))
  GPout <- plogis(-m_out_v[s] * (s_train$selfResp - 4))
  PS    <- simMat_ts[s_df$Idx, s_train$Idx]^lam_v[s]
  simW_in  <- as.numeric(PS %*% GPin)  + 1e-9
  simW_out <- as.numeric(PS %*% GPout) + 1e-9
  mean(simW_in / simW_out)
})
mcr <- data.frame(subj_idx = seq_along(uIds), subject_mcr = mcr_vals)

# ΔELPD
loo_path_sym  <- here("Fits","loo_s1_sym_lambda.rds")
loo_path_asym <- here("Fits","loo_s1_asym_lambda.rds")

delta_elpd_df <- NULL
if (file.exists(loo_path_sym) && file.exists(loo_path_asym)) {
  traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) |>
    filter(!is.na(selfResp))
  common_ids <- sort(intersect(uIds, unique(traindf$subID)))
  maxTrials <- max(fulldf$trialTotalT2, na.rm=TRUE)
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
# Two separate FDR families:
#   Family 1: computational params (m, bias/γ, lambda, subject_mcr, delta_elpd)
#   Family 2: individual-difference scales
# w is descriptive only.
param_preds <- intersect(c("m","bias","lambda","subject_mcr","delta_elpd"), names(df))
descr_preds <- character(0)  # w removed (NoW models)

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

res_params  <- run_cors_1(param_preds, df)
res_scales  <- run_cors_1(avail_scales, df)
res_params$p_fdr <- round(p.adjust(res_params$p_raw, method="BH"), 4)
res_scales$p_fdr <- round(p.adjust(res_scales$p_raw, method="BH"), 4)
res <- rbind(res_params, res_scales)
res <- res[order(res$p_raw), ]

cat("\n=== Trait Segregation ~ Computational Params (FDR-corrected) ===\n")
print(res_params[order(res_params$p_raw), ], row.names=FALSE)
cat("\n=== Trait Segregation ~ Personality Scales (FDR-corrected) ===\n")
print(res_scales[order(res_scales$p_raw), ], row.names=FALSE)
cat("\n=== FDR-significant params ===\n")
print(res_params[!is.na(res_params$p_fdr) & res_params$p_fdr < .05, ], row.names=FALSE)
cat("\n=== FDR-significant scales ===\n")
print(res_scales[!is.na(res_scales$p_fdr) & res_scales$p_fdr < .05, ], row.names=FALSE)
cat(sprintf("\nDescriptives: M = %.3f, SD = %.3f, range [%.3f, %.3f]\n",
  mean(df$groupHomoph, na.rm=TRUE), sd(df$groupHomoph, na.rm=TRUE),
  min(df$groupHomoph, na.rm=TRUE), max(df$groupHomoph, na.rm=TRUE)))

write.csv(res, here("Results","trait_segregation_s1_correlations.csv"), row.names=FALSE)
message("Saved: Results/trait_segregation_s1_correlations.csv")

# Per-Subject LOO ELPD Difference (asym_lambda - sym_lambda) — Study 1
#
# Computes a person-level model fit metric: how much better (or worse) each
# participant's data is explained by the asymmetric architecture relative to
# the symmetric one. Positive ΔELPD_subj → participant is better described by
# distinct ingroup/outgroup slopes; negative → symmetric projection is sufficient.
#
# subject_mcr is computed analytically from asym_lambda posterior medians
# (m_in, m_out, lambda) applied to the actual trial data. This is equivalent
# to the Stan generated quantity and avoids refitting the model.
#
# NOTE on padding: Stan initializes log_lik to 0 for unused trials
# (groupChoice[s,t] = 0). These contribute 0 elpd_loo (approx) and must be
# excluded when summing per-subject. We identify valid trials from the actual data.

library(tidyverse)
library(here)
library(igraph)

# ── 1. Load data (same as run_model_comparison_s1.R) ─────────────────────────
fulldf  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))  |> filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) |> filter(!is.na(selfResp))

common_ids <- sort(intersect(unique(fulldf$subID), unique(traindf$subID)))
fulldf  <- filter(fulldf,  subID %in% common_ids)
traindf <- filter(traindf, subID %in% common_ids)

# Similarity matrix (trait network)
posDf    <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
posMat   <- as.matrix(posDf)
posGraph <- graph_from_adjacency_matrix(posMat, mode = "max")
simMat   <- similarity(posGraph, method = "dice")

uIds      <- sort(common_ids)
nSubjects <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2)
maxTrain  <- 91

nTrials_vec <- sapply(uIds, function(id) nrow(filter(fulldf,  subID == id)))
nTrain_vec  <- sapply(uIds, function(id) nrow(filter(traindf, subID == id)))

cat(sprintf("Subjects: %d, maxTrials: %d, total obs: %d\n",
            nSubjects, maxTrials, nSubjects * maxTrials))

# ── 2. Load LOO objects ───────────────────────────────────────────────────────
loo_sym  <- readRDS(here("Fits", "loo_s1_sym_lambda.rds"))
loo_asym <- readRDS(here("Fits", "loo_s1_asym_lambda.rds"))

elpd_sym  <- loo_sym$pointwise[, "elpd_loo"]
elpd_asym <- loo_asym$pointwise[, "elpd_loo"]

cat(sprintf("LOO dimensions: sym=%d, asym=%d\n", length(elpd_sym), length(elpd_asym)))

# ── 3. Sum per subject (valid trials only) ────────────────────────────────────
subj_elpd <- data.frame(
  subj_idx   = 1:nSubjects,
  subID      = uIds,
  n_trials   = nTrials_vec,
  elpd_sym   = NA_real_,
  elpd_asym  = NA_real_
)

for (i in 1:nSubjects) {
  idx <- ((i - 1) * maxTrials + 1) : ((i - 1) * maxTrials + nTrials_vec[i])
  subj_elpd$elpd_sym[i]  <- sum(elpd_sym[idx])
  subj_elpd$elpd_asym[i] <- sum(elpd_asym[idx])
}

subj_elpd$delta_elpd <- subj_elpd$elpd_asym - subj_elpd$elpd_sym
# Positive delta → asym_lambda fits better; negative → sym_lambda fits better

cat(sprintf("\nΔELPD (asym - sym) per subject summary:\n"))
print(summary(subj_elpd$delta_elpd))
cat(sprintf("Subjects favoring asym_lambda (delta > 0): %d / %d\n",
            sum(subj_elpd$delta_elpd > 0), nSubjects))
cat(sprintf("Total ΔELPD = %.3f (reference ≈ -8.31)\n", sum(subj_elpd$delta_elpd)))

# ── 4. Compute subject_mcr analytically from asym_lambda posterior medians ────
# MCR_i = mean_t( simW_in[t] / simW_out[t] ) using posterior medians of
# m_in[i], m_out[i], lambda[i].  Matches the Stan generated quantity exactly.

params_al_raw <- read.csv(here("Results", "params_ind_s1_asym_lambda.csv"))

get_param_vec <- function(df, prefix) {
  rows <- df[grepl(paste0("^", prefix, "\\["), df$variable), ]
  rows$idx <- as.integer(regmatches(rows$variable, regexpr("[0-9]+", rows$variable)))
  rows$median[order(rows$idx)]
}

m_in_vec    <- get_param_vec(params_al_raw, "m_in")
m_out_vec   <- get_param_vec(params_al_raw, "m_out")
lambda_vec  <- get_param_vec(params_al_raw, "lambda")

compute_mcr <- function(s) {
  id      <- uIds[s]
  s_df    <- filter(fulldf,  subID == id)
  s_train <- filter(traindf, subID == id)
  nT      <- nrow(s_df)
  nTr     <- nrow(s_train)

  m_in_s  <- m_in_vec[s]
  m_out_s <- m_out_vec[s]
  lam_s   <- lambda_vec[s]

  GPin  <- plogis( m_in_s  * (s_train$selfResp - 4))
  GPout <- plogis(-m_out_s * (s_train$selfResp - 4))

  PS <- simMat[s_df$Idx, s_train$Idx]^lam_s  # nT × nTr

  simW_in  <- as.numeric(PS %*% GPin)  + 1e-9
  simW_out <- as.numeric(PS %*% GPout) + 1e-9

  mean(simW_in / simW_out)
}

mcr_vals <- sapply(1:nSubjects, compute_mcr)
mcr <- data.frame(subj_idx = 1:nSubjects, subject_mcr = mcr_vals)

cat(sprintf("\nsubject_mcr summary (computed analytically):\n"))
print(summary(mcr$subject_mcr))

# ── 5. Merge with structural params and MCR ───────────────────────────────────
params_sl <- read.csv(here("Results", "params_ind_s1_sym_lambda.csv"))
params_sl <- params_sl[grepl("^(m|bias|lambda)\\[", params_sl$variable), ]
params_sl$subj_idx <- as.integer(regmatches(params_sl$variable, regexpr("[0-9]+", params_sl$variable)))
params_sl$param    <- sub("\\[.*", "", params_sl$variable)
params_wide <- reshape(params_sl[, c("subj_idx","param","median")],
  idvar = "subj_idx", timevar = "param", direction = "wide")
names(params_wide) <- sub("median\\.", "", names(params_wide))

id_df <- fulldf[fulldf$subID %in% common_ids, ]
id_df <- id_df[!duplicated(id_df$subID), ]
id_df <- id_df[order(match(id_df$subID, uIds)), ]
id_df$subj_idx <- seq_len(nrow(id_df))
id_df <- id_df[, c("subj_idx","subID","DS","Proto","SCC","SI","RSE","NTB","NFC","SING.Ind","SING.Inter")]

# Merge everything
df <- merge(params_wide,                          mcr,        by = "subj_idx")
df <- merge(df, subj_elpd[, c("subj_idx","delta_elpd")],     by = "subj_idx")
df <- merge(df,                                   id_df,      by = "subj_idx")
cat(sprintf("\nMerged N: %d\n", nrow(df)))

# ── 6. Correlate params × individual difference scales ────────────────────────
# Two families with separate FDR correction:
#   Family A (mechanism-personality): m, bias/γ, lambda, subject_mcr × scales (k=36)
#   Family B (ΔELPD-personality): delta_elpd × scales — separate question (k=9)
scale_cols  <- c("DS","Proto","SCC","SI","RSE","NTB","NFC","SING.Ind","SING.Inter")
mech_params <- intersect(c("m","bias","lambda","subject_mcr"), names(df))
elpd_params <- "delta_elpd"

run_cors <- function(param_set, df, scale_cols) {
  rows <- list()
  for (p in param_set) for (s in scale_cols) {
    ct <- cor.test(df[[p]], df[[s]], use = "complete.obs")
    rows[[length(rows)+1]] <- data.frame(
      param = p, scale = s,
      r     = round(ct$estimate, 3),
      p_raw = round(ct$p.value, 4),
      n     = sum(complete.cases(df[, c(p, s)]))
    )
  }
  do.call(rbind, rows)
}

res_mech <- run_cors(mech_params, df, scale_cols)
res_elpd <- run_cors(elpd_params, df, scale_cols)

res_mech$p_fdr <- round(p.adjust(res_mech$p_raw, method = "BH"), 4)
res_elpd$p_fdr <- round(p.adjust(res_elpd$p_raw, method = "BH"), 4)

res <- rbind(res_mech, res_elpd)
res <- res[order(res$p_raw), ]

cat("\n=== ΔELPD (asym - sym) correlations with individual differences ===\n")
print(res_elpd[order(res_elpd$p_raw), ], row.names = FALSE)

cat("\n=== All nominally significant (p_raw < .05) — mechanism params ===\n")
print(res_mech[res_mech$p_raw < .05, ], row.names = FALSE)

cat("\n=== FDR-significant mechanism params (p_fdr < .05) ===\n")
print(res_mech[!is.na(res_mech$p_fdr) & res_mech$p_fdr < .05, ], row.names = FALSE)

cat("\n=== FDR-significant ΔELPD (p_fdr < .05) ===\n")
print(res_elpd[!is.na(res_elpd$p_fdr) & res_elpd$p_fdr < .05, ], row.names = FALSE)

# ── 7. Save extended individual-differences dataset ───────────────────────────
write.csv(df,  here("Results", "ind_diffs_s1_full.csv"),    row.names = FALSE)
write.csv(res, here("Results", "correlations_s1_full.csv"), row.names = FALSE)
message("\nSaved: Results/ind_diffs_s1_full.csv, Results/correlations_s1_full.csv")

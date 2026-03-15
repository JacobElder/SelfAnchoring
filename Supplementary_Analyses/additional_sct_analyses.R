# Additional SCT Analyses
# 1. Depersonalization: selfResp × SI/Proto/SING.Ind interaction in GLMM (Studies 1–3)
# 2. Positive distinctiveness: γ–outgroup warmth Spearman (Study 2, High-Status)
# 3. λ–SCC: already in correlations_s3_full.csv — documented here for reference

library(lme4)
library(lmerTest)
library(dplyr)
library(here)

# ── HELPER ────────────────────────────────────────────────────────────────────
z <- function(x) as.numeric(scale(x))

run_depersonalization <- function(df, study_label) {
  # df must have: ingChoiceN, selfResp, desirability, SI, Proto, SING.Ind, subID, trait
  df <- df %>% filter(!is.na(ingChoiceN), !is.na(selfResp))

  # Z-score within dataset
  df$selfResp.Z   <- z(df$selfResp)
  df$desir.Z      <- z(df$desirability)
  df$SI.Z         <- z(df$SI)
  df$Proto.Z      <- z(df$Proto)
  df$SINGInd.Z    <- z(df$SING.Ind)

  results <- list()
  # Try multiple optimizers in order; use first that converges cleanly
  optimizers <- c("bobyqa", "nlminbwrap", "Nelder_Mead")

  for (mod_name in c("SI.Z", "Proto.Z", "SINGInd.Z")) {
    df$mod <- df[[mod_name]]
    fit <- NULL
    for (opt in optimizers) {
      fit <- tryCatch(
        suppressWarnings(
          glmer(ingChoiceN ~ selfResp.Z * mod + desir.Z + (1 + selfResp.Z | subID) + (1 | trait),
                data = df, family = binomial,
                control = glmerControl(optimizer = opt, optCtrl = list(maxfun = 2e5)))
        ),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        warns <- fit@optinfo$conv$lme4$messages
        if (is.null(warns) || length(warns) == 0) break  # clean convergence
      }
    }
    if (!is.null(fit)) {
      coef_tbl <- summary(fit)$coefficients
      row_name  <- "selfResp.Z:mod"
      if (row_name %in% rownames(coef_tbl)) {
        r <- coef_tbl[row_name, ]
        results[[mod_name]] <- data.frame(
          study       = study_label,
          moderator   = mod_name,
          beta        = round(r["Estimate"], 4),
          se          = round(r["Std. Error"], 4),
          z_val       = round(r["z value"], 3),
          p_val       = round(r["Pr(>|z|)"], 4)
        )
      }
    }
  }
  do.call(rbind, results)
}

# ── STUDY 1 ───────────────────────────────────────────────────────────────────
message("Study 1 depersonalization interactions...")
s1_test <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))
s1_indiff <- s1_test %>% group_by(subID) %>% slice(1) %>%
  select(subID, SI, Proto, SING.Ind)
s1_df <- s1_test %>%
  select(subID, trait, ingChoiceN, selfResp, desirability) %>%
  left_join(s1_indiff, by = "subID")
res_s1 <- run_depersonalization(s1_df, "Study1")

# ── STUDY 2 ───────────────────────────────────────────────────────────────────
message("Study 2 depersonalization interactions...")
s2_test <- read.csv(here("Study 2/Cleaning/output/fullTest.csv"))
# SI: check column name
si_col <- if ("SI" %in% names(s2_test)) "SI" else grep("^SI$|MGIS", names(s2_test), value=TRUE)[1]
s2_indiff <- s2_test %>% group_by(subID) %>% slice(1) %>%
  select(subID, all_of(si_col), Proto, SING.Ind, outgroup, InUCLATherm, InCSULATherm)
if (si_col != "SI") s2_indiff <- rename(s2_indiff, SI = all_of(si_col))

s2_df <- s2_test %>%
  select(subID, trait, ingChoiceN, selfResp, desirability) %>%
  left_join(s2_indiff, by = "subID")
res_s2 <- run_depersonalization(s2_df, "Study2")

# ── STUDY 2: γ–OUTGROUP WARMTH (Positive Distinctiveness) ─────────────────────
message("Study 2: γ–outgroup warmth Spearman (High-Status condition)...")

params_s2 <- read.csv(here("Results/params_ind_s2_sym_lambda.csv"))
# bias[i] = γ for subject i
bias_rows <- params_s2 %>%
  filter(grepl("^bias\\[", variable)) %>%
  mutate(subj_idx = as.integer(gsub("bias\\[|\\]", "", variable)))

# Map subj_idx to subID using subject ordering (sorted unique IDs from fullTest)
uIds_s2 <- sort(unique(s2_test$subID))
bias_rows$subID <- uIds_s2[bias_rows$subj_idx]

# Subject-level data: condition + outgroup warmth
subj_s2 <- s2_test %>%
  group_by(subID) %>% slice(1) %>%
  select(subID, outgroup, InUCLATherm, InCSULATherm)

gamma_warmth <- bias_rows %>%
  select(subID, gamma = median) %>%
  left_join(subj_s2, by = "subID") %>%
  filter(!is.na(outgroup))

# High-Status = UCLA outgroup
hs_df <- gamma_warmth %>% filter(outgroup == "UCLA") %>%
  filter(!is.na(InUCLATherm), !is.na(gamma))

message(sprintf("High-Status N = %d", nrow(hs_df)))
sp_hs <- cor.test(hs_df$gamma, hs_df$InUCLATherm, method = "spearman", exact = FALSE)

# Low-Status = CSULA (for comparison)
ls_df <- gamma_warmth %>% filter(outgroup == "CSU LA") %>%
  filter(!is.na(InCSULATherm), !is.na(gamma))
sp_ls <- cor.test(ls_df$gamma, ls_df$InCSULATherm, method = "spearman", exact = FALSE)

# Negation condition (InUCLATherm not applicable; skip warmth for negation)
warmth_results <- data.frame(
  condition     = c("High-Status (UCLA)", "Low-Status (CSULA)"),
  n             = c(nrow(hs_df), nrow(ls_df)),
  spearman_rho  = round(c(sp_hs$estimate, sp_ls$estimate), 3),
  p_val         = round(c(sp_hs$p.value, sp_ls$p.value), 4)
)

message("\nγ–outgroup warmth Spearman results:")
print(warmth_results)

# ── STUDY 3 ───────────────────────────────────────────────────────────────────
message("Study 3 depersonalization interactions...")
s3_test <- read.csv(here("Study 3/Cleaning/output/fullTest_fixed.csv"))
s3_indiff <- s3_test %>% group_by(subID) %>% slice(1) %>%
  select(subID, SI, Proto, SING.Ind)
s3_df <- s3_test %>%
  select(subID, trait, ingChoiceN, selfResp, desirability) %>%
  left_join(s3_indiff, by = "subID")
res_s3 <- run_depersonalization(s3_df, "Study3")

# ── SAVE RESULTS ──────────────────────────────────────────────────────────────
deperso_results <- rbind(res_s1, res_s2, res_s3)
message("\nDepersonalization interaction results:")
print(deperso_results)

write.csv(deperso_results,
          here("Results/depersonalization_interactions.csv"),
          row.names = FALSE)
write.csv(warmth_results,
          here("Results/gamma_warmth_s2.csv"),
          row.names = FALSE)

# ── λ–SCC REFERENCE (Study 3, already in correlations_s3_full.csv) ────────────
message("\nλ–SCC (Study 3, from existing correlations):")
corrs3 <- read.csv(here("Results/correlations_s3_full.csv"), header = FALSE,
                   col.names = c("param", "scale", "r", "p_raw", "n", "p_fdr"))
lambda_scc <- corrs3 %>% filter(param == "lambda", scale == "SCC")
print(lambda_scc)

message("\nAll supplementary SCT analyses complete.")

# Parameter × Individual Differences Correlations — Study 1 sym_lambda
# All model parameters vs all ID measures; FDR-corrected (Benjamini-Hochberg)

library(cmdstanr); library(tidyverse); library(here); library(posterior)

# ── 1. Reconstruct subject order ─────────────────────────────────────────────
fulldf  <- read.csv(here("Study 1/Cleaning/output/fullTest.csv"))  %>% filter(!is.na(ingChoiceN))
traindf <- read.csv(here("Study 1/Cleaning/output/fullTrain.csv")) %>% filter(!is.na(selfResp))
common_ids <- intersect(unique(fulldf$subID), unique(traindf$subID))
uIds <- sort(common_ids)
n    <- length(uIds)

# ── 2. Extract posterior medians for all individual-level parameters ──────────
fit <- readRDS(here("Fits", "fit_s1_sym_lambda.rds"))

extract_median <- function(param) {
  draws <- as_draws_df(fit$draws(param))
  cols  <- paste0(param, "[", seq_len(n), "]")
  sapply(cols, function(col) median(draws[[col]]))
}

params_df <- data.frame(
  subID       = uIds,
  tau         = extract_median("tau"),
  m           = extract_median("m"),
  bias        = extract_median("bias"),
  lambda      = extract_median("lambda"),
  w           = extract_median("w"),
  subject_mcr = extract_median("subject_mcr")
)

cat("Parameters extracted. Summary:\n")
print(round(sapply(params_df[-1], summary), 3))

# ── 3. Load individual differences ───────────────────────────────────────────
# Psych scales from fullTest (one row per subject)
id_psych <- fulldf %>%
  select(subID, Proto, SING.Ind, SING.Inter, SCC, DS, SI, RSE, NTB, NFC) %>%
  group_by(subID) %>% slice(1) %>% ungroup()

# Demographics from preQs
preqs_path <- here("Study 1/Cleaning/input/Qualtrics Data/SA1preQs.csv")
id_demo <- read.csv(preqs_path) %>%
  rename(subID = id) %>%
  select(subID, Age, Gender, Polit, SStatus, GStatus, SocClass) %>%
  mutate(across(c(Age, Polit, SStatus, GStatus, SocClass), as.numeric))

# Merge everything
df <- params_df %>%
  left_join(id_psych, by = "subID") %>%
  left_join(id_demo,  by = "subID")

id_vars    <- c("Proto","SING.Ind","SING.Inter","SCC","DS","SI","RSE","NTB","NFC",
                "Age","Gender","Polit","SStatus","GStatus","SocClass")
param_vars <- c("tau","m","bias","lambda","w","subject_mcr")

cat("\nN subjects with complete data per variable:\n")
print(sapply(c(param_vars, id_vars), function(v) sum(!is.na(df[[v]]))))

# ── 4. All pairwise correlations ─────────────────────────────────────────────
results <- expand.grid(param = param_vars, id_var = id_vars, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    x        = list(df[[param]]),
    y        = list(df[[id_var]]),
    complete = sum(complete.cases(df[[param]], df[[id_var]])),
    r        = if (complete >= 5 && sd(y, na.rm=TRUE) > 0)
                 cor.test(x, y, method = "pearson")$estimate else NA_real_,
    p_raw    = if (complete >= 5 && sd(y, na.rm=TRUE) > 0)
                 cor.test(x, y, method = "pearson")$p.value  else NA_real_
  ) %>%
  ungroup() %>%
  select(-x, -y) %>%
  filter(!is.na(p_raw))

# FDR correction across all tests
results$p_fdr <- p.adjust(results$p_raw, method = "BH")

results <- results %>%
  mutate(
    sig_raw = case_when(p_raw < .001 ~ "***", p_raw < .01 ~ "**",
                        p_raw < .05  ~ "*",   p_raw < .10 ~ ".", TRUE ~ ""),
    sig_fdr = case_when(p_fdr < .001 ~ "***", p_fdr < .01 ~ "**",
                        p_fdr < .05  ~ "*",   p_fdr < .10 ~ ".", TRUE ~ "")
  ) %>%
  arrange(p_fdr)

# ── 5. Full table ─────────────────────────────────────────────────────────────
cat("\n\n=== ALL CORRELATIONS (sorted by FDR-corrected p) ===\n")
cat(sprintf("%-14s  %-14s  %6s  %7s  %5s  %7s  %5s\n",
            "Parameter", "ID Variable", "r", "p_raw", "sig", "p_fdr", "sig_fdr"))
cat(strrep("-", 70), "\n")

for (i in seq_len(nrow(results))) {
  row <- results[i, ]
  cat(sprintf("%-14s  %-14s  %6.3f  %7.4f  %5s  %7.4f  %5s\n",
              row$param, row$id_var, row$r, row$p_raw, row$sig_raw, row$p_fdr, row$sig_fdr))
}

cat("\nNote: . p<.10  * p<.05  ** p<.01  *** p<.001 (sig_fdr = BH-corrected)\n")

# ── 6. Notable findings summary ───────────────────────────────────────────────
cat("\n\n=== SURVIVES FDR CORRECTION (p_fdr < .05) ===\n")
sig <- results %>% filter(p_fdr < .05)
if (nrow(sig) == 0) {
  cat("None survive FDR at p < .05\n")
} else {
  print(sig %>% select(param, id_var, r, p_raw, p_fdr, sig_fdr))
}

cat("\n=== MARGINAL AFTER FDR (p_fdr < .10) ===\n")
marg <- results %>% filter(p_fdr >= .05, p_fdr < .10)
if (nrow(marg) == 0) cat("None\n") else
  print(marg %>% select(param, id_var, r, p_raw, p_fdr, sig_fdr))

# ── 7. Heatmap-style matrix for visual summary ────────────────────────────────
cat("\n=== CORRELATION MATRIX (r values; * = p_raw<.05, + = survives FDR<.05) ===\n")
mat_r   <- matrix(NA, nrow=length(param_vars), ncol=length(id_vars),
                  dimnames=list(param_vars, id_vars))
mat_sig <- matrix("", nrow=length(param_vars), ncol=length(id_vars),
                  dimnames=list(param_vars, id_vars))

for (i in seq_len(nrow(results))) {
  p <- results$param[i]; v <- results$id_var[i]
  mat_r[p, v]   <- results$r[i]
  marker <- if (results$p_fdr[i] < .05) "+" else if (results$p_raw[i] < .05) "*" else ""
  mat_sig[p, v] <- marker
}

for (p in param_vars) {
  cat(sprintf("%-12s", p))
  for (v in id_vars) {
    r_val <- mat_r[p, v]
    mk    <- mat_sig[p, v]
    if (is.na(r_val)) cat(sprintf(" %7s", "NA"))
    else cat(sprintf(" %5.2f%-2s", r_val, mk))
  }
  cat("\n")
}
cat(sprintf("%-12s", ""))
for (v in id_vars) cat(sprintf(" %-7s", substr(v, 1, 7)))
cat("\n* = p_raw<.05 (uncorrected)  + = p_fdr<.05 (BH-corrected)\n")

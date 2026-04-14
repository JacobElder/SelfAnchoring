# SelfAnchoring

Research investigating how people project self-perceptions onto their ingroup using relational similarity and contrast-based mechanisms. Computational models formalize self-anchoring as similarity-based generalization across a semantic trait network.

**Design:** Participants self-evaluate on 90 traits (training phase), then classify each of 148 traits as typical of their ingroup or outgroup (generalization phase).

- **Study 1:** Minimal groups (overestimator/underestimator; *N* = 61)
- **Study 2:** University groups — UCR (ingroup) vs. UCLA (high-status), CSU LA (low-status), or Not UCR (negation); *N* = 181
- **Study 3:** Racial groups — Asian/Latino students contrasting against opposite minority or White majority; *N* = 265

---

## Paper Section → Script Mapping

| Paper Section | Script(s) | Output(s) |
|:---|:---|:---|
| **Study 1 — Model Comparison (Table 2)** | `Study 1/Analysis/run_model_comparison_s1.R` | `Results/model_comparison_results_s1.csv`, `Fits/loo_s1_*.rds` |
| **Study 1 — Behavioral GLMMs** | `Study 1/Analysis/marginal_effects_s1.R` | `Results/marginal_effects_s1_ames.csv` |
| **Study 1 — Individual Differences** | `Study 1/Analysis/loo_delta_indiff_s1.R` | `Results/ind_diffs_s1_full.csv`, `Results/correlations_s1_full.csv` |
| **Study 1 — Trait Segregation (Supp K)** | `Study 1/Analysis/trait_segregation_s1.R` | `Results/trait_segregation_s1_correlations.csv` |
| **Study 1 — CCA (Supp E)** | `Study 1/Analysis/cca_s1.R` | `Results/cca_s1_correlations.csv`, `Results/cca_s1_loadings.csv` |
| **Study 2 — Model Comparison (Table 3)** | `Study 2/Analysis/run_model_comparison_s2.R` | `Results/model_comparison_results_s2.csv`, `Fits/loo_s2_*.rds` |
| **Study 2 — ELPD by Condition** | `Study 2/Analysis/elpd_by_condition_s2.R` | `Results/elpd_by_condition_s2_summary.csv` |
| **Study 2 — Parameters by Condition** | `Study 2/Analysis/param_comparison_conditions_s2.R` | `Results/param_comparison_conditions_s2.csv` |
| **Study 2 — Behavioral GLMMs** | `Study 2/Analysis/marginal_effects_s2.R` | `Results/marginal_effects_s2_ames.csv` |
| **Study 2 — Individual Differences** | `Study 2/Analysis/loo_delta_indiff_s2.R` | `Results/ind_diffs_s2_full.csv`, `Results/correlations_s2_full.csv` |
| **Study 2 — Trait Segregation (Supp K)** | `Study 2/Analysis/trait_segregation_s2.R` | `Results/trait_segregation_s2_correlations.csv` |
| **Study 2 — CCA (Supp E)** | `Study 2/Analysis/cca_s2.R` | `Results/cca_s2_correlations.csv`, `Results/cca_s2_loadings.csv` |
| **Study 3 — Model Comparison (Table 4)** | `Study 3/Analysis/run_model_comparison_s3.R` | `Results/model_comparison_results_s3.csv`, `Fits/loo_s3_*.rds` |
| **Study 3 — ELPD by Condition** | `Study 3/Analysis/elpd_by_condition_s3.R` | `Results/elpd_by_condition_s3_summary.csv` |
| **Study 3 — Parameters by Condition** | `Study 3/Analysis/param_comparison_conditions_s3.R` | — |
| **Study 3 — Behavioral GLMMs** | `Study 3/Analysis/marginal_effects_s3.R` | `Results/marginal_effects_s3_ames.csv` |
| **Study 3 — Individual Differences** | `Study 3/Analysis/loo_delta_indiff_s3.R` | `Results/ind_diffs_s3_full_enriched.csv`, `Results/correlations_s3_full.csv` |
| **Study 3 — Trait Segregation (Supp K)** | `Study 3/Analysis/trait_segregation_s3.R` | `Results/trait_segregation_s3_correlations.csv` |
| **Study 3 — CCA (Supp E)** | `Study 3/Analysis/cca_s3.R` | `Results/cca_s3_correlations.csv`, `Results/cca_s3_loadings.csv` |
| **Study 3 — Warmth Analysis (Supp H)** | `Study 3/Analysis/warmth_analysis_s3.R` | `Results/warmth_correlations_s3.csv`, `Results/warmth_lm_s3.csv` |
| **Study 3 — Race Moderation (Supp C)** | `Study 3/Analysis/race_moderation_s3.R` | `Results/race_moderation_params_s3.csv` |
| **Pooled — Model Comparison (Table 5)** | `Pooled/run_pooled_analysis.R` | `Results/pooled_model_comparison_results.csv`, `Fits/loo_pooled_*.rds` |
| **Pooled — Behavioral GLMMs (pooled)** | `Pooled/marginal_effects_pooled.R` | `Results/marginal_effects_pooled_ames.csv` |
| **Pooled — Individual Differences** | `Pooled/pooled_param_indiff_correlations.R` | `Results/correlations_pooled_standardized.csv` |
| **Pooled — CCA (Supp E)** | `Pooled/cca_pooled.R` | `Results/cca_pooled_correlations.csv`, `Results/cca_pooled_loadings.csv` |
| **Parameter Recovery (Supp B)** | `Parameter Recovery/run_parameter_recovery.R` | `Results/parameter_recovery_sym_lambda_results.csv`, `Results/parameter_recovery_asym_lambda_results.csv` |
| **Empirical Recovery (Supp B)** | `Parameter Recovery/run_parameter_recovery_empirical.R` | — |
| **Depersonalization (Supp I)** | `Supplementary_Analyses/additional_sct_analyses.R` | `Results/depersonalization_interactions.csv` |
| **Entropy (Supp J)** | `Scripts/entropy_marginal_effects.R` | `Results/entropy_ames.csv` |
| **All Figures** | `Figures/make_figures.R` | `Figures/*.tiff` (git-ignored) |
| **Forest Plot (Fig 8)** | `Figures/make_forest_plot.R` | `Figures/fig_forest_plot_pooled.tiff` |
| **Recovery Figures (Supp B)** | `Figures/make_param_recovery_fig.R` | `Figures/fig_param_recovery*.tiff` |

---

## Directory Structure

```
SelfAnchoring/
├── Manuscript/
│   └── ElderJacob_JPSP_Submission.qmd  ← source of truth; knits to .docx via apaquarto
├── Computational Models/
│   ├── S_Bias_NoW.stan                 ← 1-param baseline
│   ├── S_Symmetric_NoW.stan            ← 2-param (m, γ)
│   ├── S_Sym_Lambda_NoW.stan           ← 3-param (m, γ, λ) — WINNER all studies
│   ├── S_Asym_Lambda_NoW.stan          ← 4-param (m_in, m_out, γ, λ) — comparison
│   └── Archive/                        ← legacy models (lapse parameter w; superseded)
├── Study 1/
│   ├── Analysis/                       ← canonical scripts (see mapping above)
│   └── Cleaning/                       ← cleaning Rmd + output/fullTrain.csv, fullTest.csv
├── Study 2/
│   ├── Analysis/                       ← canonical scripts
│   └── Cleaning/                       ← cleaning Rmd + output/fullTrain_fixed.csv, fullTest.csv
├── Study 3/
│   ├── Analysis/                       ← canonical scripts
│   └── Cleaning/                       ← cleaning Rmd + output/fullTrain_fixed.csv, fullTest_fixed.csv
├── Pooled/
│   ├── run_pooled_analysis.R           ← core 3-level Stan model
│   ├── marginal_effects_pooled.R       ← behavioral GLMMs
│   ├── cca_pooled.R                    ← cross-study CCA
│   ├── pooled_param_indiff_correlations.R
│   ├── input/adjacencyMatrix_p.csv     ← 148-trait network (shared across all studies)
│   └── [legacy .Rmd notebooks]         ← pooledAnalyses, demographics, plotting, etc.
├── Parameter Recovery/
│   ├── run_parameter_recovery.R        ← random-parameter recovery (two-process)
│   ├── run_parameter_recovery_empirical.R ← empirical recovery
│   └── [legacy .Rmd notebooks]         ← parameterRecoveryScript, crossValidation, etc.
├── Figures/
│   ├── make_figures.R                  ← master figure script
│   ├── make_forest_plot.R
│   ├── make_param_recovery_fig.R
│   └── [generated .tiff files]         ← git-ignored; run make_figures.R to generate
├── Fits/                               ← LOO .rds files + pooled_stan_data.rds (retained)
├── Results/                            ← all analysis output CSVs
├── Scripts/                            ← utility scripts (entropy, export, debug)
├── Supplementary_Analyses/             ← additional SCT tests (depersonalization, etc.)
├── Shiny/app.R                         ← interactive supplement (deployable)
└── tests/
    ├── run_tests.R
    └── testthat/                       ← test-data-integrity.R, test-model-helpers.R, etc.
```

---

## Run Commands

All commands should be run from the project root (`/Users/jacobelder/Documents/GitHub/SelfAnchoring`).

### Step 1: Model Fitting

**Full pipeline (all studies, sequential, plays audio when done):**
```bash
caffeinate -i sh -c 'Rscript "Study 1/Analysis/run_model_comparison_s1.R" && Rscript "Study 2/Analysis/run_model_comparison_s2.R" && Rscript "Study 3/Analysis/run_model_comparison_s3.R" && Rscript "Pooled/run_pooled_analysis.R" && say "All updated models finally finished"'
```

**Individual studies:**
```bash
Rscript "Study 1/Analysis/run_model_comparison_s1.R"
Rscript "Study 2/Analysis/run_model_comparison_s2.R"
Rscript "Study 3/Analysis/run_model_comparison_s3.R"
```

**Pooled analysis only:**
```bash
Rscript "Pooled/run_pooled_analysis.R"
```

### Step 2: Secondary Analyses (after model fitting)

Run these after the corresponding `run_model_comparison_s*.R` has completed:

```bash
# Individual differences + correlations (each study)
Rscript "Study 1/Analysis/loo_delta_indiff_s1.R"
Rscript "Study 2/Analysis/loo_delta_indiff_s2.R"
Rscript "Study 3/Analysis/loo_delta_indiff_s3.R"

# ELPD by condition
Rscript "Study 2/Analysis/elpd_by_condition_s2.R"
Rscript "Study 3/Analysis/elpd_by_condition_s3.R"

# Parameter condition comparisons
Rscript "Study 2/Analysis/param_comparison_conditions_s2.R"
Rscript "Study 3/Analysis/param_comparison_conditions_s3.R"

# Behavioral GLMMs
Rscript "Study 1/Analysis/marginal_effects_s1.R"
Rscript "Study 2/Analysis/marginal_effects_s2.R"
Rscript "Study 3/Analysis/marginal_effects_s3.R"
Rscript "Pooled/marginal_effects_pooled.R"

# CCA
Rscript "Study 1/Analysis/cca_s1.R"
Rscript "Study 2/Analysis/cca_s2.R"
Rscript "Study 3/Analysis/cca_s3.R"
Rscript "Pooled/cca_pooled.R"
Rscript "Pooled/pooled_param_indiff_correlations.R"

# Supplementary (Study 3 specific)
Rscript "Study 3/Analysis/warmth_analysis_s3.R"
Rscript "Study 3/Analysis/race_moderation_s3.R"

# Trait segregation (all studies)
Rscript "Study 1/Analysis/trait_segregation_s1.R"
Rscript "Study 2/Analysis/trait_segregation_s2.R"
Rscript "Study 3/Analysis/trait_segregation_s3.R"
```

### Step 3: Figures

```bash
# All manuscript figures (requires all secondary analyses to be complete)
Rscript "Figures/make_figures.R"
Rscript "Figures/make_forest_plot.R"
Rscript "Figures/make_param_recovery_fig.R"
```

### Parameter Recovery

```bash
# Two-process: setup first (frees memory), then fit
caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery.R" setup && Rscript "Parameter Recovery/run_parameter_recovery.R"'

# Empirical recovery
caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery_empirical.R" setup && Rscript "Parameter Recovery/run_parameter_recovery_empirical.R"'
```

### Render Manuscript

```bash
quarto render Manuscript/ElderJacob_JPSP_Submission.qmd
```

---

## Notes

- `caffeinate -i` prevents macOS from sleeping during long-running fits
- Seeds: S1=123, S2=456, S3=789, Pooled=1234, Recovery sym=999, Recovery asym=1000
- Full fit objects (`.rds`) are NOT saved — only LOO objects, summary CSVs, and param CSVs
- **Keep:** `Fits/loo_s*.rds`, `Fits/loo_pooled_*.rds`, `Fits/pooled_stan_data.rds`
- **Never save:** `fit_s*.rds` (1–1.5 GB each; not needed; `save_object()` is disabled)
- Model outputs → `Results/`; LOO objects → `Fits/`

### Data Integrity Notes

| Study | Training data | Test data |
|:------|:-------------|:---------|
| S1 | `Study 1/Cleaning/output/fullTrain.csv` | `Study 1/Cleaning/output/fullTest.csv` |
| S2 | `Study 2/Cleaning/output/fullTrain_fixed.csv` | `Study 2/Cleaning/output/fullTest.csv` |
| S3 | `Study 3/Cleaning/output/fullTrain_fixed.csv` | `Study 3/Cleaning/output/fullTest_fixed.csv` |

S2/S3 `_fixed` files correct original binary Parquet files mislabeled as .csv. S3 test data was also corrected. Scripts hardcode the correct filenames.

---

## Legacy / Residual Files

These are kept for reference but are **not part of the active pipeline**:

| Location | Files | Status |
|:---------|:------|:-------|
| `Study 1/Analysis/` | `SAs1_Analysis.Rmd`, `analyze_latent_s1.R`, `analyze_latent_s1_lambda.R`, `PowerAnalysis_2022_4_24_JE.Rmd` | Superseded by modular .R scripts |
| `Study 2/Analysis/` | `SAs2_Analysis.Rmd`, `PowerAnalysis.Rmd`, `rescue_s2_asym_lambda.R`, `rescue_s2_minimal.R`, `get_stats_s2.R`, `a_priori_power_s2.R` | Legacy or one-time utilities |
| `Study 3/Analysis/` | `SAs3_Analysis.Rmd` | Superseded |
| `Pooled/` | `pooledAnalyses.Rmd`, `demographics.Rmd`, `plotting.Rmd`, `plotting2.Rmd`, `modelDepiction.Rmd`, `reportingAnalyses.qmd`, `prepare_pooled_data.R`, `run_pooled_compare_only.R`, `run_pooled_one_model.R`, `run_pooled_resume.R` | Legacy notebooks and utility scripts |
| `Parameter Recovery/` | `parameterRecoveryScript.Rmd`, `crossValidation.Rmd`, `parameter_recovery_template.R`, `PR_S_Logistic_1mOppose_Bias.R` | Legacy; superseded by run_parameter_recovery.R |
| `Computational Models/Archive/` | All `S_*_w.stan` and scratch scripts | Old lapse-parameter models; removed due to identifiability failure |

---

## TODO

- [ ] **Task #1 (Stan):** Refit S1 `asym_lambda` to produce Stan-native `subject_mcr[i]` output for supplementary comparisons. S_Asym_Lambda_NoW.stan was updated 2026-03-23 to output `subject_mcr`; current S1 fit predates this update. Run `Study 1/Analysis/run_asym_lambda_only_s1.R`.
- [ ] **CCA (S2/S3):** Study-level CCA results in Supp E may be stale (pre–sym-MCR). Re-run `cca_s2.R` and `cca_s3.R` with updated `ind_diffs_s2/s3_full_enriched.csv` to confirm loadings are current.
- [ ] **S1 CCA verification:** `cca_s1.R` outputs are not in a versioned CSV; verify CCA canonical correlation statistics (r = .720, p = .013 in Supp E) against current run.
- [ ] **Organize legacy files:** Consider moving all files listed in "Legacy / Residual Files" above to `Archive/` subdirectories within each study folder to keep `Analysis/` directories clean.
- [ ] **Shiny app deployment:** Confirm `Shiny/app.R` is deployed and URL in Supplementary F is live.
- [ ] **Tests:** Run `Rscript tests/run_tests.R` and fix any failures before next model refit.

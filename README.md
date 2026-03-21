# SelfAnchoring

Cleaning and Analyses for research investigating how people project self-perceptions onto their ingroup using relational similarity and contrast-based mechanisms.

Design: Self-evaluate on 90 traits. Then classify each trait as typical of ingroup or outgroup.

Study 1: Minimal groups... Randomly assign to overestimator or underestimator group.

Study 2: University of California Riverside Students (Ingroup) against UCLA (higher status), CSU LA (lower status), or Not UCLA (Negation/Control)

Study 3: Minority vs. Majority. Asian or Latino students exclusively recruited and compared against either Latino or Asian (respectively) or White.

---

## Run Commands

All commands should be run from the project root (`/Users/jacobelder/Documents/GitHub/SelfAnchoring`).

### Computational Model Fitting

**All studies (full pipeline — runs sequentially, plays audio when done):**
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

**Studies 2, 3, and Pooled only (e.g., after Study 1 is already done):**
```bash
caffeinate -i sh -c 'Rscript "Study 2/Analysis/run_model_comparison_s2.R" && Rscript "Study 3/Analysis/run_model_comparison_s3.R" && Rscript "Pooled/run_pooled_analysis.R" && say "Studies 2 3 Pooled finished"'
```

### Parameter Recovery

Parameter recovery validates that each winning model's architecture can reliably
recover known individual parameters from synthetic data. Uses Study 1 data structure
with simulated choices. Results saved to `Results/parameter_recovery_*_results.csv`.

**Run both models (S_Sym_Lambda and S_Asym_Lambda):**
```bash
caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery.R" setup && Rscript "Parameter Recovery/run_parameter_recovery.R"'
```

**Run one model only:**
```bash
caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery.R" setup && Rscript "Parameter Recovery/run_parameter_recovery.R" sym_lambda'
caffeinate -i sh -c 'Rscript "Parameter Recovery/run_parameter_recovery.R" setup && Rscript "Parameter Recovery/run_parameter_recovery.R" asym_lambda'
```

**Why two processes?** macOS OOM-kills a single R process because R never returns heap pages to the OS after `gc()`. Running `setup` first (base R + `here` only, ~20 MB) lets the OS fully reclaim memory before the fit process loads cmdstanr. The `setup` arg generates GMRF self-ratings, draws true parameters, simulates choices, and saves `Parameter Recovery/pr_cache.rds`; the fit process loads the cache and runs Stan.

### Notes
- `caffeinate -i` prevents macOS from sleeping during long-running fits
- Seeds: S1=123, S2=456, S3=789, Pooled=1234, Recovery sym=999, Recovery asym=1000
- Fit objects (`.rds`) are NOT saved — only LOO objects, summary CSVs, and param CSVs
- LOO files (`Fits/loo_*.rds`) and `Fits/pooled_stan_data.rds` are the only RDS files retained
- Model outputs go to `Results/`; LOO objects go to `Fits/`

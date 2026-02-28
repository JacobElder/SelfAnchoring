# Self-Anchoring Project Roadmap (JPSP)

This document tracks the ongoing analytical and manuscript tasks for the JPSP revision.

---

## 1. Ongoing Computational Work (External)
The following models are being fit by the user externally to avoid session timeouts:
*   **Study 2 Comparison:** `run_model_comparison_s2.R` (Symmetric, Asymmetric, Asym+Lambda, Full).
*   **Study 3 Comparison:** `run_model_comparison_s3.R` (Majority vs. Minority contexts).
*   **Pooled Hierarchical Analysis:** `run_pooled_analysis.R` (N=507, 3 studies) using `S_Pooled_Asym_Lambda.stan`.

**Task:** Once these finish, extract posterior medians and ∆LOO from the generated CSVs to populate manuscript placeholders.

---

## 2. Manuscript Refinement Strategy (ChatGPT Feedback)
We are integrating high-impact theoretical framing suggested by ChatGPT. **Mandate: Use surgical, paragraph-level edits. DO NOT perform large-scale deletions or structural overhauls.**

### Key Framing Moves:
1.  **Shepard’s Law Spine:** Position the invariance of generalization sensitivity ($\lambda$) across contexts as the primary "cognitive architecture" discovery.
2.  **Self-as-Prototype:** Frame the self not just as an exemplar, but as the *prototype* anchor in semantic space for novel/ambiguous groups.
3.  **Metacontrast Formalism:** Highlight the move from verbal theory to formal, model-based inference using network-derived relational similarity.
4.  **Parameter Defense:** Justify the low recoverability of $	au$ (temperature) as a nuisance parameter capturing choice stochasticity, which does not bias theoretical inferences ($\alpha, \lambda$).

---

## 3. Pending Analytical Innovations
*   **Parameter Recovery:** Once a winning model is selected (likely Asymmetric + Lambda), execute a formal recovery study using `parameter_recovery_template.R`.
*   **Network Hubs:** Test if trait centrality (indegree/outdegree) moderates generalization rates (Formalize: `Centrality^kappa`).
*   **Latent Marginal Effects:** Use posterior draws to calculate Average Marginal Effects (AME) for the impact of Social Identification on projection rates.

---

## 4. Technical Environment State
*   **Fixed Data:** Study 2 and 3 `fullTrain_fixed.csv` must be used (original `.csv` files were actually binary Parquet).
*   **Stan:** All models in `Computational Models/` are updated to modern `array[]` syntax and vectorized for matrix-vector multiplication.
*   **Portability:** Comparison scripts now export clean CSV summaries of medians and diagnostics to avoid RDS loading issues.

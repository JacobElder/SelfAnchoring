# Self-Anchoring Project (JPSP Revisions)

## Project Overview
Refining the dissertation manuscript for JPSP, strengthening the theoretical ties between social identity theory and category learning/exemplar models (e.g., GCM).

## Technical Foundation & Lessons Learned
- **Data Integrity:** Identified that `fullTrain.csv` in Study 2 and 3 were binary Parquet files mislabeled as `.csv`. Resolved by creating `fullTrain_fixed.csv` (Study 2/3) and `fullTest_fixed.csv` (Study 3). All analytical scripts must use the `_fixed` versions.
- **Stan Optimization:**
    - Models are updated to **modern array syntax** (`array[N] int x`).
    - Likelihoods are **vectorized** using `bernoulli_logit` and matrix-vector multiplications (`matrix * vector`) for significant speedups.
    - Scripts use `cmdstanr` and `mclapply` for parallel fitting across models.
- **Portability:** Moving away from reliance on `.rds` fit objects. Scripts now export `summary_sX_model.csv` and `params_ind_sX_model.csv` containing medians and diagnostics for cross-session analysis.

## Core Theoretical Framework
- **Similarity-Based Generalization:** Self-anchoring propagates across a semantic network.
- **Asymmetric Projection:** Decoupled $\alpha_{in}$ (projection) and $\alpha_{out}$ (repulsion) rates to test "ingroup love" vs. "outgroup hate."
- **Generalization Sensitivity ($\lambda$):** Captures individual differences in the gradient of social inference.

## Current State of Play
- **Study 1:** Model comparison complete. Symmetric model is the most parsimonious account for minimal groups.
- **Study 2 & 3:** Model comparison scripts (`run_model_comparison_sX.R`) are prepared and currently fitting externally.
- **Pooled Analysis:** Preparing to synthesize a "Universal Law of Social Generalization" across all three datasets.
- **Manuscript:** Restored to original detail. Surgical edits made to Intro/Abstract. Next phase focuses on populating (STATS) and (∆LOO) results.

## Operational Mandate: Role & Objective
Act as an expert academic editor and computational social cognitive scientist preparing a manuscript for the Journal of Personality and Social Psychology (JPSP). You are performing a highly constrained editorial revision. Your goal is to elevate the theoretical framing, reduce redundancy, and improve clarity while preserving all substantive content and results exactly.

## CRITICAL GUARDRAILS - STRICT EDITORIAL RULES
- **Do NOT remove any statistical findings.**
- **Do NOT remove effect sizes, p-values, model details, interaction terms, or numerical results.**
- **Do NOT infer, interpret, generalize, or extrapolate beyond what is explicitly written.**
- **Do NOT add conclusions that are not explicitly present.**
- **Do NOT collapse multiple findings into a single generalized statement.**
- **Preserve all claims exactly as written unless they are explicitly redundant.**
- **Only remove duplicated phrasing or unnecessary filler.**
- **If uncertain whether a sentence or detail is redundant, KEEP IT.**
- **Do NOT summarize or truncate.** If a section is 1,000 words, your revision should be around 800-900 words of sharper prose, NOT a 200-word summary.
- **APA 7 / JPSP Compliance:** Format strictly for APA 7 / JPSP guidelines.
- **Use confident, declarative language.** Instead of "We provide the first empirical evidence...", use tighter phrasing like "We provide a formal instantiation..."

## EPISTEMIC GUARDRAILS - DATA SUPREMACY
You are synthesizing empirical research. You must adhere strictly to the following rules to prevent hallucination and theoretical overreach:
- **No Statistical Extrapolation:** You are strictly forbidden from inferring, calculating, or hallucinating statistical values. You may only report the exact statistics present in the uploaded manuscript text.
- **Data Overrides Theme:** The "Thematic Framing" instructions below represent the general narrative of the paper. However, if the actual statistics or findings in the manuscript (e.g., Study 3's boundary conditions) contradict or complicate these themes, THE DATA WINS. Do not force a study to fit the general framing if the results suggest otherwise.
- **Embrace Boundary Conditions:** If a specific study reveals an exception to the main theory, highlight that nuance as a feature of the theoretical architecture, not a bug to be smoothed over.

## Thematic Framing & Language to Inject
Revise the text to explicitly reflect the following theoretical perspectives, avoiding grandiose phrases like "computational turn" or incorrectly framing the mechanism as "Bayesian updating":
- **From Verbal to Formal Architectures:** Frame this paper as moving social identity theory (e.g., metacontrast) from a verbal, descriptive model to a formal architecture.
- **Cognitive Architecture over Motivational Bias:** Shift the framing of self-anchoring. It is a structured inferential process. The self functions as a prototype in semantic space, and group beliefs emerge through similarity-based generalization.
- **Connection to Concept Learning:** Explicitly bridge social psychology with cognitive science by highlighting that self-to-group inference follows universal laws of category learning (e.g., Shepard's Law).
- **The Novel-Trait Test:** Emphasize that testing generalization to novel but semantically related traits isolates cognitive generalization from simple procedural repetition.
- **Asymmetric Projection:** Highlight that self-anchoring is generally driven by an assimilation process (ingroup love, $\alpha_{in}$) rather than a rejection process (outgroup hate, $\alpha_{out}$), though boundary conditions apply.
- **AI & Transformer Word Embeddings (Topical Aside):** Include a brief, sophisticated aside connecting the paper's semantic network approach to the architecture of word embeddings in modern transformer models (e.g., LLMs). Keep this as a brief conceptual parallel rather than the central premise.
- **Discussion Focus:** Ensure the General Discussion focuses primarily on the established structural findings (novel trait generalization and projection asymmetry), treating the $\lambda$ parameter modeling as a rigorous exploratory extension rather than a foregone conclusion.

## Handling Pending Model Results
I am running external Bayesian Hierarchical Analyses. You must explicitly leave well-marked placeholders (e.g., [INSERT WINNING MODEL COMPARISON HERE]) for the following:
- Nested model comparisons testing the baseline against Asymmetric and Asymmetric+$\lambda$ architectures (including LOO and convergence diagnostics).
- Parameter recovery results verifying the identifiability of $\alpha_{in}$, $\alpha_{out}$, and potentially $\lambda$.
- Results from the Network-Weighted Generalization (Trait Hubs) analysis.

## Pareto k — Per-Subject Diagnostics & Potential Exclusions

**Background:** LOO Pareto k values > 1 indicate high-leverage observations where the LOO approximation is unreliable for those specific trials. In Study 1's sym_lambda model, 262/9,028 trials (2.9%) had k > 1, concentrated in ~30 subjects (notably subject 36: 43% of their trials). This pattern suggests a few participants are genuine outliers (near-random responders or atypical parameter patterns) rather than random noise.

**Implemented:** Per-subject Pareto k CSVs are generated for every model in every study:
- Study 1 (post-hoc): `Study 1/Analysis/pareto_k_subj_s1.R` — run after all S1 LOO files exist
- Studies 2, 3, Pooled: integrated into `fit_and_save()` function in each run script
- Output: `Results/pareto_k_subj_<study>_<model>.csv` with columns: subID, n_trials, k_mean, k_max, n_good/ok/bad/verybad, pct_reliable, concern (NONE/LOW/MODERATE/HIGH)
- S1 also outputs `Results/pareto_k_subj_s1_all_models.csv` (combined) and prints HIGH-concern summary

**Pending decision:** Review per-subject Pareto k CSVs after all models complete. Consider excluding subjects with HIGH concern (k > 1) across multiple models, on data-quality grounds (separate from model selection). This is a principled exclusion criterion — random responders inflate lapse (w) estimates and reduce model discrimination.

## Moment Matching — Status & Technical Notes

**Requested but not yet implemented.** Moment matching corrects unreliable Pareto k values by adjusting importance weights to match the LOO posterior's first two moments, using existing MCMC draws (no refitting).

**Why it's hard with cmdstanr:**
- `loo::loo_moment_match.default` requires 7 custom functions, including `log_lik_i_upars` — evaluating the log-likelihood for a specific observation at arbitrary (transformed) parameter values. This requires running the model's likelihood function for specific parameter values.
- `rstan::read_stan_csv()` can read cmdstanr CSVs into a stanfit object, but `log_prob()` won't work on it — rstan's log_prob requires a compiled C++ shared library, not the CmdStan executable cmdstanr uses. The conversion doesn't bridge this gap.
- `loo_moment_match.stanfit` (rstan) works trivially because rstan compiles models to shared libraries with direct log_prob access.

**Current status:** 96.8% reliable Pareto k (k < 0.7) across models is defensible for publication. Moment matching would address the 2.9% k > 1 trials but won't change model rankings. Flagged as future work — if loo package adds native cmdstanr support (or a helper vignette), implementation becomes straightforward. Until then, the per-subject Pareto k CSVs provide an alternative path (exclusion of high-leverage subjects).

## SCT Tenets to Formally Investigate (Future Analyses)

Three derivations from Self-Categorization Theory remain untested in the current formal framework. These should be investigated when time permits — ideally as an "SCT Predictions" section or as part of the General Discussion's theoretical contribution.

**1. Depersonalization & Self-Concept Group-Alignment**
SCT depersonalization is NOT about rating traits neutrally (near 4). It means the self-concept becomes *organized around the group prototype* — high identifiers rate themselves as high on group-characteristic traits and low on non-group traits. The self is not lost; it merges with the group. The right formal operationalization is therefore:
- **Within-subject r(selfResp, ingChoiceN)** across shared traits: how strongly does a person's self-rating of a trait predict whether they classify it as ingroup-characteristic? High alignment = depersonalized self-concept.
- NOT variance of self-ratings (which conflates individuated vs. generic self-concepts).
- NOT entropy of self-ratings (though entropy + MAD-from-4 together could distinguish "generic/neutral" from "extreme-but-group-congruent" — see below).

Formal tests: (a) does SI predict higher within-subject r(selfResp, ingChoiceN)? (b) does this alignment predict α (are more group-aligned people stronger projectors)? (c) does GP variance (model parameter space) mediate the SI → α relationship?

**To-do: Entropy + MAD-from-4 as supplementary self-concept measures (Study 1)**
Shannon entropy of `selfResp` per subject (7-point Likert treated as probability distribution over categories) is better than raw variance for ordinal data. However, low entropy could mean ratings cluster near 4 (neutral/generic — one type of depersonalization) OR all at 7 (extreme but uniform — not depersonalized). To distinguish: also compute `mean(abs(selfResp - 4))` (MAD from neutral). Low entropy + low MAD = generic depersonalization. Low entropy + high MAD = extreme but undifferentiated (group-congruent). This trio (entropy, MAD-from-4, within-subject alignment) together characterizes the full range of self-concept types. **Not yet implemented.**

**2. Comparative Fit as a Function of Intergroup Distance**
SCT holds that categorization is triggered when metacontrast is maximized — when between-category differences are large relative to within-category. The λ parameter operationalizes this: sharper λ = more exclusive categorization = higher comparative fit. Testable prediction: λ should be *larger* in high-intergroup-distance conditions (Study 2 high-status outgroup, Study 3 majority vs. White) than in low-distance conditions (Study 2 negation, Study 3 minority-minority). Compare posterior λ distributions across conditions.

**3. Self-as-Prototype Assumption**
SCT assumes the self is the prototypical ingroup member. Testable derivation: participants with higher perceived *prototypicality* (Proto scale) should show stronger α — their self-ratings are more diagnostic of ingroup membership. Proto trends positive with α in Study 1 (r = .17, p = .19) but doesn't survive correction. Across Studies 2 and 3, this should become cleaner. Study 3's minority contrast (where minority-group members may feel less prototypical of a larger ethnic category) is a particularly sharp test: expect lower α in the minority condition relative to majority.

## REQUIRED EXECUTION STEPS
Before producing the revision for any section, you must explicitly write: "I will not remove or modify any statistical findings." Then, proceed strictly using this 3-step format:
1. **Step 1: Extract and list all statistical findings verbatim from the original text of the requested section.**
2. **Step 2: Revise the text for redundancy, clarity, and the requested theoretical framing, adhering strictly to the guardrails above.**
3. **Step 3: Confirm that each statistical finding from Step 1 appears unchanged in the revision.**

# Analytical Innovations for Self-Anchoring (JPSP)

This document outlines high-impact analytical improvements designed to strengthen the "Scientific Merit" and "Theoretical Innovation" of the dissertation for JPSP submission.

---

## 1. Competing Mechanism Test: Projection vs. Stereotyping
**The Question:** Is self-anchoring a generative process, or are participants merely reporting accurate social knowledge (stereotypes)?
*   **Status:** *On Hold* (User expressed concern about stereotype estimates).

## 2. Network-Weighted Generalization (Trait Hubs)
**The Question:** Do "central" traits in the semantic network facilitate faster or broader generalization?
*   **Implementation:** Weight similarity by Trait Centrality.
*   **Status:** *Drafting model logic*.

## 3. Decomposing $\lambda$ (Generalization Sharpness)
**The Question:** Does social identification increase the *rate* of projection or the *sharpness* of the category?
*   **Status:** *Pending final model parameters*.

## 4. Bayesian "Evidence Synthesis" (Pooled Hierarchical Analysis)
**The Question:** Is there a universal law governing self-to-group generalization across different contexts?
*   **Status:** *Script prepared (`run_pooled_analysis.R`)*. User running externally.

## 5. Latent Marginal Effects for Interpretation
**The Question:** What is the real-world impact of a 1-unit increase in social identification on social behavior?
*   **Status:** *Integrated in GD sections*.

## 6. MANDATORY: Parameter Recovery for New Models
**The Requirement:** Any model architecture that wins the comparison (e.g., Asymmetric + Lambda) must undergo a formal parameter recovery study.
*   **Procedure:**
    1.  Fit model to real data to get plausible group-level means/covariances.
    2.  Simulate new "synthetic" participants using those parameters.
    3.  Fit the model to synthetic data.
    4.  Correlation between "True" and "Recovered" parameters must be reported in Appendix.
*   **Focus:** Specifically verify $\alpha_{in}, \alpha_{out}$, and $\lambda$.

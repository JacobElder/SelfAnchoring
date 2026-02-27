# Self-Anchoring Project (JPSP Revisions)

## Project Overview
This project investigates the cognitive mechanisms of **self-anchoring**—the process by which individuals project their self-beliefs onto their ingroups. The primary goal is to refine the dissertation manuscript for JPSP, strengthening the theoretical ties between social identity theory and category learning/exemplar models (e.g., Generalized Context Model).

## Core Theoretical Framework
- **Similarity-Based Generalization:** Moves beyond one-to-one mapping (e.g., "I am outgoing, so my group is outgoing") to show that self-anchoring generalizes to semantically related but novel traits (e.g., "I am sociable, so my group is likely to be fun and witty").
- **Relational Similarity:** Leverages a pre-defined trait network (Elder et al., 2023) and Dice similarity to quantify the overlap of shared features between traits.
- **Contrastive Mechanisms:** Explores how intergroup dynamics (status, majority/minority size, and negation) augment self-projection rates.

## Study Summary
- **Study 1:** Minimal group paradigm (overestimators/underestimators). Established basic similarity-based generalization.
- **Study 2:** University status. Ingroup vs. Higher-Status Outgroup, Lower-Status Outgroup, or Negation ("Not Ingroup"). Found that higher-status outgroup contrasts can attenuate self-projection rates.
- **Study 3:** Racial group size. Asian/Latino Ingroups vs. Minority/Majority (White) Outgroups. Investigates distinctiveness and "common fate" as moderators of projection.

## Computational Modeling
- **Advanced Model Implementation:** `Computational Models/S_Logistic_Asym_Lambda.stan`
    - **Theoretical Rationale:** This model decouples the rate of **ingroup projection** ($\alpha_{in}$) and **outgroup repulsion** ($\alpha_{out}$), allowing for an empirical test of the "ingroup love vs. outgroup hate" hypothesis within a unified framework. It also introduces a **generalization sensitivity parameter** ($\lambda$) to model the rate of semantic decay, consistent with the Generalized Context Model.
    - **Key Parameters:**
        - `m_in`: Rate of projecting self-beliefs to the ingroup.
        - `m_out`: Rate of repelling self-beliefs from the outgroup.
        - `lambda`: Generalization sharpness/decay across semantic space.
        - `bias`: Ingroup/outgroup choice bias.
        - `tau`: Choice stochasticity (inverse temperature).
- **Previous Model:** `S_Logistic_1mOppose_Bias.stan`. Uses a single `m` parameter for both projection and rejection (symmetric assumption).

## Development Mandates
- **JPSP Focus:** Prioritize scientific merit, coherence, and the formalization of mechanisms that were previously only verbal descriptions in social identity literature.
- **Iterative Refinement:** Strategy should follow **Research -> Strategy -> Execution (Plan-Act-Validate)**.
- **Validation:** All code changes must be verified with existing or new R/Stan scripts to ensure convergence and parameter recoverability.

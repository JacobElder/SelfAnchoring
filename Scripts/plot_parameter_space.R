# Parameter Space Visualization
# Explores three implementations of how lambda and alpha shape social inference:
#   A. Generalization gradient curves (lambda effect on similarity weighting)
#   B. Projection gradient curves (alpha effect on self-rating -> group belief)
#   C. Venn diagram grid (self-group semantic overlap across parameter combinations)

library(tidyverse)
library(ggforce)
library(patchwork)
library(here)

# ============================================================
# PANEL A: Lambda — Generalization Breadth
# Shows S^lambda as a function of similarity for different lambda values.
# The key insight: lambda controls HOW FAR inference travels from the self.
# Low lambda = activation bleeds broadly across the semantic network.
# High lambda = activation is restricted to close semantic neighbors.
# ============================================================

sim_seq <- seq(0.001, 1, length.out = 300)

lambda_df <- expand.grid(S = sim_seq, lambda = c(0.5, 1, 2, 4)) %>%
  mutate(
    weight = S^lambda,
    lambda_label = factor(lambda,
      labels = c("λ = 0.5  (Analogical Thinker)",
                 "λ = 1.0  (Baseline)",
                 "λ = 2.0",
                 "λ = 4.0  (Compartmentalizer)"))
  )

p_lambda <- ggplot(lambda_df, aes(x = S, y = weight, color = lambda_label)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c("#1971C2", "#74C0FC", "#F08C00", "#E03131"),
                     name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\n(Unrelated)", ".25", ".50", ".75", "1\n(Identical)")) +
  labs(
    x = "Semantic Similarity to Self-Concept",
    y = expression("Generalization Weight  " ~ (S^lambda)),
    title = expression(bold("A.  ") * lambda * "  — Generalization Breadth")
  ) +
  annotate("text", x = 0.55, y = 0.92, hjust = 0,
           label = "Broad reach: distant traits still\ninfluence group beliefs",
           color = "#1971C2", size = 3.2) +
  annotate("text", x = 0.55, y = 0.12, hjust = 0,
           label = "Narrow reach: only close\nneighbors influence beliefs",
           color = "#E03131", size = 3.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.text = element_text(size = 9))

# ============================================================
# PANEL B: Alpha — Projection Force
# Shows G_in(E) = sigma(alpha * (E - 4)) for different alpha values.
# The key insight: alpha controls HOW STRONGLY self-rating drives group belief.
# Low alpha = self stays separate from group perception (objective observer).
# High alpha = self aggressively dictates group perception (egocentric anchor).
# ============================================================

E_seq <- seq(1, 7, length.out = 300)

alpha_df <- expand.grid(E = E_seq, alpha = c(0.5, 1.5, 4, 8)) %>%
  mutate(
    G_in  = plogis(alpha * (E - 4)),
    G_out = plogis(-alpha * (E - 4)),
    alpha_label = factor(alpha,
      labels = c("α = 0.5  (Objective Observer)",
                 "α = 1.5  (Moderate)",
                 "α = 4.0",
                 "α = 8.0  (Egocentric Anchor)"))
  )

p_alpha <- ggplot(alpha_df, aes(x = E, y = G_in, color = alpha_label)) +
  geom_line(linewidth = 1.3) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  scale_color_manual(values = c("#1971C2", "#74C0FC", "#F08C00", "#E03131"),
                     name = NULL) +
  scale_x_continuous(breaks = 1:7,
                     labels = c("1\nNot\nat all","2","3","4\nNeutral","5","6","7\nExtremely")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x = "Self-Rating on Trait",
    y = expression("Ingroup Projection Weight  " ~ G["in"](E)),
    title = expression(bold("B.  ") * alpha["in"] * "  — Projection Force (Ingroup)")
  ) +
  annotate("text", x = 6.7, y = 0.62, hjust = 1, size = 3.2, color = "#E03131",
           label = "High α: any self-descriptive\ntrait strongly projected") +
  annotate("text", x = 6.7, y = 0.72, hjust = 1, size = 3.2, color = "#1971C2",
           label = "Low α: self-description barely\ninfluences group perception") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.text = element_text(size = 9))

# ============================================================
# PANEL C: Venn Diagram Grid — Self-Group Semantic Overlap
# 2x2 grid: rows = Low/High lambda, cols = Low/High alpha_in
# Lambda controls the RADIUS of the self-concept circle (semantic reach).
# Alpha_in controls the OVERLAP between self and ingroup circles.
# ============================================================

make_venn_panel <- function(lambda_val, alpha_in_val, alpha_out_val = 3) {
  # Self radius: lower lambda = broader semantic reach = larger circle
  r_self <- 1.6 / sqrt(lambda_val)

  # Ingroup distance from self: higher alpha_in = stronger projection = more overlap
  # Map alpha_in (0-10) to distance: high alpha_in means centers are close
  d_in <- 2.2 - (alpha_in_val / 10) * 1.6

  # Outgroup distance: higher alpha_out = pushed further away
  d_out <- 2.2 + (alpha_out_val / 10) * 1.2

  data.frame(
    x0    = c(0,    d_in, -d_out),
    y0    = c(0,    0,     0),
    r     = c(r_self, 1.2,  1.2),
    label = c("Self", "Ingroup", "Outgroup"),
    fill  = c("#339AF0", "#FF6B6B", "#ADB5BD")
  )
}

venn_params <- list(
  list(lambda = 0.5, alpha_in = 1,  title = "Low λ, Low α\n(Broad reach, weak projection)"),
  list(lambda = 0.5, alpha_in = 8,  title = "Low λ, High α\n(Broad reach, strong projection)"),
  list(lambda = 3.5, alpha_in = 1,  title = "High λ, Low α\n(Narrow reach, weak projection)"),
  list(lambda = 3.5, alpha_in = 8,  title = "High λ, High α\n(Narrow reach, strong projection)")
)

venn_plots <- map(venn_params, function(p) {
  df <- make_venn_panel(p$lambda, p$alpha_in)

  ggplot(df) +
    geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = label),
                alpha = 0.35, color = "white", linewidth = 1.2) +
    geom_text(aes(x = x0, y = y0, label = label),
              size = 3.2, fontface = "bold", color = "gray20") +
    scale_fill_manual(values = c("Self" = "#339AF0",
                                 "Ingroup" = "#FF6B6B",
                                 "Outgroup" = "#ADB5BD")) +
    coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-3, 3)) +
    labs(title = p$title) +
    theme_void(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 9, color = "gray30"))
})

p_venn <- wrap_plots(venn_plots, nrow = 2) +
  plot_annotation(
    title = expression(bold("C.  ") * "Self-Group Semantic Overlap as a Function of  " * lambda * "  and  " * alpha["in"])
  )

# ============================================================
# COMBINE AND SAVE
# ============================================================

final_plot <- (p_lambda | p_alpha) / p_venn +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    caption = paste0(
      "Note. Panel A: Generalization weight S^λ as a function of semantic similarity for four λ values. ",
      "Panel B: Ingroup projection gradient G_in(E) as a function of self-rating for four α_in values. ",
      "Panel C: Conceptual Venn diagrams of self-ingroup-outgroup semantic overlap across λ × α_in combinations. ",
      "All other parameters held at representative values (τ = 3, γ = 0.5, w = 0.05)."
    ),
    theme = theme(plot.caption = element_text(size = 8, color = "gray40",
                                               hjust = 0, margin = margin(t = 8)))
  )

ggsave(here("Figures", "parameter_space_visualization.pdf"),
       final_plot, width = 12, height = 11)
ggsave(here("Figures", "parameter_space_visualization.png"),
       final_plot, width = 12, height = 11, dpi = 300)

message("Saved: Figures/parameter_space_visualization.pdf/.png")

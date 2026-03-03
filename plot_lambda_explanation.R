# Refined Visualization of how Lambda enables Alpha Identifiability (Staggered Labels)
library(tidyverse)

# Parameters
lambda_vals <- c(0.5, 1, 2, 4)
alpha <- 1.5 
E <- 7       
G <- 1 / (1 + exp(-alpha * (E - 4))) 

# Generate Similarity space
sim_space <- seq(0, 1, length.out = 100)

# Create curve data
df_curves <- expand.grid(Similarity = sim_space, lambda = lambda_vals) %>%
  mutate(
    Weight = (Similarity^lambda) * G,
    lambda_label = paste0("λ = ", lambda)
  )

# Define Traits to highlight
traits_to_show <- tibble(
  Trait = c("Rated Trait", "Novel Trait"),
  Similarity = c(1.0, 0.6),
  Rating = c("7", "?")
)

# Calculate points on the curves for these traits
trait_points <- expand.grid(Similarity = traits_to_show$Similarity, lambda = lambda_vals) %>%
  left_join(traits_to_show, by = "Similarity") %>%
  mutate(
    Weight = (Similarity^lambda) * G
  )

# STAGGERED Label positions for lambda
# We pick different X positions where the curves are naturally separated
label_pos <- tibble(
  lambda = c(0.5, 1, 2, 4),
  x_pos = c(0.1, 0.2, 0.35, 0.5)
) %>%
  mutate(
    y_pos = (x_pos^lambda) * G + 0.08
  )

# Plot
p <- ggplot() +
  # Curves
  geom_line(data = df_curves, aes(x = Similarity, y = Weight, group = as.factor(lambda)), 
            color = "black", linewidth = 0.8) +
  # Staggered Lambda labels in white boxes
  geom_label(data = label_pos, aes(x = x_pos, y = y_pos, label = paste0("λ = ", lambda)),
            size = 3.5, fontface = "bold", family = "serif", fill = "white", linewidth = 0.2) +
  # Points for the rated trait
  geom_point(data = trait_points %>% filter(Trait == "Rated Trait"), 
             aes(x = Similarity, y = Weight), size = 3, color = "black") +
  # Points for the novel trait
  geom_point(data = trait_points %>% filter(Trait == "Novel Trait"), 
             aes(x = Similarity, y = Weight), size = 3, color = "black", shape = 21, fill = "white") +
  # Vertical line connecting the novel trait points
  geom_segment(aes(x = 0.6, xend = 0.6, y = (0.6^4)*G, yend = (0.6^0.5)*G), 
               linetype = "dashed", alpha = 0.3) +
  # Boxed annotations
  annotate("label", x = 1.0, y = G + 0.08, label = "Self-Rating: 7", 
           fontface = "italic", family = "serif", fill = "white", linewidth = 0.5) +
  annotate("label", x = 0.6, y = (0.6^0.5)*G + 0.08, label = "Novel Trait", 
           fontface = "italic", family = "serif", fill = "white", linewidth = 0.5) +
  annotate("label", x = 0.75, y = 0.2, label = "Higher λ flattens\ninfluence of novel traits", 
           hjust = 0, size = 3.2, family = "serif", fill = "#f8f8f8", linewidth = 0.2) +
  labs(
    x = "Semantic Similarity to Rated Trait (S)",
    y = "Generalization Strength (S^λ * Belief)"
  ) +
  scale_x_continuous(expand = c(0, 0.05), limits = c(0, 1.1)) +
  scale_y_continuous(expand = c(0, 0.05), limits = c(0, 1.1)) +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 12, family = "serif"),
    axis.text = element_text(size = 10, family = "serif"),
    legend.position = "none"
  )

# Save
ggsave("Figure_Lambda_Identifiability.png", p, width = 7, height = 5, dpi = 300)

message("Figure_Lambda_Identifiability.png has been generated (staggered labels).")

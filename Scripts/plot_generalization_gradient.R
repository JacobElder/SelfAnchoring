# Visualization of Social Generalization Gradient (Shepard's Law)
library(tidyverse)

# Study 1 Parameters (Asymmetric + Lambda)
m_in <- 0.0755
m_out <- 0.640
lambda <- 2.64

# Define the Latent Belief function (inv_logit)
# Assuming a maximal self-evaluation of 7 (7-4 = 3) to show the full gradient
G_in_max <- 1 / (1 + exp(-m_in * 3))
G_out_max <- 1 / (1 + exp(m_out * 3)) # Note: -m_out * (E-4) in Stan code

# Generate Similarity values from 0 to 1
sim_space <- seq(0, 1, length.out = 100)

# Calculate Generalization Strength (S^lambda * G)
df <- tibble(
  Similarity = sim_space,
  Distance = 1 - sim_space,
  Ingroup_Assimilation = (Similarity^lambda) * G_in_max,
  Outgroup_Repulsion = (Similarity^lambda) * G_out_max
)

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = c(Ingroup_Assimilation, Outgroup_Repulsion),
               names_to = "Mechanism",
               values_to = "Generalization_Weight")

# Plot
p <- ggplot(df_long, aes(x = Similarity, y = Generalization_Weight, color = Mechanism)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c("Ingroup_Assimilation" = "#2c7bb6", "Outgroup_Repulsion" = "#d7191c")) +
  labs(
    title = "Universal Law of Social Generalization (Study 1)",
    subtitle = paste0("Generalization Sensitivity (lambda) = ", round(lambda, 2)),
    x = "Semantic Similarity (S)",
    y = "Generalization Strength (S^lambda * Latent Belief)",
    caption = "Curves represent the decay of self-to-group inference across semantic space for a self-evaluation of '7'."
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Save the plot
ggsave("Figure_Generalization_Gradient_S1.png", p, width = 7, height = 5, dpi = 300)

message("Figure_Generalization_Gradient_S1.png has been generated.")

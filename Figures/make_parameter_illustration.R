# Concrete Visualization: From Self-Ratings to Ingroup Classifications
# Shows how the same self-rating profile produces different predictions
# for two contrasting parameter profiles (low vs high alpha, narrow vs broad lambda).
#
# Output: Figures/parameter_illustration.pdf (and .png)
# Panel A: Lollipop of self-rating profile for 8 training traits
# Panel B: Similarity structure (heatmap of Dice similarities)
# Panel C: Predicted P(ingroup) for trained + novel traits — low-alpha participant
# Panel D: Predicted P(ingroup) for trained + novel traits — high-alpha participant

library(tidyverse)
library(patchwork)
library(here)

set.seed(42)

# ── 1. Define training traits and self-ratings ────────────────────────────────
# Self-ratings are on a discrete 1–7 Likert scale (integer values only)
training_traits <- c("Sociable","Warm","Outgoing","Friendly",
                     "Organized","Accurate","Quiet","Reserved")
self_ratings    <- c(7L, 6L, 7L, 6L, 4L, 4L, 2L, 2L)   # 1–7 integer Likert scale

novel_traits    <- c("Fun","Witty","Helpful","Calm","Careful","Shy")

all_traits      <- c(training_traits, novel_traits)
n_train         <- length(training_traits)
n_all           <- length(all_traits)

# ── 2. Synthetic Dice similarity matrix ──────────────────────────────────────
# Social cluster (Sociable, Warm, Outgoing, Friendly, Fun, Witty, Helpful) vs
# Quiet cluster (Organized, Accurate, Quiet, Reserved, Calm, Careful, Shy)
social_idx <- c(1,2,3,4, 9,10,11)     # indices in all_traits
quiet_idx  <- c(5,6,7,8, 12,13,14)

sim_mat <- matrix(0.05, nrow=n_all, ncol=n_all)
diag(sim_mat) <- 1.0
# Within-cluster similarities
for (i in social_idx) for (j in social_idx) if (i != j) sim_mat[i,j] <- runif(1, 0.40, 0.70)
for (i in quiet_idx)  for (j in quiet_idx)  if (i != j) sim_mat[i,j] <- runif(1, 0.35, 0.65)
# Cross-cluster (low)
for (i in social_idx) for (j in quiet_idx) {
  sim_mat[i,j] <- sim_mat[j,i] <- runif(1, 0.02, 0.12)
}
# Ensure symmetry and clip
sim_mat <- (sim_mat + t(sim_mat)) / 2
diag(sim_mat) <- 1.0
sim_mat <- pmin(pmax(sim_mat, 0), 1)
rownames(sim_mat) <- colnames(sim_mat) <- all_traits

# ── 3. Model prediction function ──────────────────────────────────────────────
predict_ingroup <- function(test_trait_idx, sim_mat, self_ratings,
                            alpha, lambda, bias, lapse = 0.3) {
  # GP from self-ratings on training traits
  GP_in  <- plogis( alpha * (self_ratings - 4))    # ingroup belief
  GP_out <- plogis(-alpha * (self_ratings - 4))    # outgroup belief

  # Similarity-weighted sums (Shepard's law applied to training→test Dice)
  sims <- sim_mat[1:length(training_traits), test_trait_idx] ^ lambda

  simW_in  <- sum(sims * GP_in)
  simW_out <- sum(sims * GP_out)

  # Choice probability
  if (simW_in + simW_out < 1e-9) return(bias)
  p_model <- simW_in / (simW_in + simW_out)
  # Bias shift (logit scale)
  p_biased <- plogis(qlogis(p_model) + qlogis(bias))
  # Lapse
  p_final <- lapse * 0.5 + (1 - lapse) * p_biased
  return(p_final)
}

# ── 4. Compute predictions for two parameter profiles ────────────────────────
profiles <- list(
  "Low projection (α=0.5, λ=1.0)" = list(alpha=0.5, lambda=1.0, bias=0.5, lapse=0.3, col="#4E9A9A"),
  "High projection (α=3.0, λ=3.5)" = list(alpha=3.0, lambda=3.5, bias=0.5, lapse=0.3, col="#C7522A")
)

pred_df <- map_dfr(names(profiles), function(prof_name) {
  p <- profiles[[prof_name]]
  map_dfr(seq_len(n_all), function(tidx) {
    tibble(
      trait        = all_traits[tidx],
      trait_type   = if (tidx <= n_train) "Training" else "Novel",
      p_ingroup    = predict_ingroup(tidx, sim_mat, self_ratings,
                                     p$alpha, p$lambda, p$bias, p$lapse),
      profile      = prof_name
    )
  })
})

# Preserve trait order
pred_df$trait <- factor(pred_df$trait, levels=all_traits)

# ── 5. Panel A: Self-rating lollipop ─────────────────────────────────────────
rating_df <- tibble(
  trait  = factor(training_traits, levels=rev(training_traits)),
  rating = self_ratings,
  cluster = c(rep("Social",4), rep("Quiet/Neutral",4))
)

p_A <- ggplot(rating_df, aes(x=trait, y=rating, color=cluster)) +
  geom_segment(aes(xend=trait, y=1, yend=rating), linewidth=0.8) +
  geom_point(size=4) +
  coord_flip() +
  scale_y_continuous(limits=c(1,7), breaks=1:7, name="Self-Rating (1–7)") +
  scale_color_manual(values=c("Social"="#4E9A9A","Quiet/Neutral"="#8C8C8C"),
                     name="Trait cluster") +
  labs(title="(A) Training trait\nself-ratings", x=NULL) +
  geom_hline(yintercept=4, linetype="dashed", color="grey60", linewidth=0.5) +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom",
        plot.title=element_text(face="bold", size=11),
        panel.grid.minor=element_blank())

# ── 6. Panel B: Similarity heatmap (training traits only) ────────────────────
sim_long <- as.data.frame(sim_mat[1:n_train, 1:n_train]) |>
  rownames_to_column("from") |>
  pivot_longer(-from, names_to="to", values_to="similarity") |>
  mutate(from = factor(from, levels=training_traits),
         to   = factor(to,   levels=training_traits))

p_B <- ggplot(sim_long, aes(x=to, y=from, fill=similarity)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient(low="white", high="#2B5C8A", limits=c(0,1),
                      name="Dice\nsimilarity") +
  labs(title="(B) Semantic similarity\namong training traits", x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        axis.text.y=element_text(size=9),
        plot.title=element_text(face="bold", size=11),
        legend.position="right")

# ── 7. Panels C & D: Predicted P(ingroup) for each profile ───────────────────
col_map <- setNames(
  c(profiles[[1]]$col, profiles[[2]]$col),
  names(profiles)
)

make_pred_panel <- function(prof_name, panel_label) {
  sub_df <- pred_df |> filter(profile == prof_name) |>
    mutate(trait = factor(trait, levels=rev(all_traits)))
  ggplot(sub_df, aes(x=trait, y=p_ingroup, fill=trait_type)) +
    geom_col(width=0.6, alpha=0.85) +
    geom_hline(yintercept=0.5, linetype="dashed", color="grey40", linewidth=0.5) +
    coord_flip() +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.25),
                       labels=scales::percent_format(1),
                       name="P(Ingroup classification)") +
    scale_fill_manual(values=c("Training"="#2B5C8A","Novel"="#A8C8E8"),
                      name="Trait type") +
    labs(title=paste0("(", panel_label, ") ", prof_name), x=NULL) +
    theme_minimal(base_size=11) +
    theme(legend.position="bottom",
          plot.title=element_text(face="bold", size=10),
          panel.grid.minor=element_blank())
}

p_C <- make_pred_panel(names(profiles)[1], "C")
p_D <- make_pred_panel(names(profiles)[2], "D")

# ── 8. Compose layout ─────────────────────────────────────────────────────────
# No embedded title/subtitle: per APA 7 convention, figure title goes in the
# figure caption in the manuscript, not inside the image.
fig <- (p_A | p_B) / (p_C | p_D)

# ── 9. Save ───────────────────────────────────────────────────────────────────
ggsave(here("Figures","parameter_illustration.pdf"), fig, width=13, height=10)
ggsave(here("Figures","parameter_illustration.png"), fig, width=13, height=10, dpi=200)
ggsave(here("Figures","fig04_parameter_illustration.tiff"), fig, width=13, height=10, dpi=300)
message("Saved: Figures/parameter_illustration.pdf/.png + fig04_parameter_illustration.tiff")

# ── 10. Also print summary table for quick inspection ─────────────────────────
cat("\n=== Predicted P(ingroup) by profile and trait ===\n")
pred_df |>
  select(profile, trait, trait_type, p_ingroup) |>
  mutate(p_ingroup = round(p_ingroup, 3)) |>
  arrange(profile, trait) |>
  print(n=Inf)

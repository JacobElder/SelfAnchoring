# make_trait_segregation_fig.R
# Conceptual diagram: High vs. Low Trait Segregation in a semantic trait network.
# Uses ACTUAL traits and ACTUAL edges from Pooled/input/adjacencyMatrix_p.csv.
# Two-panel figure: same network topology, same node positions, same trait labels.
# Only node coloring differs — high segregation clusters by semantic neighborhood,
# low segregation scatters ingroup/outgroup labels across semantic structure.
#
# Selected clusters (all verified in allPosCents.csv):
#   Knowledgeable: Precise, Smart, Rational, Capable, Clever
#   Respectful:    Friendly, Respectful, Peaceful, Positive, Good-natured
#
# Output: Figures/fig_trait_segregation.tiff

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})

TIFF_DPI   <- 300
TIFF_UNITS <- "in"

# ── Load network ──────────────────────────────────────────────────────────────
traits_df <- read.csv(here("Pooled/input/allPosCents.csv"))
posDf     <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))

# ── Selected traits ───────────────────────────────────────────────────────────
sel_k <- c("Precise", "Smart", "Rational", "Capable", "Clever")
sel_r <- c("Friendly", "Respectful", "Peaceful", "Positive", "Good-natured")
sel   <- c(sel_k, sel_r)

idx     <- match(sel, traits_df$trait)
sub_adj <- as.matrix(posDf)[idx, idx]
rownames(sub_adj) <- colnames(sub_adj) <- sel

# ── Build edge list from actual adjacency matrix ──────────────────────────────
edge_pairs <- which(sub_adj == 1 & upper.tri(sub_adj), arr.ind = TRUE)
edges_real <- data.frame(from = edge_pairs[, 1], to = edge_pairs[, 2])

# ── Fixed node positions (manual, ensures clear cluster separation) ────────────
# Knowledgeable cluster left; Respectful cluster right.
# Order: Precise, Smart, Rational, Capable, Clever,
#        Friendly, Respectful, Peaceful, Positive, Good-natured
positions <- data.frame(
  x = c(-2.8, -3.8, -2.2, -3.8, -1.6,   # Knowledgeable (left)
          2.8,  1.8,  3.8,  2.2,  3.8),   # Respectful (right)
  y = c( 0.2,  1.1, -1.1, -1.0,  0.9,   # Knowledgeable
         0.2,  1.2,  1.2, -1.0, -0.5)    # Respectful
)

# ── Node color assignments ─────────────────────────────────────────────────────
# HIGH: Knowledgeable cluster = Ingroup (blue), Respectful = Outgroup (pink)
colors_high <- c(rep("Ingroup", 5), rep("Outgroup", 5))

# LOW: labels scattered across semantic neighborhoods
colors_low  <- c("Ingroup",  "Outgroup", "Outgroup", "Ingroup",  "Outgroup",
                 "Outgroup", "Ingroup",  "Ingroup",  "Outgroup", "Ingroup")

# ── Panel-building function ────────────────────────────────────────────────────
node_palette <- c("Ingroup" = "#0072B2", "Outgroup" = "#CC79A7")

make_panel <- function(node_colors, title_txt) {
  node_df <- data.frame(
    id    = seq_along(sel),
    label = sel,
    group = node_colors,
    x     = positions$x,
    y     = positions$y
  )
  node_df$label[node_df$label == "Good-natured"] <- "Good-\nnatured"

  edge_df <- edges_real |>
    mutate(
      x    = positions$x[from],
      y    = positions$y[from],
      xend = positions$x[to],
      yend = positions$y[to]
    )

  ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey55", linewidth = 0.75, alpha = 0.75) +
    geom_point(data = node_df,
               aes(x = x, y = y, fill = group),
               size = 15, shape = 21, stroke = 0.5, color = "white") +
    scale_fill_manual(values = node_palette, name = NULL) +
    geom_text(data = node_df,
              aes(x = x, y = y, label = label),
              size = 2.7, fontface = "bold", color = "black",
              lineheight = 0.85) +
    labs(title = title_txt) +
    coord_cartesian(xlim = c(-5.2, 5.2), ylim = c(-2.8, 2.8), clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5,
                                     margin = margin(b = 6)),
      legend.position = "bottom",
      legend.text     = element_text(size = 10.5),
      legend.key.size = unit(1.1, "lines"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(8, 12, 4, 12)
    )
}

p_high <- make_panel(colors_high, "High Trait Segregation")
p_low  <- make_panel(colors_low,  "Low Trait Segregation")

fig_ts <- (p_high | p_low) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("Figures", "fig_trait_segregation.tiff"),
       fig_ts, width = 10, height = 5, dpi = TIFF_DPI, units = TIFF_UNITS,
       compression = "lzw")

message("Saved: Figures/fig_trait_segregation.tiff")
message(sprintf("Traits: %s", paste(sel, collapse = ", ")))
message(sprintf("Edges: %d (all from actual adjacency matrix)", nrow(edges_real)))

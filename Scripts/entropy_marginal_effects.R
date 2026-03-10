# entropy_marginal_effects.R
# Self-concept entropy (H_g) as a predictor of ingroup classification.
#
# H_g = similarity-weighted Shannon entropy of self-evaluation distribution
#       over the semantic neighborhood of each test trait g.
# High H_g: mixed/ambiguous self-ratings near g → noisy ingroup signal.
# Low H_g:  concentrated self-ratings near g → coherent ingroup signal.
#
# Theoretical motivation: Hogg's self-uncertainty theory (2000, 2007).
# Prediction: higher H_g → lower P(ingroup), because the self-concept
# signal in that semantic neighborhood is too diffuse to drive projection.
#
# OUTPUT:
#   Results/entropy_ames.csv            — AMEs per study
#   Results/entropy_predictions.csv     — dense-grid predictions for plotting
#   Figures/fig_entropy_effects.tiff    — 3-panel marginal effects plot

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(marginaleffects)
  library(igraph)
  library(here)
  library(patchwork)
})

TIFF_DPI   <- 300
TIFF_UNITS <- "in"

# ── Load trait network once ───────────────────────────────────────────────────
message("Loading trait network...")
posDf    <- read.csv(here("Pooled/input/adjacencyMatrix_p.csv"))
posGraph <- graph_from_adjacency_matrix(as.matrix(posDf), mode = "max")
simMat   <- similarity(posGraph, method = "dice")

# ── Core computation: H_g per (subject, test trait) ──────────────────────────
# H_g = -Σ_e p_e * log2(p_e)
# p_e = Σ_{train t: E_t = e} sim(g, t) / Σ_t sim(g, t)
compute_entropy_df <- function(fullTest, traindf) {
  uIds  <- sort(intersect(unique(fullTest$subID), unique(traindf$subID)))
  rows  <- vector("list", length(uIds))

  for (i in seq_along(uIds)) {
    s_test  <- filter(fullTest, subID == uIds[i])
    s_train <- filter(traindf,  subID == uIds[i])
    if (nrow(s_test) == 0 || nrow(s_train) == 0) next

    H <- numeric(nrow(s_test))
    for (k in seq_len(nrow(s_test))) {
      g    <- s_test$Idx[k]
      sims <- simMat[g, s_train$Idx]
      tot  <- sum(sims, na.rm = TRUE)
      if (is.na(tot) || tot < 1e-9) { H[k] <- NA; next }
      p_e  <- sapply(1:7, function(e) {
        idx_e <- which(s_train$selfResp == e)
        if (length(idx_e) == 0) return(0)
        sum(sims[idx_e], na.rm = TRUE) / tot
      })
      p_nz <- p_e[p_e > 0]
      H[k] <- -sum(p_nz * log2(p_nz))
    }
    rows[[i]] <- s_test |>
      mutate(entropy = H, subID = uIds[i])
  }
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ── Load data per study ───────────────────────────────────────────────────────
load_study <- function(test_path, train_path) {
  tp <- here(test_path); trp <- here(train_path)
  if (!file.exists(tp))  { message("Missing: ", test_path);  return(NULL) }
  if (!file.exists(trp)) { message("Missing: ", train_path); return(NULL) }
  list(
    test  = read.csv(tp)  |> filter(!is.na(ingChoiceN)),
    train = read.csv(trp) |> filter(!is.na(selfResp))
  )
}

s1 <- load_study("Study 1/Cleaning/output/fullTest.csv",
                 "Study 1/Cleaning/output/fullTrain.csv")
s2 <- load_study("Study 2/Cleaning/output/fullTest.csv",
                 "Study 2/Cleaning/output/fullTrain_fixed.csv")
s3 <- load_study("Study 3/Cleaning/output/fullTest_fixed.csv",
                 "Study 3/Cleaning/output/fullTrain_fixed.csv")

# ── Fit model and compute AMEs per study ─────────────────────────────────────
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

run_entropy_study <- function(study_data, study_label) {
  if (is.null(study_data)) return(NULL)

  message(sprintf("\n── %s: computing entropy ──", study_label))
  df <- compute_entropy_df(study_data$test, study_data$train)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Winsorise entropy at 1st/99th percentile, then Z-score
  q <- quantile(df$entropy, c(0.01, 0.99), na.rm = TRUE)
  df <- df |>
    filter(!is.na(entropy), entropy >= q[1], entropy <= q[2]) |>
    mutate(
      entropy.Z    = scale(entropy)[, 1],
      ingChoiceN   = as.integer(ingChoiceN),
      desirability.Z = scale(desirability)[, 1]
    )

  message(sprintf("  %d obs, entropy mean=%.2f SD=%.2f",
                  nrow(df), mean(df$entropy), sd(df$entropy)))

  # M_entropy: H_g + desirability control (same structure as M2/M3 in main scripts)
  m <- tryCatch(
    glmer(
      ingChoiceN ~ entropy.Z + desirability.Z +
        (entropy.Z | subID) + (1 | trait),
      data = df, family = binomial, control = ctrl, nAGQ = 1
    ),
    error = function(e) {
      message("  Random slopes failed, retrying with intercept-only: ", e$message)
      tryCatch(
        glmer(ingChoiceN ~ entropy.Z + desirability.Z +
                (1 | subID) + (1 | trait),
              data = df, family = binomial, control = ctrl, nAGQ = 1),
        error = function(e2) { message("  glmer failed: ", e2$message); NULL }
      )
    }
  )
  if (is.null(m)) return(NULL)

  beta <- fixef(m)["entropy.Z"]
  message(sprintf("  beta(entropy.Z) = %.3f", beta))

  # AME
  ame <- tryCatch(
    avg_comparisons(m, variables = "entropy.Z"),
    error = function(e) { message("  AME failed: ", e$message); NULL }
  )

  # Dense-grid predictions for plotting
  x_seq <- seq(min(df$entropy.Z, na.rm = TRUE),
               max(df$entropy.Z, na.rm = TRUE), length.out = 50)
  preds_raw <- tryCatch(
    avg_predictions(m, variables = list(entropy.Z = x_seq)),
    error = function(e) { message("  avg_predictions failed: ", e$message); NULL }
  )

  list(
    study  = study_label,
    data   = df,
    model  = m,
    ame    = ame,
    preds  = preds_raw
  )
}

res_s1 <- run_entropy_study(s1, "Study 1 (Minimal Groups)")
res_s2 <- run_entropy_study(s2, "Study 2 (University Status)")
res_s3 <- run_entropy_study(s3, "Study 3 (Racial Groups)")

# ── Save AMEs ─────────────────────────────────────────────────────────────────
ame_list <- Filter(Negate(is.null), lapply(
  list(res_s1, res_s2, res_s3),
  function(r) {
    if (is.null(r) || is.null(r$ame)) return(NULL)
    as.data.frame(r$ame) |> mutate(study = r$study)
  }
))

if (length(ame_list) > 0) {
  ame_df <- do.call(rbind, ame_list)
  write.csv(ame_df, here("Results", "entropy_ames.csv"), row.names = FALSE)
  message("\nSaved: Results/entropy_ames.csv")
  print(ame_df[, intersect(c("study","term","estimate","std.error","conf.low","conf.high","p.value"),
                            names(ame_df))], digits = 3, row.names = FALSE)
}

# ── Save predictions ──────────────────────────────────────────────────────────
pred_list <- Filter(Negate(is.null), lapply(
  list(res_s1, res_s2, res_s3),
  function(r) {
    if (is.null(r) || is.null(r$preds)) return(NULL)
    as.data.frame(r$preds) |>
      rename(predictor_val = entropy.Z) |>
      mutate(study = r$study, predictor = "entropy")
  }
))

if (length(pred_list) > 0) {
  pred_df <- do.call(rbind, pred_list)
  write.csv(pred_df, here("Results", "entropy_predictions.csv"), row.names = FALSE)
  message("Saved: Results/entropy_predictions.csv")
}

# ── Build figure ──────────────────────────────────────────────────────────────
theme_dissert <- function(base_size = 11) {
  theme(
    panel.background  = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_blank(),
    axis.ticks        = element_line(colour = "black", linewidth = 0.4),
    axis.text         = element_text(size = base_size - 1, colour = "black"),
    axis.title        = element_text(size = base_size, face = "bold", colour = "black"),
    plot.title        = element_text(size = base_size, face = "bold"),
    plot.background   = element_rect(fill = "white", colour = NA)
  )
}

make_entropy_panel <- function(res, title_txt, color = "black") {
  if (is.null(res) || is.null(res$preds)) {
    return(
      ggplot() + theme_void() +
        annotate("rect", xmin=0, xmax=1, ymin=0, ymax=1,
                 fill="grey95", color="grey70", linewidth=0.5) +
        annotate("text", x=0.5, y=0.5, hjust=0.5, size=3, colour="grey50",
                 label=paste0(gsub("\n"," ",title_txt), "\n[Pending]")) +
        theme(plot.background = element_rect(fill="white", colour=NA))
    )
  }

  pd <- as.data.frame(res$preds) |> rename(predictor_val = entropy.Z)

  # AME label
  ame_lab <- if (!is.null(res$ame)) {
    a <- as.data.frame(res$ame)
    sprintf("AME = %.3f\n95%% CI [%.3f, %.3f]\np = %.3f",
            a$estimate, a$conf.low, a$conf.high, a$p.value)
  } else { "" }

  ggplot(pd, aes(x = predictor_val, y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                fill = color, alpha = 0.15, colour = NA) +
    geom_line(colour = color, linewidth = 1.0) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey55", linewidth = 0.45) +
    annotate("label", x = Inf, y = Inf,
             label = ame_lab, hjust = 1.05, vjust = 1.2,
             size = 2.8, colour = "grey20", fill = "white", linewidth = 0.3) +
    scale_y_continuous(labels = scales::percent_format(1),
                       limits = c(NA, NA),
                       expand = expansion(mult = 0.05)) +
    labs(
      title = title_txt,
      x     = "Self-Concept Entropy, H\u2092 (Z)",
      y     = "P(Ingroup Classification)"
    ) +
    theme_dissert()
}

p1 <- make_entropy_panel(res_s1, "Study 1\n(Minimal Groups)",    "#2B5C8A")
p2 <- make_entropy_panel(res_s2, "Study 2\n(University Status)", "#E69F00")
p3 <- make_entropy_panel(res_s3, "Study 3\n(Racial Groups)",     "#CC79A7")

fig_ent <- (p1 | p2 | p3) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

ggsave(here("Figures", "fig_entropy_effects.tiff"),
       fig_ent, width = 11, height = 4.5,
       dpi = TIFF_DPI, units = TIFF_UNITS, compression = "lzw")

message("Saved: Figures/fig_entropy_effects.tiff")
message("\nDone. Bayesian β values for reference:")
message("  S1: β = -.739 (pd = 98.7%, p = .026)")
message("  S2: β = -.182 (p = .122, not significant)")
message("  S3: β = -.187 (p = .242, not significant)")
message("The frequentist AMEs above should match directionally.")

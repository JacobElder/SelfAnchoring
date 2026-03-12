# app.R
# Self-Anchoring Generalization Model — Interactive Supplement
# Elder et al. (2026), JPSP
#
# Training traits are randomly sampled with self-ratings moderately correlated
# with semantic position (~r = .70). Click "Resample" to draw a new
# configuration and see how the prediction curve changes.
#
# Deploy: rsconnect::deployApp("Shiny/")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(scales)
})

# ── Helpers ───────────────────────────────────────────────────────────────────

inv_logit <- function(x) 1 / (1 + exp(-x))
logit     <- function(p) log(p / (1 - p))

# Sample 20 training traits.
# pos  ~ Uniform(0, 1): semantic position (0 = distant from prototype,
#                                          1 = close to prototype)
# self_rating ~ correlated with pos at r_target, clamped to {1,...,7}
# r_target = 0 → ratings are uniformly random (no self-anchoring structure)
# r_target = 1 → ratings increase perfectly with position
sample_traits <- function(n = 90, r_target = 0.70) {
  pos        <- runif(n)
  signal_var <- 25 / 12          # Var(pos * 5) for pos ~ U(0,1)
  if (r_target < 0.02) {
    sr <- sample(1L:7L, n, replace = TRUE)
  } else {
    noise_sd <- sqrt(signal_var * (1 / pmin(r_target, 0.999)^2 - 1))
    raw      <- pos * 5 + rnorm(n, sd = noise_sd)
    sr       <- pmax(1L, pmin(7L, round(raw + 1.5)))
  }
  data.frame(pos = pos, self_rating = sr)
}

# Symmetric + λ: one α converts self-ratings into group beliefs
compute_sym <- function(x, traits, alpha, lambda, bias, w) {
  sim_raw  <- pmax(0, 1 - abs(x - traits$pos))
  G        <- inv_logit(alpha * (traits$self_rating - 4))
  simW_in  <- sum(sim_raw^lambda * G)       + 1e-9
  simW_out <- sum(sim_raw^lambda * (1 - G)) + 1e-9
  p_model  <- inv_logit(logit(bias) + log(simW_in) - log(simW_out))
  w * 0.5 + (1 - w) * p_model
}

# Asymmetric + λ: separate α_in (ingroup assimilation) and α_out (outgroup repulsion)
# component: "combined" | "ingroup_only" (α_out→0) | "outgroup_only" (α_in→0)
compute_asym <- function(x, traits, alpha_in, alpha_out, lambda, bias, w,
                          component = "combined") {
  sim_raw <- pmax(0, 1 - abs(x - traits$pos))
  Gin  <- if (component == "outgroup_only") rep(0.5, nrow(traits)) else
            inv_logit( alpha_in  * (traits$self_rating - 4))
  Gout <- if (component == "ingroup_only")  rep(0.5, nrow(traits)) else
            inv_logit(-alpha_out * (traits$self_rating - 4))
  simW_in  <- sum(sim_raw^lambda * Gin)  + 1e-9
  simW_out <- sum(sim_raw^lambda * Gout) + 1e-9
  p_model  <- inv_logit(logit(bias) + log(simW_in) - log(simW_out))
  w * 0.5 + (1 - w) * p_model
}

# ── Shared theme & constants ──────────────────────────────────────────────────

theme_app <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93"),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(size = 11, colour = "grey30")
    )
}

S_seq        <- seq(0.01, 0.99, length.out = 300)
COL_SYM      <- "#2B5C8A"
COL_ASYM_IN  <- "#2B5C8A"
COL_ASYM_OUT <- "#CC79A7"

# ── Scatter helper (shared by both tabs) ─────────────────────────────────────

make_scatter <- function(traits, accent = "#2B5C8A") {
  ggplot(traits, aes(x = pos, y = self_rating)) +
    geom_smooth(method = "lm", se = FALSE,
                colour = "grey70", linewidth = 0.5, linetype = "dashed",
                formula = y ~ x) +
    geom_point(aes(fill = self_rating), shape = 21, size = 3.8,
               colour = "white", stroke = 0.3) +
    scale_fill_gradient(low = "#ccdff0", high = accent,
                        limits = c(1, 7), guide = "none") +
    scale_x_continuous(
      "Semantic Similarity",
      limits = c(0, 1), breaks = c(0, 0.5, 1),
      labels = c("0\n(distant)", "0.5", "1\n(close)")
    ) +
    scale_y_continuous(
      "Self-rating", limits = c(0.5, 7.5), breaks = c(1, 4, 7),
      labels = c("1\n(Not at all)", "4", "7\n(Extremely)")
    ) +
    labs(title = "Training traits  (N = 90)") +
    theme_app()
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_navbar(
  title = span("Self-Anchoring Generalization Model", style = "font-weight:600;"),
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2B5C8A",
    base_font  = font_google("Source Sans Pro")
  ),
  bg = "#2B5C8A", inverse = TRUE,

  # ── Tab 1: Symmetric + λ ─────────────────────────────────────────────────
  nav_panel(
    title = "Symmetric + \u03bb",
    layout_sidebar(
      sidebar = sidebar(
        width = 270, open = TRUE,
        h6("Model Parameters", class = "text-primary fw-bold mt-1"),
        sliderInput("sym_alpha",
                    HTML("Projection rate (\u03b1)"),
                    min = 0.1, max = 10, value = 2.8, step = 0.1),
        sliderInput("sym_lambda",
                    HTML("Generalization sensitivity (\u03bb)"),
                    min = 0.1, max = 5, value = 3.6, step = 0.1),
        sliderInput("sym_bias",
                    HTML("Ingroup bias (\u03b3)"),
                    min = 0.01, max = 0.99, value = 0.36, step = 0.01),
        sliderInput("sym_w", "Lapse rate (w)",
                    min = 0, max = 0.99, value = 0.56, step = 0.01),
        hr(),
        actionButton("sym_reset",    "Reset to Study 1 estimates",
                     class = "btn-sm btn-outline-primary w-100 mb-2"),
        hr(),
        h6("Training Trait Structure", class = "text-primary fw-bold"),
        sliderInput("sym_r",
                    HTML("Self-rating coherence (<em>r</em>)"),
                    min = 0, max = 0.95, value = 0.70, step = 0.05),
        actionButton("sym_resample", "Resample training traits",
                     class = "btn-sm btn-outline-secondary w-100 mt-1")
      ),
      layout_columns(
        col_widths = c(5, 7),
        card(
          card_header("Training Traits"),
          plotOutput("sym_scatter", height = "390px"),
          card_footer(HTML(
            "<small>Each dot is a training trait. Dots are colored by self-rating.
             The dashed line shows the overall trend. Use the
             <em>Self-rating coherence (r)</em> slider to control how structured
             the self-concept is; click <em>Resample</em> to draw another set
             at the same coherence level.</small>"
          ))
        ),
        card(
          full_screen = TRUE,
          card_header("P(Ingroup Classification) \u2014 Symmetric + \u03bb"),
          plotOutput("sym_plot", height = "390px"),
          card_footer(HTML(
            "<small>
             <b>X-axis:</b> semantic position of the test trait relative to the
             ingroup prototype (0 = maximally distant, 1 = maximally close).
             This reflects the self-evaluation-weighted similarity mean:
             P(Ingroup) rises as the test trait falls closer to training traits
             on which the participant rated themselves highly.
             <b>\u03bb</b> sharpens the gradient; <b>\u03b1</b> converts self-ratings
             into group beliefs; <b>w</b> compresses toward chance (0.50).
             </small>"
          ))
        )
      )
    )
  ),

  # ── Tab 2: Asymmetric + λ ────────────────────────────────────────────────
  nav_panel(
    title = "Asymmetric + \u03bb",
    layout_sidebar(
      sidebar = sidebar(
        width = 270, open = TRUE,
        h6("Model Parameters", class = "fw-bold mt-1", style = "color:#CC79A7"),
        sliderInput("asym_alpha_in",
                    HTML("Ingroup projection (\u03b1<sub>in</sub>)"),
                    min = 0.1, max = 10, value = 4.6, step = 0.1),
        sliderInput("asym_alpha_out",
                    HTML("Outgroup repulsion (\u03b1<sub>out</sub>)"),
                    min = 0.1, max = 10, value = 3.2, step = 0.1),
        sliderInput("asym_lambda",
                    HTML("Generalization sensitivity (\u03bb)"),
                    min = 0.1, max = 5, value = 3.5, step = 0.1),
        sliderInput("asym_bias",
                    HTML("Ingroup bias (\u03b3)"),
                    min = 0.01, max = 0.99, value = 0.35, step = 0.01),
        sliderInput("asym_w", "Lapse rate (w)",
                    min = 0, max = 0.99, value = 0.58, step = 0.01),
        hr(),
        actionButton("asym_reset",    "Reset to Study 1 estimates",
                     class = "btn-sm btn-outline-primary w-100 mb-2"),
        hr(),
        h6("Training Trait Structure", class = "fw-bold", style = "color:#CC79A7"),
        sliderInput("asym_r",
                    HTML("Self-rating coherence (<em>r</em>)"),
                    min = 0, max = 0.95, value = 0.70, step = 0.05),
        actionButton("asym_resample", "Resample training traits",
                     class = "btn-sm btn-outline-secondary w-100 mt-1")
      ),
      layout_columns(
        col_widths = c(5, 7),
        card(
          card_header("Training Traits"),
          plotOutput("asym_scatter", height = "390px"),
          card_footer(HTML(
            "<small>Each dot is a training trait. Dots are colored by self-rating.
             The dashed line shows the overall trend. Use the
             <em>Self-rating coherence (r)</em> slider to control how structured
             the self-concept is; click <em>Resample</em> to draw another set
             at the same coherence level.</small>"
          ))
        ),
        card(
          full_screen = TRUE,
          card_header("P(Ingroup Classification) \u2014 Asymmetric + \u03bb"),
          plotOutput("asym_plot", height = "390px"),
          card_footer(HTML(
            "<small>
             <b>X-axis:</b> semantic position of the test trait relative to the
             ingroup prototype (0 = maximally distant, 1 = maximally close),
             reflecting the self-evaluation-weighted similarity mean.
             <b>Dark line:</b> combined model.
             <b>Blue dashed:</b> ingroup-love only (\u03b1<sub>out</sub>\u00a0\u2192\u00a00).
             <b>Pink dashed:</b> outgroup-hate only (\u03b1<sub>in</sub>\u00a0\u2192\u00a00).
             Comparing the three lines reveals whether classification is driven
             by attraction to the ingroup or repulsion from the outgroup.
             </small>"
          ))
        )
      )
    )
  ),

  # ── Tab 3: About ─────────────────────────────────────────────────────────
  nav_panel(
    title = "About",
    card(
      card_header("About this App"),
      card_body(
        h4("Self-Anchoring as Similarity-Based Generalization"),
        p(HTML(
          "This interactive supplement accompanies <b>Elder et al. (2026)</b>,
           <em>TBD</em>."
        )),
        p("The app illustrates two cognitive architectures for self-anchoring:"),
        tags$ul(
          tags$li(
            strong("Symmetric + \u03bb:"),
            " A single projection rate (\u03b1) converts self-evaluations into
              group beliefs. Shepard\u2019s generalization gradient (\u03bb) controls
              how steeply inferences decay with semantic distance from the test trait."
          ),
          tags$li(
            strong("Asymmetric + \u03bb:"),
            " Separate rates for ingroup assimilation (\u03b1\u1d35\u2099) and
              outgroup repulsion (\u03b1\u2092\u1d64\u209c) test whether group identity
              is driven by \u2018ingroup love\u2019 or \u2018outgroup hate\u2019."
          )
        ),
        hr(),
        h5("Training Traits"),
        p(HTML(
          "Each tab samples 90 training traits at random — matching the number
           rated by participants in the actual studies. Each trait has a
           <em>semantic position</em> (0\u20131; higher = more similar to the
           ingroup prototype) and a <em>self-rating</em> (1 = not at all like me;
           7 = very much like me). The <em>Self-rating coherence (r)</em> slider
           controls how strongly self-ratings correlate with semantic position.
           At r\u00a0=\u00a00, ratings are random and the ingroup classification curve
           is flat; at r\u00a0=\u00a00.95, the self-concept is highly structured and
           the effect is maximal. The default (r\u00a0=\u00a00.70) reflects the
           empirical pattern in which people rate themselves more highly on traits
           conceptually close to traits they already see as self-descriptive.
           Adjusting r shows how self-concept coherence modulates the ingroup
           classification gradient; clicking <em>Resample</em> draws another
           random configuration at the same level."
        )),
        hr(),
        h5("Parameter Guide"),
        tags$dl(
          tags$dt(HTML("\u03b1 / \u03b1<sub>in</sub>, \u03b1<sub>out</sub>
                        \u2014 Projection rate(s)")),
          tags$dd("How strongly self-evaluations are converted into group beliefs.
                   Higher \u03b1 produces sharper, more extreme predictions."),
          tags$dt(HTML("\u03bb \u2014 Generalization sensitivity")),
          tags$dd("Steepness of the Shepard-style generalization gradient.
                   High \u03bb = narrow projection (only the most similar traits matter);
                   low \u03bb = broad projection (even distant traits contribute)."),
          tags$dt(HTML("\u03b3 \u2014 Ingroup bias")),
          tags$dd("Baseline probability of ingroup classification independent of
                   self-knowledge. Shifts the entire curve up or down."),
          tags$dt("w \u2014 Lapse rate"),
          tags$dd("Probability of a random response (P = 0.50). Compresses all
                   predictions toward chance, modeling inattention or indecision.")
        ),
        hr(),
        p(em(
          "Default parameter values are Study 1 group-level posterior medians
           from the hierarchical Bayesian model."
        ))
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Traits resample when r slider changes OR resample button is clicked.
  # Referencing both inputs creates the dependency automatically.
  sym_traits <- reactive({
    input$sym_resample                          # re-run on button click
    sample_traits(r_target = input$sym_r)       # re-run when r changes
  })
  asym_traits <- reactive({
    input$asym_resample
    sample_traits(r_target = input$asym_r)
  })

  # Reset parameters (not traits)
  observeEvent(input$sym_reset, {
    updateSliderInput(session, "sym_alpha",  value = 2.8)
    updateSliderInput(session, "sym_lambda", value = 3.6)
    updateSliderInput(session, "sym_bias",   value = 0.36)
    updateSliderInput(session, "sym_w",      value = 0.56)
  })
  observeEvent(input$asym_reset, {
    updateSliderInput(session, "asym_alpha_in",  value = 4.6)
    updateSliderInput(session, "asym_alpha_out", value = 3.2)
    updateSliderInput(session, "asym_lambda",    value = 3.5)
    updateSliderInput(session, "asym_bias",      value = 0.35)
    updateSliderInput(session, "asym_w",         value = 0.58)
  })

  # ── Scatter plots ─────────────────────────────────────────────────────────

  output$sym_scatter <- renderPlot({
    make_scatter(sym_traits(), accent = COL_SYM)
  }, res = 110)

  output$asym_scatter <- renderPlot({
    make_scatter(asym_traits(), accent = COL_ASYM_OUT)
  }, res = 110)

  # ── Symmetric P curve ─────────────────────────────────────────────────────

  output$sym_plot <- renderPlot({
    traits <- sym_traits()
    P <- sapply(S_seq, compute_sym,
                traits  = traits,
                alpha   = input$sym_alpha,
                lambda  = input$sym_lambda,
                bias    = input$sym_bias,
                w       = input$sym_w)

    ggplot(data.frame(x = S_seq, P = P), aes(x, P)) +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 colour = "grey55", linewidth = 0.45) +
      geom_ribbon(aes(ymin = 0.5, ymax = P), fill = COL_SYM, alpha = 0.12) +
      geom_line(colour = COL_SYM, linewidth = 1.5) +
      scale_y_continuous(labels = percent_format(1),
                         limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      labs(
        x = "Test Trait Semantic Position\n(0 = Distant from Prototype,  1 = Close)",
        y = "P(Ingroup Classification)",
        title = sprintf(
          "\u03b1 = %.1f  |  \u03bb = %.1f  |  \u03b3 = %.2f  |  w = %.2f",
          input$sym_alpha, input$sym_lambda, input$sym_bias, input$sym_w
        )
      ) +
      theme_app()
  }, res = 110)

  # ── Asymmetric P curve ────────────────────────────────────────────────────

  output$asym_plot <- renderPlot({
    traits <- asym_traits()

    run <- function(comp)
      sapply(S_seq, compute_asym,
             traits    = traits,
             alpha_in  = input$asym_alpha_in,
             alpha_out = input$asym_alpha_out,
             lambda    = input$asym_lambda,
             bias      = input$asym_bias,
             w         = input$asym_w,
             component = comp)

    lab_combined <- "Combined"
    lab_in       <- "Ingroup love only"
    lab_out      <- "Outgroup hate only"

    df <- rbind(
      data.frame(x = S_seq, P = run("combined"),      Line = lab_combined),
      data.frame(x = S_seq, P = run("ingroup_only"),  Line = lab_in),
      data.frame(x = S_seq, P = run("outgroup_only"), Line = lab_out)
    )
    df$Line <- factor(df$Line, levels = c(lab_combined, lab_in, lab_out))

    ggplot(df, aes(x, P, colour = Line, linewidth = Line, linetype = Line)) +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 colour = "grey55", linewidth = 0.45) +
      geom_line() +
      scale_colour_manual(values = c(
        "Combined"           = "grey20",
        "Ingroup love only"  = COL_ASYM_IN,
        "Outgroup hate only" = COL_ASYM_OUT
      )) +
      scale_linewidth_manual(values = c(
        "Combined"           = 1.5,
        "Ingroup love only"  = 0.85,
        "Outgroup hate only" = 0.85
      )) +
      scale_linetype_manual(values = c(
        "Combined"           = "solid",
        "Ingroup love only"  = "dashed",
        "Outgroup hate only" = "dashed"
      )) +
      scale_y_continuous(labels = percent_format(1),
                         limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      labs(
        x = "Test Trait Semantic Position\n(0 = Distant from Prototype,  1 = Close)",
        y = "P(Ingroup Classification)",
        title = sprintf(
          "\u03b1\u1d35\u2099 = %.1f  |  \u03b1\u2092\u1d64\u209c = %.1f  |  \u03bb = %.1f  |  \u03b3 = %.2f  |  w = %.2f",
          input$asym_alpha_in, input$asym_alpha_out,
          input$asym_lambda, input$asym_bias, input$asym_w
        )
      ) +
      guides(linewidth = "none") +
      theme_app() +
      theme(legend.position = "bottom")
  }, res = 110)
}

shinyApp(ui, server)

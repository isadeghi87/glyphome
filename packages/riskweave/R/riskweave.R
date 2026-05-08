## riskweave: Multi-Objective Treatment Decision Frontier
##
## Visualizes treatment selection as a Pareto problem with three axes:
##   - efficacy   (predicted response probability)
##   - toxicity   (predicted adverse event severity)
##   - confidence (model + evidence integration)
##
## Each candidate treatment is plotted on the efficacy-toxicity plane with
## confidence encoded as a ring radius. Pareto-dominant options are highlighted.
## A second panel shows modality-contribution rays per option, and a third
## panel shows decision sensitivity to single-modality perturbation.
##
## Figures generated:
##   A. Decision frontier        (efficacy x toxicity, confidence ring, Pareto)
##   B. Modality-contribution    (rays decomposing each option's recommendation)
##   C. Sensitivity drift        (how the recommended option shifts under perturbation)

# ── Palettes ─────────────────────────────────────────────────────────────────
riskweave_modality_pal <- function() {
  c(Genomic   = "#264653",
    Transcript= "#2a9d8f",
    Pathology = "#bc6c25",
    Imaging   = "#6a4c93",
    Clinical  = "#e76f51",
    Liquid    = "#457b9d")
}

riskweave_status_pal <- function() {
  c(Pareto = "#1d3557", Dominated = "#bdbdbd",
    Recommended = "#7a1f1f", Alternative = "#f4a261")
}

# ── Data generators ──────────────────────────────────────────────────────────
riskweave_options_data <- function(seed = 53) {
  set.seed(seed)
  treatments <- c("Chemo doublet","Chemo + IO","IO mono","TKI A","TKI B",
                  "PARP-i","Combo TKI+IO","Bispecific","Watchful wait")
  n <- length(treatments)
  efficacy <- c(0.55, 0.72, 0.40, 0.65, 0.50, 0.30, 0.80, 0.62, 0.18)
  toxicity <- c(0.65, 0.55, 0.25, 0.40, 0.45, 0.30, 0.78, 0.50, 0.05)
  confidence <- c(0.80, 0.75, 0.55, 0.85, 0.60, 0.50, 0.45, 0.40, 0.95)
  data.frame(
    treatment = treatments,
    efficacy = efficacy + stats::rnorm(n, 0, 0.02),
    toxicity = toxicity + stats::rnorm(n, 0, 0.02),
    confidence = confidence + stats::rnorm(n, 0, 0.02),
    stringsAsFactors = FALSE
  )
}

riskweave_modality_contributions <- function(options_df, seed = 67) {
  set.seed(seed)
  mods <- names(riskweave_modality_pal())
  do.call(rbind, lapply(seq_len(nrow(options_df)), function(i) {
    raw <- abs(stats::rnorm(length(mods), 0, 1))
    raw <- raw / sum(raw)
    data.frame(
      treatment = options_df$treatment[i],
      modality = mods,
      contribution = raw,
      direction = sample(c(1, -1), length(mods), replace = TRUE,
                         prob = c(0.7, 0.3)),
      stringsAsFactors = FALSE
    )
  }))
}

riskweave_sensitivity_data <- function(options_df, seed = 71) {
  set.seed(seed)
  mods <- names(riskweave_modality_pal())
  scenarios <- c("Baseline",
                 paste0("Drop ", mods),
                 "Add liquid biopsy")
  do.call(rbind, lapply(scenarios, function(sc) {
    eff <- options_df$efficacy
    tox <- options_df$toxicity
    if (grepl("^Drop ", sc)) {
      drop_mod <- sub("^Drop ", "", sc)
      # Drop a modality => reduce confidence and shift efficacy modestly
      bias <- stats::rnorm(nrow(options_df), 0, 0.04)
      eff <- pmin(1, pmax(0, eff + bias))
      tox <- pmin(1, pmax(0, tox + stats::rnorm(nrow(options_df), 0, 0.03)))
    } else if (sc == "Add liquid biopsy") {
      eff <- pmin(1, pmax(0, eff + stats::rnorm(nrow(options_df), 0.03, 0.02)))
      tox <- pmin(1, pmax(0, tox + stats::rnorm(nrow(options_df), -0.02, 0.02)))
    }
    score <- eff - 0.6 * tox
    rec <- options_df$treatment[which.max(score)]
    data.frame(scenario = sc,
               treatment = options_df$treatment,
               efficacy = eff, toxicity = tox,
               score = score,
               recommended = options_df$treatment == rec,
               stringsAsFactors = FALSE)
  }))
}

# Pareto front in (efficacy higher = better, toxicity lower = better)
riskweave_pareto <- function(df) {
  is_pareto <- logical(nrow(df))
  for (i in seq_len(nrow(df))) {
    dominated <- any(
      df$efficacy >= df$efficacy[i] &
      df$toxicity <= df$toxicity[i] &
      (df$efficacy > df$efficacy[i] | df$toxicity < df$toxicity[i])
    )
    is_pareto[i] <- !dominated
  }
  is_pareto
}

# ── Plot A: Decision frontier ────────────────────────────────────────────────
riskweave_plot_frontier <- function() {
  opts <- riskweave_options_data()
  opts$pareto <- riskweave_pareto(opts)
  opts$score <- opts$efficacy - 0.6 * opts$toxicity
  opts$recommended <- opts$treatment == opts$treatment[which.max(opts$score)]
  status_pal <- riskweave_status_pal()

  # Build the Pareto front line (sorted)
  pf <- opts[opts$pareto, ]
  pf <- pf[order(-pf$efficacy), ]

  ggplot(opts, aes(x = toxicity, y = efficacy)) +
    # iso-utility contours (efficacy - 0.6*toxicity = const)
    geom_abline(slope = 1 / 0.6,
                intercept = seq(-0.2, 1.2, by = 0.1) / 0.6 * 0.6,
                colour = "grey92", linewidth = 0.2) +
    # Pareto front
    geom_step(data = pf, direction = "vh",
              colour = status_pal["Pareto"], linewidth = 1.0,
              linetype = "dashed") +
    # confidence ring
    geom_point(aes(size = confidence * 18),
               shape = 21, fill = NA, colour = "grey60", stroke = 0.6) +
    # core point
    geom_point(aes(fill = ifelse(recommended, "Recommended",
                                 ifelse(pareto, "Pareto", "Dominated"))),
               shape = 21, colour = "grey15", size = 4.2, stroke = 0.4) +
    geom_text(aes(label = treatment),
              hjust = -0.15, vjust = -0.4, size = 3.0, colour = "grey20") +
    scale_fill_manual(values = status_pal, name = "Status") +
    scale_size_identity(guide = "none") +
    scale_x_continuous(limits = c(-0.05, 1.05)) +
    scale_y_continuous(limits = c(-0.05, 1.15)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10)) +
    labs(title = "A. Treatment decision frontier",
         subtitle = "Ring radius = model confidence; dashed step = Pareto frontier; recommended = highest efficacy − 0.6·toxicity",
         x = "Predicted toxicity", y = "Predicted efficacy")
}

# ── Plot B: Modality-contribution rays ───────────────────────────────────────
riskweave_plot_contributions <- function() {
  opts <- riskweave_options_data()
  cont <- riskweave_modality_contributions(opts)
  pal <- riskweave_modality_pal()

  # Order treatments by composite score (descending)
  opts$score <- opts$efficacy - 0.6 * opts$toxicity
  opts <- opts[order(-opts$score), ]
  ord <- opts$treatment
  cont$treatment <- factor(cont$treatment, levels = rev(ord))
  cont$signed <- cont$contribution * cont$direction

  ggplot(cont, aes(y = treatment, x = signed, fill = modality)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.35) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal, name = "Modality") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(-1, 1, by = 0.25)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10)) +
    labs(title = "B. Modality contributions to each option",
         subtitle = "Right = supports recommendation; left = argues against. Bar length = relative weight.",
         x = "Signed modality contribution",
         y = NULL)
}

# ── Plot C: Sensitivity drift ────────────────────────────────────────────────
riskweave_plot_sensitivity <- function() {
  opts <- riskweave_options_data()
  sens <- riskweave_sensitivity_data(opts)

  # Trace the *recommended* treatment per scenario as a path
  rec <- sens[sens$recommended, c("scenario","treatment","efficacy","toxicity")]

  baseline_rec <- rec$treatment[rec$scenario == "Baseline"]
  rec$is_baseline <- rec$scenario == "Baseline"
  rec$shifted <- rec$treatment != baseline_rec

  status_pal <- riskweave_status_pal()

  p_arrow <- ggplot(rec, aes(x = toxicity, y = efficacy)) +
    # all options as light backdrop
    geom_point(data = opts, aes(x = toxicity, y = efficacy),
               colour = "grey80", size = 2.5) +
    geom_text(data = opts, aes(label = treatment),
              colour = "grey70", size = 2.5, hjust = -0.1, vjust = -0.4) +
    # scenario points
    geom_point(aes(colour = shifted), size = 4.5) +
    geom_text(aes(label = scenario, colour = shifted),
              hjust = -0.15, vjust = 1.4, size = 3.0, fontface = "bold") +
    scale_colour_manual(values = c(`FALSE` = status_pal[["Pareto"]],
                                   `TRUE` = status_pal[["Recommended"]]),
                       name = "Recommendation shifted?") +
    scale_x_continuous(limits = c(-0.05, 1.05)) +
    scale_y_continuous(limits = c(-0.05, 1.15)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "C1. Sensitivity drift of the recommendation",
         x = "Predicted toxicity", y = "Predicted efficacy")

  # Score lollipop per scenario, only for the recommended option
  rec$scenario <- factor(rec$scenario, levels = rev(unique(rec$scenario)))
  p_lolli <- ggplot(rec, aes(y = scenario, x = efficacy - 0.6 * toxicity)) +
    geom_segment(aes(yend = scenario, x = 0,
                     xend = efficacy - 0.6 * toxicity,
                     colour = shifted), linewidth = 0.8) +
    geom_point(aes(colour = shifted), size = 3.5) +
    geom_text(aes(label = treatment), hjust = -0.1, size = 3.0,
              colour = "grey20") +
    scale_colour_manual(values = c(`FALSE` = status_pal[["Pareto"]],
                                   `TRUE` = status_pal[["Recommended"]]),
                       guide = "none") +
    scale_x_continuous(limits = c(0, NA), expand = c(0, 0.02)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "C2. Composite score under each perturbation",
         x = "Score = efficacy − 0.6·toxicity",
         y = NULL)

  patchwork::wrap_plots(p_arrow, p_lolli, nrow = 1, widths = c(1.3, 1)) +
    patchwork::plot_annotation(
      title = "C. Decision sensitivity to evidence perturbation",
      subtitle = "Each scenario drops one modality or adds liquid-biopsy evidence; red = the recommendation flips.",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13),
                             plot.subtitle = ggplot2::element_text(colour = "grey35", size = 10))
    )
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
riskweave_plot_demo <- function() riskweave_plot_frontier()

riskweave_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, "riskweave_01_frontier.png"),
                  plot = riskweave_plot_frontier(),
                  width = 10, height = 8, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "riskweave_02_contributions.png"),
                  plot = riskweave_plot_contributions(),
                  width = 11, height = 7, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "riskweave_03_sensitivity.png"),
                  plot = riskweave_plot_sensitivity(),
                  width = 14, height = 7, dpi = 300, bg = "white")
}

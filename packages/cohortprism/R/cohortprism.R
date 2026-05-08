## cohortprism: Spectral Patient Prisms
##
## Each patient is rendered as a vertical multi-band prism stacking:
##   - continuous molecular signature scores (color-banded)
##   - discrete mutation tick marks (overlaid pins)
##   - temporal clinical events (lane below the prism: relapse, switch, etc.)
##
## Replaces the OncoPrint discrete grid with continuous spectra plus a
## temporal event lane in one compact glyph. Designed for tumor board /
## molecular review where molecular state and clinical history must be
## inspected together.
##
## Figures generated:
##   A. Cohort prism wall   (patients ordered by similarity, response ribbon)
##   B. Single-patient drill (annotated prism with full event timeline)
##   C. Cohort projection   (case-vs-cohort multimodal signature scatter)

# ── Palettes ─────────────────────────────────────────────────────────────────
cohortprism_signature_pal <- function() {
  list(
    Proliferation = c("#f1faee", "#7a1f1f"),
    IFNg          = c("#f1faee", "#1d3557"),
    Stromal       = c("#f1faee", "#bc6c25"),
    Hypoxia       = c("#f1faee", "#6a4c93"),
    EMT           = c("#f1faee", "#264653"),
    AntigenPres   = c("#f1faee", "#2a9d8f")
  )
}

cohortprism_response_pal <- function() {
  c(CR = "#1d3557", PR = "#457b9d", SD = "#e9c46a", PD = "#e76f51")
}

cohortprism_event_pal <- function() {
  c(Diagnosis = "#264653", Therapy_start = "#2a9d8f",
    Therapy_switch = "#f4a261", Relapse = "#e76f51", Death = "#7a1f1f")
}

cohortprism_mutation_pal <- function() {
  c(Missense = "#264653", Truncating = "#e76f51",
    Inframe = "#f4a261", Splice = "#6a4c93", Fusion = "#2a9d8f")
}

# ── Data generators ──────────────────────────────────────────────────────────
cohortprism_cohort_data <- function(seed = 23) {
  set.seed(seed)
  n_pat <- 32L
  sigs <- c("Proliferation","IFNg","Stromal","Hypoxia","EMT","AntigenPres")
  # Three latent subgroups drive the signature profile
  groups <- sample(c("Hot","Cold","Mixed"), n_pat, replace = TRUE,
                   prob = c(0.4, 0.35, 0.25))
  group_means <- list(
    Hot   = c(Proliferation=0.7, IFNg=0.85, Stromal=0.30, Hypoxia=0.30, EMT=0.30, AntigenPres=0.80),
    Cold  = c(Proliferation=0.55,IFNg=0.20, Stromal=0.65, Hypoxia=0.55, EMT=0.55, AntigenPres=0.25),
    Mixed = c(Proliferation=0.50,IFNg=0.50, Stromal=0.50, Hypoxia=0.50, EMT=0.50, AntigenPres=0.50)
  )
  M <- do.call(rbind, lapply(seq_len(n_pat), function(i) {
    mu <- group_means[[groups[i]]]
    pmin(1, pmax(0, mu + stats::rnorm(length(mu), 0, 0.10)))
  }))
  M <- as.data.frame(M); colnames(M) <- sigs
  M$patient <- sprintf("P%02d", seq_len(n_pat))
  M$group <- groups
  # Response biased by IFNg vs Stromal axis
  resp_score <- M$IFNg - M$Stromal + stats::rnorm(n_pat, 0, 0.15)
  M$response <- cut(resp_score,
                    breaks = c(-Inf, -0.3, 0.0, 0.3, Inf),
                    labels = c("PD","SD","PR","CR"))
  M$response <- factor(as.character(M$response), levels = c("CR","PR","SD","PD"))
  M
}

cohortprism_mutations_data <- function(cohort, seed = 31) {
  set.seed(seed)
  genes <- c("TP53","KRAS","PIK3CA","EGFR","BRAF","ARID1A","STK11","PTEN")
  do.call(rbind, lapply(cohort$patient, function(p) {
    nm <- sample(0:4, 1, prob = c(0.10, 0.30, 0.30, 0.20, 0.10))
    if (nm == 0) return(data.frame(patient=character(0), gene=character(0), type=character(0)))
    g <- sample(genes, nm)
    t <- sample(c("Missense","Truncating","Inframe","Splice","Fusion"), nm,
                replace = TRUE, prob = c(0.55, 0.20, 0.10, 0.10, 0.05))
    data.frame(patient = p, gene = g, type = t, stringsAsFactors = FALSE)
  }))
}

cohortprism_events_data <- function(cohort, seed = 37) {
  set.seed(seed)
  do.call(rbind, lapply(cohort$patient, function(p) {
    rsp <- as.character(cohort$response[cohort$patient == p])
    grp <- cohort$group[cohort$patient == p]
    ev <- data.frame(patient = p, event = "Diagnosis", time = 0)
    ev <- rbind(ev, data.frame(patient = p, event = "Therapy_start",
                                time = stats::runif(1, 0.5, 2)))
    if (rsp %in% c("PD","SD")) {
      ev <- rbind(ev, data.frame(patient = p, event = "Therapy_switch",
                                  time = stats::runif(1, 4, 8)))
    }
    if (rsp == "PD" && stats::runif(1) < 0.6) {
      ev <- rbind(ev, data.frame(patient = p, event = "Relapse",
                                  time = stats::runif(1, 8, 14)))
    }
    if (rsp == "PD" && stats::runif(1) < 0.3) {
      ev <- rbind(ev, data.frame(patient = p, event = "Death",
                                  time = stats::runif(1, 14, 24)))
    }
    ev$group <- grp
    ev
  }))
}

# Order patients by response then molecular similarity within group
cohortprism_order <- function(cohort) {
  sigs <- c("Proliferation","IFNg","Stromal","Hypoxia","EMT","AntigenPres")
  M <- as.matrix(cohort[, sigs])
  rownames(M) <- cohort$patient
  ord_vec <- character(0)
  for (rsp in c("CR","PR","SD","PD")) {
    sub <- M[as.character(cohort$response) == rsp, , drop = FALSE]
    if (nrow(sub) <= 1) { ord_vec <- c(ord_vec, rownames(sub)); next }
    d <- stats::dist(sub)
    h <- stats::hclust(d, method = "average")
    ord_vec <- c(ord_vec, rownames(sub)[h$order])
  }
  ord_vec
}

# ── Plot A: Cohort prism wall ────────────────────────────────────────────────
cohortprism_plot_wall <- function() {
  cohort <- cohortprism_cohort_data()
  muts   <- cohortprism_mutations_data(cohort)
  ord    <- cohortprism_order(cohort)
  cohort$patient <- factor(cohort$patient, levels = ord)
  muts$patient   <- factor(muts$patient, levels = ord)

  sigs <- c("Proliferation","IFNg","Stromal","Hypoxia","EMT","AntigenPres")
  pal_list <- cohortprism_signature_pal()

  # Build per-band tile data
  band_h <- 1
  bands <- do.call(rbind, lapply(seq_along(sigs), function(i) {
    sg <- sigs[i]
    data.frame(patient = cohort$patient,
               band = sg,
               band_y = (length(sigs) - i + 1) * band_h,  # top-most = first
               value = cohort[[sg]],
               stringsAsFactors = FALSE)
  }))
  bands$band <- factor(bands$band, levels = sigs)

  # Build a single multi-segment "spectrum" using one fill per band by mapping
  # the value to its band-specific colour ramp at draw time.
  bands$fill <- vapply(seq_len(nrow(bands)), function(k) {
    rmp <- pal_list[[as.character(bands$band[k])]]
    grDevices::colorRampPalette(rmp)(101)[round(bands$value[k] * 100) + 1]
  }, character(1))

  # Mutation pins overlaid on the prism — placed above the prism in a small lane
  pin_lane_y <- length(sigs) * band_h + 0.6
  mut_pal <- cohortprism_mutation_pal()
  pins <- muts
  pins$y <- pin_lane_y + as.integer(factor(pins$gene)) * 0.18

  # Response ribbon below
  resp_y <- -0.8
  resp_pal <- cohortprism_response_pal()
  resp_df <- data.frame(patient = cohort$patient,
                        y = resp_y,
                        response = cohort$response)

  p <- ggplot() +
    geom_tile(data = bands, aes(x = patient, y = band_y, fill = fill),
              colour = "white", linewidth = 0.15, height = band_h * 0.95) +
    scale_fill_identity(guide = "none") +
    # band labels on the left
    annotate("text", x = 0.3, y = (length(sigs):1) * band_h,
             label = sigs, hjust = 1, size = 3.2, fontface = "bold",
             colour = "grey25") +
    # mutation pins
    geom_point(data = pins, aes(x = patient, y = y, colour = type),
               shape = 18, size = 2.3) +
    scale_colour_manual(values = mut_pal, name = "Mutation") +
    # gene labels for pins
    annotate("text", x = 0.3, y = pin_lane_y + (1:length(unique(muts$gene))) * 0.18,
             label = sort(unique(muts$gene)), hjust = 1, size = 2.5, colour = "grey30") +
    # response ribbon
    ggnewscale::new_scale_fill() +
    geom_tile(data = resp_df, aes(x = patient, y = y, fill = response),
              height = 0.5, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = resp_pal, name = "Response", drop = FALSE) +
    annotate("text", x = 0.3, y = resp_y, label = "Response",
             hjust = 1, size = 3.0, fontface = "bold", colour = "grey25") +
    coord_cartesian(xlim = c(-3.5, length(levels(cohort$patient)) + 0.5),
                    clip = "off") +
    theme_void(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          legend.position = "right",
          plot.margin = ggplot2::margin(10, 12, 10, 70),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "grey35", size = 10)) +
    labs(title = "A. Cohort prism wall",
         subtitle = "Each column = one patient. Vertical bands = continuous signatures; pins = mutations; bottom ribbon = response")
  p
}

# ── Plot B: Single-patient drill ─────────────────────────────────────────────
cohortprism_plot_single <- function(patient_id = "P05") {
  cohort <- cohortprism_cohort_data()
  muts   <- cohortprism_mutations_data(cohort)
  events <- cohortprism_events_data(cohort)
  if (!patient_id %in% cohort$patient) patient_id <- cohort$patient[1]
  pcoh <- cohort[cohort$patient == patient_id, , drop = FALSE]
  pmut <- muts[muts$patient == patient_id, , drop = FALSE]
  pev  <- events[events$patient == patient_id, , drop = FALSE]

  sigs <- c("Proliferation","IFNg","Stromal","Hypoxia","EMT","AntigenPres")
  pal_list <- cohortprism_signature_pal()

  band_df <- data.frame(
    band = factor(sigs, levels = sigs),
    band_y = (length(sigs):1),
    value = as.numeric(pcoh[, sigs]),
    stringsAsFactors = FALSE
  )
  band_df$fill <- vapply(seq_len(nrow(band_df)), function(k) {
    rmp <- pal_list[[as.character(band_df$band[k])]]
    grDevices::colorRampPalette(rmp)(101)[round(band_df$value[k] * 100) + 1]
  }, character(1))

  # Cohort distribution per band as background
  cohort_long <- do.call(rbind, lapply(sigs, function(sg) {
    data.frame(band = sg, value = cohort[[sg]],
               band_y = which(sigs == sg))
  }))
  cohort_long$band_y <- length(sigs) - cohort_long$band_y + 1

  p_prism <- ggplot() +
    # cohort jitter density behind
    geom_jitter(data = cohort_long,
                aes(x = value, y = band_y),
                width = 0, height = 0.10,
                colour = "grey75", size = 1.3, alpha = 0.6) +
    # the patient's value as a tall colored bar
    geom_segment(data = band_df,
                 aes(x = 0, xend = value, y = band_y, yend = band_y, colour = fill),
                 linewidth = 6) +
    geom_point(data = band_df, aes(x = value, y = band_y),
               shape = 21, fill = "white", colour = "grey15", size = 3.2, stroke = 0.6) +
    scale_colour_identity(guide = "none") +
    scale_y_continuous(breaks = length(sigs):1, labels = sigs) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 12)) +
    labs(title = sprintf("B1. Patient %s — signature prism", patient_id),
         x = "Signature score (vs cohort grey points)", y = NULL)

  # Mutation strip
  if (nrow(pmut) > 0) {
    pmut$y <- 1
    p_mut <- ggplot(pmut, aes(x = gene, y = y, colour = type)) +
      geom_point(shape = 18, size = 5) +
      geom_text(aes(label = gene), vjust = -1.5, size = 3.2, colour = "grey20") +
      scale_colour_manual(values = cohortprism_mutation_pal(), name = "Mutation") +
      scale_y_continuous(limits = c(0.5, 1.6)) +
      theme_void(base_size = 10) +
      theme(legend.position = "right",
            plot.title = element_text(face = "bold", size = 11)) +
      labs(title = "B2. Driver mutations")
  } else {
    p_mut <- ggplot() + theme_void() +
      annotate("text", x = 0, y = 0, label = "No driver mutations called",
               size = 3.5, colour = "grey40") +
      labs(title = "B2. Driver mutations")
  }

  # Timeline
  pev$event <- factor(pev$event,
                      levels = c("Diagnosis","Therapy_start","Therapy_switch",
                                 "Relapse","Death"))
  pev$y <- 1
  p_time <- ggplot(pev, aes(x = time, y = y, colour = event)) +
    geom_segment(aes(xend = time, yend = 0.6), linewidth = 0.5) +
    geom_point(size = 4) +
    geom_text(aes(label = event), angle = 35, hjust = 0, vjust = -0.4,
              size = 2.9, colour = "grey25", nudge_y = 0.05) +
    scale_colour_manual(values = cohortprism_event_pal(),
                        name = "Event", drop = FALSE) +
    scale_y_continuous(limits = c(0.4, 1.6)) +
    scale_x_continuous(name = "Months from diagnosis", limits = c(0, NA)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(title = "B3. Clinical event timeline")

  patchwork::wrap_plots(p_prism, p_mut, p_time, ncol = 1,
                        heights = c(2, 0.8, 1.2)) +
    patchwork::plot_annotation(
      title = sprintf("Single-patient drill-down: %s (response = %s)",
                      patient_id, as.character(pcoh$response)),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13))
    )
}

# ── Plot C: Cohort projection (case-vs-cohort) ──────────────────────────────
cohortprism_plot_projection <- function(highlight = "P05") {
  cohort <- cohortprism_cohort_data()
  sigs <- c("Proliferation","IFNg","Stromal","Hypoxia","EMT","AntigenPres")
  M <- as.matrix(cohort[, sigs])
  pc <- stats::prcomp(M, scale. = TRUE)
  proj <- as.data.frame(pc$x[, 1:2])
  colnames(proj) <- c("PC1","PC2")
  proj$patient <- cohort$patient
  proj$response <- cohort$response
  proj$is_case <- cohort$patient == highlight

  loadings <- as.data.frame(pc$rotation[, 1:2] * 4)
  loadings$signature <- rownames(loadings)
  colnames(loadings)[1:2] <- c("PC1","PC2")

  resp_pal <- cohortprism_response_pal()

  ggplot(proj, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(aes(fill = response, size = is_case),
               shape = 21, colour = "grey20", stroke = 0.3, alpha = 0.85) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = grid::arrow(length = grid::unit(2, "mm")),
                 colour = "grey35", linewidth = 0.5) +
    geom_text(data = loadings, aes(x = PC1 * 1.12, y = PC2 * 1.12, label = signature),
              size = 3.0, colour = "grey25", fontface = "bold") +
    geom_text(data = proj[proj$is_case, ],
              aes(label = patient), nudge_y = 0.4,
              size = 4, fontface = "bold", colour = "#7a1f1f") +
    scale_fill_manual(values = resp_pal, name = "Response", drop = FALSE) +
    scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 7), guide = "none") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10)) +
    labs(title = "C. Case-vs-cohort projection",
         subtitle = sprintf("PCA over signatures; highlighted case = %s. Arrows = signature loadings.",
                            highlight),
         x = "PC1", y = "PC2")
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
cohortprism_plot_demo <- function() cohortprism_plot_wall()

cohortprism_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, "cohortprism_01_wall.png"),
                  plot = cohortprism_plot_wall(),
                  width = 14, height = 7, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "cohortprism_02_single.png"),
                  plot = cohortprism_plot_single("P05"),
                  width = 11, height = 9, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "cohortprism_03_projection.png"),
                  plot = cohortprism_plot_projection("P05"),
                  width = 10, height = 8, dpi = 300, bg = "white")
}

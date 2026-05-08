## uncertainr: Uncertainty-Native Bioinformatics Visualization
##
## Visual grammar: Every biological mark is a compound glyph ("certainty card")
## that fuses result + three uncertainty dimensions into one atomic visual object.
##
## Uncertainty types encoded:
##   1. Measurement noise     — sequencing depth, dropout rate
##   2. Model uncertainty     — posterior variance, bootstrap SD
##   3. Annotation ambiguity  — label entropy, multi-mapping score
##
## Figures generated:
##   A. Cell-type UMAP: certainty halos + decomposed panels + distribution
##   B. Branching trajectory: confidence tubes + decomposition along pseudotime
##   C. Pathway enrichment: lollipop + heatmap decomposition + fragility scatter

# ── Palettes ─────────────────────────────────────────────────────────────────
uncertainr_celltype_palette <- function() {
  c("CD8 T" = "#264653", "CD4 T" = "#2a9d8f", "Treg" = "#e9c46a",
    "B cell" = "#577590", "Plasma" = "#f4a261", "Monocyte" = "#e76f51",
    "cDC" = "#bc6c25", "pDC" = "#6a994e", "NK" = "#d62828",
    "Cycling" = "#7b2cbf", "Progenitor" = "#8ecae6", "Exhausted" = "#6d597a",
    "Memory" = "#1d3557", "Effector" = "#fb8500", "Naive" = "#90be6d")
}

uncertainr_unc_scales <- function() {
  list(
    measurement = scales::col_numeric(c("#f8f9fa", "#e63946"), domain = c(0, 1)),
    model       = scales::col_numeric(c("#f8f9fa", "#1d3557"), domain = c(0, 1)),
    annotation  = scales::col_numeric(c("#f8f9fa", "#7b2cbf"), domain = c(0, 1))
  )
}

# ── Data generators ──────────────────────────────────────────────────────────
uncertainr_celltype_data <- function(seed = 42) {
  set.seed(seed)
  cc <- data.frame(
    ct = c("CD8 T","CD4 T","Treg","B cell","Plasma","Monocyte","cDC","pDC","NK","Cycling"),
    cx = c(-4.2,-2.8,-1.5,3.2,4.8,1.0,2.5,3.8,-5.5,0),
    cy = c(2,3.5,4.2,-1.5,-2.8,-3.5,-1,-0.2,0.5,5.5),
    sp = c(.55,.65,.40,.70,.50,.80,.45,.35,.60,.50),
    n  = c(80L,90L,30L,70L,40L,85L,35L,20L,60L,40L), stringsAsFactors = FALSE)

  do.call(rbind, lapply(seq_len(nrow(cc)), function(i) {
    r <- cc[i,]; nn <- r$n
    dx <- stats::rnorm(nn, 0, r$sp); dy <- stats::rnorm(nn, 0, r$sp)
    d_from_center <- sqrt(dx^2 + dy^2) / r$sp
    data.frame(
      UMAP_1 = r$cx + dx, UMAP_2 = r$cy + dy, celltype = r$ct,
      measurement_noise  = pmin(1, pmax(0, 0.12 + 0.30 * stats::runif(nn) + 0.18 * d_from_center)),
      model_uncertainty  = pmin(1, pmax(0, stats::rbeta(nn, ifelse(nn < 40, 2, 5), ifelse(nn < 40, 5, 12)))),
      annotation_ambiguity = pmin(1, pmax(0, ifelse(r$ct %in% c("Treg","Cycling","cDC"),
                                                     stats::rbeta(nn, 3, 4), stats::rbeta(nn, 2, 10)))),
      stringsAsFactors = FALSE)
  }))
}

uncertainr_trajectory_data <- function(seed = 7) {
  set.seed(seed)
  n <- 120L; pt <- seq(0, 1, length.out = n)
  tx <- pt * 6 - 3; ty <- sin(pt * pi * 0.8) * 1.5; bp <- round(0.55 * n)
  make_seg <- function(x, y, t, branch, state) {
    nn <- length(x)
    data.frame(x = x + stats::rnorm(nn, 0, .08), y = y + stats::rnorm(nn, 0, .08),
               pseudotime = t, branch = branch, state = state,
               measurement_noise  = pmin(1, pmax(0, .10 + .20*t + stats::runif(nn, 0, .15))),
               model_uncertainty  = pmin(1, pmax(0, .08 + .60*(t-.4)^2 + stats::runif(nn, 0, .12))),
               annotation_ambiguity = pmin(1, pmax(0, ifelse(t > .45, .2 + .50*(t-.45), .05 + stats::runif(nn, 0, .10)))),
               stringsAsFactors = FALSE)
  }
  trunk_st <- ifelse(pt[1:bp] < .25, "Progenitor", ifelse(pt[1:bp] < .45, "Naive", "Memory"))
  bax <- tx[bp:n] + (pt[bp:n] - .55) * 2.5; bay <- ty[bp] + (pt[bp:n] - .55) * 3.2
  bbx <- tx[bp:n] + (pt[bp:n] - .55) * 1.8; bby <- ty[bp] - (pt[bp:n] - .55) * 2.8
  rbind(make_seg(tx[1:bp], ty[1:bp], pt[1:bp], "Trunk", trunk_st),
        make_seg(bax, bay, pt[bp:n], "Effector", rep("Effector", length(bax))),
        make_seg(bbx, bby, pt[bp:n], "Exhausted", rep("Exhausted", length(bbx))))
}

uncertainr_pathway_data <- function(seed = 99) {
  set.seed(seed)
  pw <- c("IFN-alpha","IFN-gamma","TNF/NF-kB","OxPhos","MYC V1","E2F targets",
          "G2M checkpoint","p53 pathway","Apoptosis","IL6-JAK-STAT3","Complement",
          "Coagulation","EMT","Angiogenesis","Hypoxia","DNA repair",
          "UPR","mTORC1","Glycolysis","FA metabolism")
  nes <- c(3.2,2.8,2.1,-1.9,1.7,1.5,1.3,-2.4,-1.1,2.5,.9,.4,-.8,1.0,-1.6,1.8,.6,1.4,-.7,-1.2)
  data.frame(pathway = factor(pw, levels = rev(pw)), NES = nes,
             measurement_noise = pmin(1,pmax(0, abs(nes)*.08+stats::runif(20,.05,.30))),
             model_uncertainty = pmin(1,pmax(0, 1-abs(nes)*.22+stats::runif(20,-.10,.20))),
             annotation_ambiguity = pmin(1,pmax(0, ifelse(abs(nes)<1.2, stats::runif(20,.35,.75), stats::runif(20,.02,.25)))),
             padj = 10^(-abs(nes)*1.5+stats::runif(20,-.5,.5)),
             direction = ifelse(nes>0,"Up","Down"), stringsAsFactors = FALSE)
}

# ── Plot A: Cell type UMAP ───────────────────────────────────────────────────
uncertainr_plot_celltype <- function() {
  d <- uncertainr_celltype_data(); pal <- uncertainr_celltype_palette()
  d$certainty <- 1 - (.35*d$measurement_noise + .40*d$model_uncertainty + .25*d$annotation_ambiguity)
  d$halo_sz <- 3 + d$annotation_ambiguity * 6
  d$celltype_f <- factor(d$celltype, levels = names(pal)[names(pal) %in% unique(d$celltype)])

  # -- Main UMAP with halos + cell-type legend --
  p_ov <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP_1, y = UMAP_2)) +
    ggplot2::geom_point(ggplot2::aes(size = halo_sz), color = "#7b2cbf", alpha = .12, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(fill = celltype_f), shape = 21, size = 2.2, stroke = .2, color = "#2b2d42") +
    ggplot2::scale_fill_manual(values = pal, name = "Cell type") +
    ggplot2::scale_size_identity() +
    ggplot2::coord_equal() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.title = ggplot2::element_text(face = "bold", size = 10),
                   legend.text = ggplot2::element_text(size = 9),
                   plot.title = ggplot2::element_text(face = "bold", size = 12),
                   plot.subtitle = ggplot2::element_text(size = 10, color = "#555555")) +
    ggplot2::labs(title = "A. Cell-type UMAP with certainty halos",
                  subtitle = "Core = assigned type; halo radius = annotation ambiguity")

  # -- Decomposed uncertainty mini-panels with proper gradient legends --
  mkp <- function(col, title, lo_col, hi_col) {
    ggplot2::ggplot(d, ggplot2::aes(x = UMAP_1, y = UMAP_2)) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[col]]), size = 1.4, alpha = .8) +
      ggplot2::scale_color_gradient(low = "#f0f0f0", high = hi_col,
                                    limits = c(0, 1), name = NULL,
                                    breaks = c(0, 0.5, 1), labels = c("Low", "", "High")) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 10, face = "bold", color = hi_col),
                     legend.position = "bottom",
                     legend.key.width = ggplot2::unit(0.8, "cm"),
                     legend.key.height = ggplot2::unit(0.25, "cm"),
                     legend.text = ggplot2::element_text(size = 7)) +
      ggplot2::labs(title = title)
  }
  p_m <- mkp("measurement_noise", "Measurement noise", "#f8f9fa", "#e63946")
  p_o <- mkp("model_uncertainty", "Model uncertainty", "#f8f9fa", "#1d3557")
  p_a <- mkp("annotation_ambiguity", "Annotation ambiguity", "#f8f9fa", "#7b2cbf")

  # -- Boxplot with matching colors + jittered points --
  d$ct_f <- factor(d$celltype, levels = names(sort(tapply(d$certainty, d$celltype, stats::median))))
  p_dist <- ggplot2::ggplot(d, ggplot2::aes(x = ct_f, y = certainty, fill = ct_f)) +
    ggplot2::geom_boxplot(alpha = .65, outlier.shape = NA, linewidth = .4, width = .7) +
    ggplot2::geom_jitter(size = .4, alpha = .25, width = .15, color = "#2b2d42") +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", size = 11),
                   axis.text.y = ggplot2::element_text(size = 9)) +
    ggplot2::labs(x = NULL, y = "Composite certainty (1 = fully certain)", title = "B. Certainty by cell type")

  patchwork::wrap_plots(
    p_ov,
    patchwork::wrap_plots(p_m, p_o, p_a, nrow = 1),
    p_dist,
    nrow = 3, heights = c(3, 1.4, 1.5)
  ) + patchwork::plot_annotation(
    title = "uncertainr: Uncertainty-native cell type assignment",
    subtitle = "Each cell carries a tri-dimensional uncertainty signature (measurement \u00b7 model \u00b7 annotation)",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#444444")
    )
  )
}

# ── Plot B: Trajectory ───────────────────────────────────────────────────────
uncertainr_plot_trajectory <- function() {
  d <- uncertainr_trajectory_data(); pal <- uncertainr_celltype_palette()
  d$total_unc <- .35*d$measurement_noise + .40*d$model_uncertainty + .25*d$annotation_ambiguity
  d$tube_w <- 1.5 + d$total_unc * 4
  d$state_f <- factor(d$state, levels = c("Progenitor","Naive","Memory","Effector","Exhausted"))
  d$branch_f <- factor(d$branch, levels = c("Trunk","Effector","Exhausted"))

  # -- Compute smoothed backbone per branch for path ribbon --
  branches <- split(d, d$branch)
  backbone <- do.call(rbind, lapply(branches, function(bd) {
    bd <- bd[order(bd$pseudotime), ]
    sm_x <- stats::loess(x ~ pseudotime, data = bd, span = .4)
    sm_y <- stats::loess(y ~ pseudotime, data = bd, span = .4)
    data.frame(x = stats::fitted(sm_x), y = stats::fitted(sm_y),
               branch = bd$branch[1], pseudotime = bd$pseudotime,
               stringsAsFactors = FALSE)
  }))

  p_main <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    # Backbone path
    ggplot2::geom_path(data = backbone, ggplot2::aes(group = branch),
                       color = "#adb5bd", linewidth = 1.2, alpha = .5) +
    # Uncertainty halos
    ggplot2::geom_point(ggplot2::aes(size = tube_w), color = "#9d4edd", alpha = .08, show.legend = FALSE) +
    # Core points colored by state
    ggplot2::geom_point(ggplot2::aes(fill = state_f, size = tube_w * 0.5),
                        shape = 21, stroke = .2, color = "#2b2d42", alpha = .75) +
    ggplot2::scale_fill_manual(values = pal, name = "State",
                               drop = FALSE) +
    ggplot2::scale_size_identity() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold", size = 12),
                   plot.subtitle = ggplot2::element_text(size = 10, color = "#555555")) +
    ggplot2::labs(x = "Trajectory dimension 1", y = "Trajectory dimension 2",
                  title = "A. Branching trajectory with confidence tubes",
                  subtitle = "Point size = total uncertainty; color = differentiation state; grey backbone = smoothed path")

  dec <- data.frame(
    pseudotime = rep(d$pseudotime, 3),
    branch = factor(rep(d$branch, 3), levels = c("Trunk","Effector","Exhausted")),
    utype = rep(c("Measurement","Model","Annotation"), each = nrow(d)),
    value = c(d$measurement_noise, d$model_uncertainty, d$annotation_ambiguity))
  p_dec <- ggplot2::ggplot(dec, ggplot2::aes(x = pseudotime, y = value, color = utype)) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, alpha = .15, linewidth = .9, span = .4) +
    ggplot2::scale_color_manual(values = c(Measurement = "#e63946", Model = "#1d3557", Annotation = "#7b2cbf")) +
    ggplot2::facet_wrap(~branch, nrow = 1) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "top",
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", size = 11),
                   plot.title = ggplot2::element_text(face = "bold", size = 11)) +
    ggplot2::labs(x = "Pseudotime", y = "Uncertainty", color = NULL,
                  title = "B. Uncertainty decomposition along pseudotime")

  patchwork::wrap_plots(p_main, p_dec, nrow = 2, heights = c(2, 1.2)) +
    patchwork::plot_annotation(
      title = "uncertainr: Trajectory inference with structured uncertainty",
      subtitle = "Confidence tubes swell where branches diverge; decomposition reveals dominant uncertainty type",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 11, color = "#444444")
      )
    )
}

# ── Plot C: Pathway enrichment ───────────────────────────────────────────────
uncertainr_plot_pathway <- function() {
  d <- uncertainr_pathway_data()
  d$sig_lbl <- ifelse(d$padj < .05, "\u2605", "")
  d$direction_f <- factor(d$direction, levels = c("Up","Down"))

  # -- Lollipop with proper segments + up/down legend --
  p_main <- ggplot2::ggplot(d, ggplot2::aes(x = NES, y = pathway)) +
    ggplot2::geom_vline(xintercept = 0, color = "#adb5bd", linewidth = .5) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, y = pathway, yend = pathway),
                          color = "#dee2e6", linewidth = .6) +
    ggplot2::geom_point(ggplot2::aes(fill = direction_f, size = -log10(padj)),
                        shape = 21, stroke = .3, color = "#2b2d42") +
    ggplot2::geom_text(ggplot2::aes(label = sig_lbl),
                       nudge_x = ifelse(d$NES > 0, .2, -.2), size = 3.5, color = "#d00000") +
    ggplot2::scale_fill_manual(values = c(Up = "#2a9d8f", Down = "#e76f51"), name = "Direction") +
    ggplot2::scale_size_continuous(range = c(2, 5.5), name = expression(-log[10]*"(padj)")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.box = "horizontal",
                   plot.title = ggplot2::element_text(face = "bold", size = 11)) +
    ggplot2::labs(x = "NES", y = NULL, title = "A. Pathway enrichment lollipop",
                  subtitle = "Point size = significance; \u2605 = padj < 0.05")

  # -- Uncertainty heatmap with better color range --
  ul <- data.frame(
    pathway = rep(d$pathway, 3),
    type = factor(rep(c("Measurement","Model","Annotation"), each = nrow(d)),
                  levels = c("Measurement","Model","Annotation")),
    value = c(d$measurement_noise, d$model_uncertainty, d$annotation_ambiguity))
  ul$pathway <- factor(ul$pathway, levels = levels(d$pathway))
  p_unc <- ggplot2::ggplot(ul, ggplot2::aes(x = type, y = pathway, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = .4) +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1,
                                  limits = c(0, 1), name = "Uncertainty",
                                  breaks = c(0, .5, 1)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 8),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(.6, "cm"),
                   plot.title = ggplot2::element_text(face = "bold", size = 11)) +
    ggplot2::labs(title = "B. Uncertainty decomposition")

  # -- Fragility scatter with non-overlapping labels --
  d$fragility <- d$model_uncertainty * (1 - abs(d$NES) / max(abs(d$NES)))
  d$abs_nes <- abs(d$NES)
  d$frag_lbl <- ifelse(d$fragility > .35, as.character(d$pathway), "")

  p_frag <- ggplot2::ggplot(d, ggplot2::aes(x = abs_nes, y = fragility)) +
    ggplot2::geom_hline(yintercept = .35, linetype = "dashed", color = "#d00000", linewidth = .5, alpha = .6) +
    ggplot2::annotate("text", x = max(d$abs_nes), y = .37, label = "Fragility threshold",
                      hjust = 1, size = 2.8, color = "#d00000", fontface = "italic") +
    ggplot2::geom_point(ggplot2::aes(fill = direction_f), shape = 21, size = 3, stroke = .2, color = "#343a40") +
    ggplot2::geom_text(ggplot2::aes(label = frag_lbl), size = 2.5, nudge_y = .025,
                       check_overlap = TRUE, color = "#333333") +
    ggplot2::scale_fill_manual(values = c(Up = "#2a9d8f", Down = "#e76f51"), guide = "none") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 11)) +
    ggplot2::labs(x = "|NES|", y = "Fragility index",
                  title = "C. Fragility: pathways at risk of flipping")

  patchwork::wrap_plots(
    patchwork::wrap_plots(p_main, p_unc, nrow = 1, widths = c(2.5, 1)),
    p_frag, nrow = 2, heights = c(2.5, 1.2)
  ) + patchwork::plot_annotation(
    title = "uncertainr: Pathway enrichment with tri-uncertainty decomposition",
    subtitle = "Heatmap shows where each result is vulnerable; scatter plots fragility index",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#444444")
    )
  )
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
uncertainr_plot_demo <- function() { uncertainr_plot_celltype() }

uncertainr_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir,"uncertainr_01_celltype_umap.png"),
                  plot=uncertainr_plot_celltype(), width=11,height=12,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"uncertainr_02_trajectory.png"),
                  plot=uncertainr_plot_trajectory(), width=13,height=9,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"uncertainr_03_pathway.png"),
                  plot=uncertainr_plot_pathway(), width=13,height=10,dpi=300,bg="white")
}

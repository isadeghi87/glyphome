## scalefusion: Multi-Resolution Bioinformatics Visualization
##
## Visual grammar: "Scale ladders" that show how a biological signal behaves when
## viewed at different resolutions—from individual cells → subclusters →
## cell types → compartments → tissue regions → patient-level summaries.
##
## Key idea: expose "semantic drift" — when aggregating reveals or hides biology.
## Each observation carries a glyph that encodes its identity at 3+ resolution levels.
##
## Figures generated:
##   A. Gene expression scale ladder: expression across 4 resolution levels with
##      semantic continuity ribbons
##   B. Spatial composition zoom: concentric view from spot → neighborhood →
##      region → section
##   C. Cross-scale disagreement heatmap: which genes change behavior at coarser
##      resolution

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Palettes ─────────────────────────────────────────────────────────────────
scalefusion_res_pal <- function() {
  c("Cell"="#264653","Subcluster"="#2a9d8f","CellType"="#e9c46a","Compartment"="#e76f51")
}
scalefusion_ct_pal <- function() {
  c("CD8 T"="#264653","CD4 T"="#2a9d8f","Treg"="#e9c46a","B cell"="#577590",
    "Monocyte"="#e76f51","NK"="#d62828","Cycling"="#7b2cbf","DC"="#bc6c25")
}
scalefusion_spatial_pal <- function() {
  c("Tumor core"="#e76f51","Tumor margin"="#f4a261","Stroma"="#2a9d8f",
    "Immune hub"="#264653","Necrotic"="#6c757d","Normal"="#8ecae6")
}

# ── Data: Gene expression across scales ──────────────────────────────────────
scalefusion_expr_data <- function(seed = 42) {
  set.seed(seed)
  celltypes <- c("CD8 T","CD4 T","Treg","B cell","Monocyte","NK","Cycling","DC")
  compartments <- c("Lymphoid","Lymphoid","Lymphoid","Lymphoid","Myeloid","Lymphoid","Cycling","Myeloid")
  genes <- c("CD3E","CD8A","FOXP3","MS4A1","CD14","NKG7","MKI67","HLA-DRA",
             "GZMA","PDCD1","IL2RA","CD79A","LYZ","KLRB1","TOP2A","CST3")
  ng <- length(genes); nct <- length(celltypes)

  # True expression means per gene per celltype
  expr_mat <- matrix(stats::rnorm(ng*nct, 1, .5), ng, nct,
                     dimnames=list(genes, celltypes))
  # markers get high expression in their cell type
  marker_map <- list("CD8 T"=c("CD3E","CD8A","GZMA","PDCD1"),
                     "CD4 T"=c("CD3E","IL2RA"),
                     "Treg"=c("FOXP3","IL2RA"),
                     "B cell"=c("MS4A1","CD79A"),
                     "Monocyte"=c("CD14","LYZ","CST3"),
                     "NK"=c("NKG7","KLRB1"),
                     "Cycling"=c("MKI67","TOP2A"),
                     "DC"=c("HLA-DRA","CST3"))
  for (ct in celltypes) for (g in marker_map[[ct]]) expr_mat[g,ct] <- expr_mat[g,ct]+3

  # generate cells
  cells_per_ct <- c(80,90,25,60,70,50,30,25)
  rows <- list()
  for (i in seq_along(celltypes)) {
    nc <- cells_per_ct[i]; ct <- celltypes[i]; comp <- compartments[i]
    scl <- sample(paste0(ct,"_sc",1:3), nc, replace=TRUE)
    for (g in genes) {
      mu <- expr_mat[g,ct]; expr_vals <- pmax(0, stats::rnorm(nc, mu, .6))
      rows[[length(rows)+1]] <- data.frame(
        cell_id=paste0(ct,"_",1:nc), gene=g, expr=expr_vals,
        celltype=ct, subcluster=scl, compartment=comp, stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

# ── Data: Spatial multi-scale ────────────────────────────────────────────────
scalefusion_spatial_data <- function(seed = 7) {
  set.seed(seed)
  regions <- c("Tumor core","Tumor margin","Stroma","Immune hub","Necrotic","Normal")
  nx <- 25; ny <- 25; n <- nx * ny
  gx <- rep(1:nx, ny); gy <- rep(1:ny, each=nx)
  # region assignment: spatial blocks
  region <- character(n)
  region[gx<=10 & gy<=10] <- "Tumor core"; region[gx<=14 & gx>10 & gy<=10] <- "Tumor margin"
  region[gx>14 & gy<=10] <- "Stroma"; region[gy>10 & gy<=18] <- "Stroma"
  region[gx<=8 & gy>18] <- "Immune hub"; region[gx>8 & gy>10 & gy<=18] <- "Tumor margin"
  region[gx>18 & gy>10] <- "Normal"; region[region==""] <- "Necrotic"
  # cell type composition per spot
  ct_names <- c("CD8 T","CD4 T","Monocyte","Tumor","Stroma","NK","DC","B cell")
  comp <- matrix(0, n, length(ct_names), dimnames=list(NULL,ct_names))
  base_comp <- list("Tumor core"=c(.05,.02,.05,.70,.10,.02,.01,.05),
                     "Tumor margin"=c(.15,.08,.12,.35,.15,.05,.03,.07),
                     "Stroma"=c(.05,.05,.08,.10,.55,.03,.04,.10),
                     "Immune hub"=c(.25,.15,.15,.02,.05,.15,.08,.15),
                     "Necrotic"=c(.02,.01,.02,.50,.30,.01,.01,.13),
                     "Normal"=c(.05,.05,.05,.05,.60,.05,.05,.10))
  for (i in 1:n) {
    bc <- base_comp[[region[i]]]; noise <- pmax(0, bc + stats::rnorm(length(bc), 0, .03))
    comp[i,] <- noise / sum(noise)
  }
  # neighborhood averaging (3x3)
  nb_comp <- comp
  for (i in 1:n) {
    xi <- gx[i]; yi <- gy[i]
    nbrs <- which(abs(gx-xi)<=1 & abs(gy-yi)<=1)
    nb_comp[i,] <- colMeans(comp[nbrs,,drop=FALSE])
  }
  list(spots=data.frame(x=gx,y=gy,region=region,stringsAsFactors=FALSE),
       comp=comp, nb_comp=nb_comp, ct_names=ct_names)
}

# ── Plot A: Gene expression scale ladder ─────────────────────────────────────
scalefusion_plot_expr <- function() {
  d <- scalefusion_expr_data()
  genes_show <- c("CD8A","FOXP3","CD14","NKG7","MKI67","HLA-DRA")
  ds <- d[d$gene %in% genes_show, ]
  res_pal <- scalefusion_res_pal()

  # Cell-level: violin per gene per celltype
  ds$ct_gene <- paste0(ds$celltype, "|", ds$gene)
  p_cell <- ggplot(ds, aes(x = gene, y = expr, fill = celltype)) +
    geom_violin(scale = "width", linewidth = .2, alpha = .55, position = position_dodge(width = .8)) +
    scale_fill_manual(values = scalefusion_ct_pal()) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = NULL, y = "Expression", fill = NULL,
         title = "A1. Cell-level: full distribution per cell type")

  # Subcluster level: mean + SE
  sc_agg <- ds %>% group_by(gene,celltype,subcluster) %>%
    summarise(mean_e=mean(expr),se_e=stats::sd(expr)/sqrt(n()),.groups="drop")
  p_sc <- ggplot(sc_agg, aes(x = gene, y = mean_e, color = celltype, group = subcluster)) +
    geom_point(position = position_jitterdodge(jitter.width = .1, dodge.width = .6), size = 1.8, alpha = .7) +
    geom_errorbar(aes(ymin = mean_e - se_e, ymax = mean_e + se_e),
                  position = position_jitterdodge(jitter.width = .1, dodge.width = .6), width = .1, linewidth = .3, alpha = .5) +
    scale_color_manual(values = scalefusion_ct_pal()) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = NULL, y = "Mean expr \u00b1 SE", title = "A2. Subcluster-level: means + SE")

  # CellType level: single dot per celltype per gene
  ct_agg <- ds %>% group_by(gene,celltype) %>% summarise(mean_e=mean(expr),.groups="drop")
  p_ct <- ggplot(ct_agg, aes(x = gene, y = mean_e, fill = celltype)) +
    geom_col(position = position_dodge(width = .7), width = .6, color = "white", linewidth = .2) +
    scale_fill_manual(values = scalefusion_ct_pal()) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = NULL, y = "Mean expression", title = "A3. Cell-type-level: aggregated bars")

  # Compartment level: lose individual cell types
  ds$compartment <- factor(ds$compartment)
  comp_agg <- ds %>% group_by(gene,compartment) %>% summarise(mean_e=mean(expr),.groups="drop")
  p_comp <- ggplot(comp_agg, aes(x = gene, y = mean_e, fill = compartment)) +
    geom_col(position = position_dodge(width = .6), width = .5, color = "white", linewidth = .2) +
    scale_fill_manual(values = c("Lymphoid" = "#264653", "Myeloid" = "#e76f51", "Cycling" = "#7b2cbf")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = NULL, y = "Mean expression", fill = NULL, title = "A4. Compartment-level: high-level summary")

  patchwork::wrap_plots(p_cell, p_sc, p_ct, p_comp, nrow = 4, heights = c(1.5, 1, 1, 1)) +
    patchwork::plot_annotation(
      title = "scalefusion: Gene expression across 4 resolution levels",
      subtitle = "Each row zooms out \u2014 detail lost vs. pattern gained visible at a glance",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "#444444")
      )
    )
}

# ── Plot B: Spatial multi-scale ──────────────────────────────────────────────
scalefusion_plot_spatial <- function() {
  sp <- scalefusion_spatial_data()
  spots <- sp$spots; comp <- sp$comp; nb <- sp$nb_comp; ct <- sp$ct_names
  pal <- scalefusion_spatial_pal()

  # spot-level composition: dominant cell type
  dom_ct <- ct[apply(comp, 1, which.max)]
  spots$dominant <- dom_ct
  spots$entropy <- -rowSums(comp * log2(comp + 1e-9))
  spots$CD8_frac <- comp[,"CD8 T"]

  p_region <- ggplot(spots, aes_string(x="x",y="y",fill="region")) +
    geom_tile(color="white",linewidth=.1) + scale_fill_manual(values=pal) +
    coord_equal() + theme_void(base_size=9) +
    theme(legend.position="right") + labs(title="B1. Region annotation",fill="Region")

  p_dom <- ggplot(spots, aes_string(x="x",y="y",fill="dominant")) +
    geom_tile(color="white",linewidth=.1) + scale_fill_manual(values=scales::hue_pal()(length(ct))) +
    coord_equal() + theme_void(base_size=9) +
    theme(legend.position="right") + labs(title="B2. Spot-level dominant cell type",fill="Cell type")

  p_ent <- ggplot(spots, aes_string(x="x",y="y",fill="entropy")) +
    geom_tile(color="white",linewidth=.1) +
    scale_fill_viridis_c(option="magma",name="Entropy") +
    coord_equal() + theme_void(base_size=9) +
    theme(legend.position="right") + labs(title="B3. Spot composition entropy (heterogeneity)")

  # neighborhood-averaged CD8 fraction
  spots$CD8_nb <- nb[,"CD8 T"]
  p_cd8 <- ggplot(spots, aes_string(x="x",y="y",fill="CD8_nb")) +
    geom_tile(color="white",linewidth=.1) +
    scale_fill_gradient(low="#f8f9fa",high="#d62828",name="CD8 frac") +
    coord_equal() + theme_void(base_size=9) +
    theme(legend.position="right") + labs(title="B4. Neighborhood-smoothed CD8 T fraction")

  patchwork::wrap_plots(p_region, p_dom, p_ent, p_cd8, nrow=2) +
    patchwork::plot_annotation(title="scalefusion: Spatial composition across zoom levels",
                               subtitle="Region → dominant type → entropy → smoothed marker: each reveals different biology")
}

# ── Plot C: Cross-scale disagreement ─────────────────────────────────────────
scalefusion_plot_disagreement <- function() {
  d <- scalefusion_expr_data()

  # Compute eta-squared (proportion of total variance explained) at each scale
  # eta² = SS_between / SS_total = 1 means perfectly explained by grouping
  eta_sq <- function(expr, group) {
    grand_mean <- mean(expr)
    ss_total <- sum((expr - grand_mean)^2)
    if (ss_total < 1e-12) return(0)
    gm <- tapply(expr, group, mean)
    gn <- tapply(expr, group, length)
    ss_between <- sum(gn * (gm - grand_mean)^2)
    ss_between / ss_total
  }

  genes_all <- unique(d$gene)
  vr <- data.frame(gene = genes_all, stringsAsFactors = FALSE)
  vr$eta_ct <- sapply(genes_all, function(g) {
    sub <- d[d$gene == g, ]; eta_sq(sub$expr, sub$celltype)
  })
  vr$eta_comp <- sapply(genes_all, function(g) {
    sub <- d[d$gene == g, ]; eta_sq(sub$expr, sub$compartment)
  })
  vr$info_lost_ct <- 1 - vr$eta_ct
  vr$info_lost_comp <- 1 - vr$eta_comp
  vr$gene <- factor(vr$gene, levels = vr$gene[order(vr$info_lost_comp)])

  vr_long <- vr %>%
    dplyr::select(gene, info_lost_ct, info_lost_comp) %>%
    pivot_longer(-gene, names_to = "level", values_to = "frac_lost")
  vr_long$level <- ifelse(vr_long$level == "info_lost_ct",
                           "Cell \u2192 CellType", "Cell \u2192 Compartment")
  vr_long$level <- factor(vr_long$level, levels = c("Cell \u2192 CellType", "Cell \u2192 Compartment"))

  p_dis <- ggplot(vr_long, aes(x = gene, y = frac_lost, fill = level)) +
    geom_col(position = position_dodge(width = .6), width = .5) +
    geom_hline(yintercept = .5, linetype = "dashed", color = "#d00000", linewidth = .5) +
    scale_fill_manual(values = c("Cell \u2192 CellType" = "#2a9d8f",
                                  "Cell \u2192 Compartment" = "#e76f51")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    coord_flip() +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "#666666")) +
    labs(x = NULL, y = "Information lost by aggregation",
         fill = NULL,
         title = "Cross-scale disagreement: information lost by aggregation",
         subtitle = "Dashed line = 50% information loss; genes above this are scale-sensitive")

  p_dis + patchwork::plot_annotation(
    title = "scalefusion: Quantifying semantic drift across resolution levels",
    subtitle = "Some genes are robust to aggregation, others lose their discriminative signal",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "#444444")
    )
  )
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
scalefusion_plot_demo <- function() { scalefusion_plot_expr() }

scalefusion_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ggplot2::ggsave(file.path(output_dir,"scalefusion_01_expr_ladder.png"),
                  plot=scalefusion_plot_expr(), width=12,height=14,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"scalefusion_02_spatial.png"),
                  plot=scalefusion_plot_spatial(), width=13,height=10,dpi=300,bg="white")
}

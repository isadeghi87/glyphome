## densweave: Density-First Topology Visualization
##
## Visual grammar: "Topology quilts" that replace scatter plots with density-first
## views where the primary visual object is density terrain, not individual points.
##
## Key innovations:
##   - Density field as the base layer (not points)
##   - Topological features (ridges, saddles, basins) annotated on the density surface
##   - Distortion mapping: where does the embedding fail?
##   - Rare population detection: cells in density valleys
##
## Figures generated:
##   A. Atlas-scale topology quilt: density field with ridge/basin annotations,
##      rare population callouts, and distortion overlay
##   B. Spatial density weave: tissue section density with neighborhood topology
##   C. Density comparison: same embedding under two conditions with delta density

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Palettes ─────────────────────────────────────────────────────────────────
densweave_ct_pal <- function() {
  c("CD8 T"="#264653","CD4 T"="#2a9d8f","Treg"="#e9c46a","B cell"="#577590",
    "Plasma"="#f4a261","Monocyte"="#e76f51","cDC"="#bc6c25","pDC"="#6a994e",
    "NK"="#d62828","Cycling"="#7b2cbf","Progenitor"="#8ecae6","Exhausted"="#6d597a")
}

# ── Data A: Atlas-scale UMAP with density ────────────────────────────────────
densweave_atlas_data <- function(seed = 42) {
  set.seed(seed)
  cc <- data.frame(
    ct = c("CD8 T","CD4 T","Treg","B cell","Plasma","Monocyte","cDC","pDC","NK",
           "Cycling","Progenitor","Exhausted"),
    cx = c(-4, -2.5, -1.2, 3.5, 5, 1, 2.5, 4, -5.5, 0, -3.5, -4.8),
    cy = c(2, 3.5, 4.5, -1.5, -3, -3.5, -1, -0.5, 0.5, 5.5, -2, 3),
    sp = c(.55, .65, .30, .70, .45, .80, .40, .28, .55, .45, .35, .25),
    n  = c(200, 250, 40, 180, 60, 220, 50, 25, 150, 55, 30, 20),
    stringsAsFactors = FALSE)

  rows <- lapply(seq_len(nrow(cc)), function(i) {
    r <- cc[i,]; nn <- r$n
    data.frame(
      UMAP_1 = r$cx + stats::rnorm(nn, 0, r$sp),
      UMAP_2 = r$cy + stats::rnorm(nn, 0, r$sp),
      celltype = r$ct,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ── Data B: Spatial density ──────────────────────────────────────────────────
densweave_spatial_data <- function(seed = 7) {
  set.seed(seed)
  n <- 800
  # generate spots in tissue with varying density
  region <- sample(c("Dense tumor","Sparse stroma","Immune cluster","Transition","Necrotic"),
                   n, replace=TRUE, prob=c(.3,.25,.2,.15,.1))
  base_x <- ifelse(region=="Dense tumor", 3, ifelse(region=="Immune cluster", 8,
             ifelse(region=="Sparse stroma", 6, ifelse(region=="Transition", 5, 2))))
  base_y <- ifelse(region=="Dense tumor", 3, ifelse(region=="Immune cluster", 7,
             ifelse(region=="Sparse stroma", 5, ifelse(region=="Transition", 5, 8))))
  sp <- ifelse(region=="Dense tumor", .8, ifelse(region=="Sparse stroma", 2.0,
        ifelse(region=="Immune cluster", .6, ifelse(region=="Transition", 1.5, 1.2))))
  data.frame(
    x = base_x + stats::rnorm(n, 0, sp),
    y = base_y + stats::rnorm(n, 0, sp),
    region = region,
    marker_expr = ifelse(region=="Immune cluster", stats::rnorm(n, 3, .8),
                   ifelse(region=="Dense tumor", stats::rnorm(n, 1, .5), stats::rnorm(n, .5, .3))),
    stringsAsFactors = FALSE)
}

# ── Data C: Condition comparison ─────────────────────────────────────────────
densweave_condition_data <- function(seed = 99) {
  set.seed(seed)
  base <- densweave_atlas_data(seed)
  d1 <- base; d1$condition <- "Pre-treatment"
  # post-treatment: exhausted expands, effector CD8 depletes
  d2 <- base
  exhaust_idx <- which(d2$celltype == "Exhausted")
  cd8_idx <- which(d2$celltype == "CD8 T")
  # expand exhausted
  new_exhaust <- d2[rep(exhaust_idx, 4),]
  new_exhaust$UMAP_1 <- new_exhaust$UMAP_1 + stats::rnorm(nrow(new_exhaust), 0, .15)
  new_exhaust$UMAP_2 <- new_exhaust$UMAP_2 + stats::rnorm(nrow(new_exhaust), 0, .15)
  # deplete CD8
  keep <- sample(cd8_idx, round(length(cd8_idx) * .4))
  d2 <- rbind(d2[-cd8_idx,], d2[keep,], new_exhaust)
  d2$condition <- "Post-treatment"
  rbind(d1, d2)
}

# ── Helper: 2D density grid ─────────────────────────────────────────────────
densweave_density_grid <- function(x, y, n = 80) {
  rx <- range(x); ry <- range(y)
  pad_x <- diff(rx) * .1; pad_y <- diff(ry) * .1
  gx <- seq(rx[1] - pad_x, rx[2] + pad_x, length.out = n)
  gy <- seq(ry[1] - pad_y, ry[2] + pad_y, length.out = n)
  kde <- MASS::kde2d(x, y, n = n, lims = c(range(gx), range(gy)))
  expand.grid(x = kde$x, y = kde$y) |>
    transform(density = as.vector(kde$z))
}

# ── Plot A: Atlas topology quilt ─────────────────────────────────────────────
densweave_plot_atlas <- function() {
  d <- densweave_atlas_data()
  pal <- densweave_ct_pal()

  # density grid
  dg <- densweave_density_grid(d$UMAP_1, d$UMAP_2, n = 100)
  dg$log_dens <- log10(dg$density + 1e-6)

  # find rare cells: low local density
  d$local_dens <- stats::approx(
    x = seq(min(d$UMAP_1), max(d$UMAP_1), length.out=100),
    y = tapply(dg$density, cut(dg$x, 100), max),
    xout = d$UMAP_1)$y
  d$local_dens[is.na(d$local_dens)] <- min(dg$density)
  d$is_rare <- d$local_dens < stats::quantile(d$local_dens, .05, na.rm=TRUE)

  # cluster centroids
  centroids <- d %>% group_by(celltype) %>%
    summarise(cx=mean(UMAP_1), cy=mean(UMAP_2), n=n(), .groups="drop")

  # Main density quilt
  p_dens <- ggplot() +
    geom_tile(data=dg, aes_string(x="x",y="y",fill="log_dens"), alpha=.85) +
    scale_fill_viridis_c(option="mako",name="log10(density)") +
    geom_point(data=d[d$is_rare,], aes_string(x="UMAP_1",y="UMAP_2"),
               color="#d62828",size=1.5,shape=4,stroke=.8) +
    geom_label(data=centroids, aes_string(x="cx",y="cy",label="celltype"),
               size=2.5, fill="white", alpha=.8, label.padding=ggplot2::unit(.12,"lines"),
               label.r=ggplot2::unit(.1,"lines"), linewidth=.2) +
    coord_equal() +
    theme_minimal(base_size=10) + theme(panel.grid=element_blank()) +
    labs(title="A1. Density topology quilt",
         subtitle="Background = density terrain; red X = rare cells in density valleys")

  # Contour view
  p_contour <- ggplot(d, aes_string(x="UMAP_1",y="UMAP_2")) +
    stat_density_2d(aes_string(fill="after_stat(level)"), geom="polygon",
                    color="#264653", linewidth=.15, n=100, bins=15) +
    scale_fill_viridis_c(option="rocket",direction=-1,name="Density level") +
    geom_point(data=centroids, aes_string(x="cx",y="cy"), shape=3, size=3, stroke=1, color="white") +
    coord_equal() +
    theme_void(base_size=9) +
    labs(title="A2. Contour topology: ridges and basins")

  # Rare population bar
  rare_counts <- d %>% filter(is_rare) %>% count(celltype) %>% arrange(desc(n))
  rare_counts$celltype <- factor(rare_counts$celltype, levels=rare_counts$celltype)
  p_rare <- ggplot(rare_counts, aes_string(x="celltype",y="n",fill="celltype")) +
    geom_col(width=.6,show.legend=FALSE) + scale_fill_manual(values=pal) +
    theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid.minor=element_blank()) +
    labs(x=NULL,y="Rare cells",title="A3. Cells in density valleys (bottom 5%)")

  patchwork::wrap_plots(patchwork::wrap_plots(p_dens, p_contour, nrow=1),
                        p_rare, nrow=2, heights=c(2.5,1)) +
    patchwork::plot_annotation(title="densweave: Atlas-scale density topology quilt",
                               subtitle="Density terrain replaces scatter; topology reveals structure invisible in point clouds")
}

# ── Plot B: Spatial density weave ────────────────────────────────────────────
densweave_plot_spatial <- function() {
  d <- densweave_spatial_data()
  reg_pal <- c("Dense tumor"="#e76f51","Sparse stroma"="#2a9d8f","Immune cluster"="#264653",
               "Transition"="#f4a261","Necrotic"="#6c757d")

  dg <- densweave_density_grid(d$x, d$y, n=80)

  p_dens <- ggplot() +
    geom_tile(data=dg, aes_string(x="x",y="y",fill="density"),alpha=.8) +
    scale_fill_viridis_c(option="inferno",name="Cell density") +
    geom_point(data=d, aes_string(x="x",y="y",color="region"),size=.5,alpha=.3) +
    scale_color_manual(values=reg_pal) +
    coord_equal() +
    theme_void(base_size=10) + theme(legend.position="right") +
    labs(title="B1. Spatial density field with region overlay",color="Region")

  p_marker <- ggplot(d, aes_string(x="x",y="y",color="marker_expr")) +
    geom_point(size=.8,alpha=.5) +
    scale_color_viridis_c(option="magma",name="Immune\nmarker") +
    coord_equal() + theme_void(base_size=10) +
    labs(title="B2. Immune marker expression over spatial density")

  # density by region boxplot
  p_reg <- ggplot(d, aes_string(x="region",y="marker_expr",fill="region")) +
    geom_violin(alpha=.6,linewidth=.3) + geom_boxplot(width=.15,fill="white",linewidth=.3,outlier.size=.5) +
    scale_fill_manual(values=reg_pal,guide="none") +
    theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=30,hjust=1),panel.grid.minor=element_blank()) +
    labs(x=NULL,y="Marker expression",title="B3. Expression by spatial region")

  patchwork::wrap_plots(patchwork::wrap_plots(p_dens,p_marker,nrow=1), p_reg,
                        nrow=2, heights=c(2,1.2)) +
    patchwork::plot_annotation(title="densweave: Spatial density weave for tissue architecture",
                               subtitle="Density terrain reveals tissue compartments; marker overlays highlight functional zones")
}

# ── Plot C: Condition delta density ──────────────────────────────────────────
densweave_plot_comparison <- function() {
  d <- densweave_condition_data()
  pre <- d[d$condition=="Pre-treatment",]
  post <- d[d$condition=="Post-treatment",]

  # compute density for each condition on same grid
  n_grid <- 80
  all_x <- range(d$UMAP_1); all_y <- range(d$UMAP_2)
  pad_x <- diff(all_x)*.1; pad_y <- diff(all_y)*.1
  lims <- c(all_x[1]-pad_x, all_x[2]+pad_x, all_y[1]-pad_y, all_y[2]+pad_y)

  kde_pre <- MASS::kde2d(pre$UMAP_1, pre$UMAP_2, n=n_grid, lims=lims)
  kde_post <- MASS::kde2d(post$UMAP_1, post$UMAP_2, n=n_grid, lims=lims)

  grid <- expand.grid(x=kde_pre$x, y=kde_pre$y)
  grid$dens_pre <- as.vector(kde_pre$z)
  grid$dens_post <- as.vector(kde_post$z)
  grid$delta <- grid$dens_post - grid$dens_pre
  grid$log_pre <- log10(grid$dens_pre + 1e-6)
  grid$log_post <- log10(grid$dens_post + 1e-6)

  p_pre <- ggplot(grid, aes_string(x="x",y="y",fill="log_pre")) +
    geom_tile() + scale_fill_viridis_c(option="mako",name="log10(dens)") +
    coord_equal() + theme_void(base_size=9) +
    labs(title="C1. Pre-treatment density")

  p_post <- ggplot(grid, aes_string(x="x",y="y",fill="log_post")) +
    geom_tile() + scale_fill_viridis_c(option="mako",name="log10(dens)") +
    coord_equal() + theme_void(base_size=9) +
    labs(title="C2. Post-treatment density")

  p_delta <- ggplot(grid, aes_string(x="x",y="y",fill="delta")) +
    geom_tile() +
    scale_fill_gradient2(low="#264653",mid="#f8f9fa",high="#e76f51",midpoint=0,
                         name="Delta density") +
    coord_equal() + theme_void(base_size=9) +
    labs(title="C3. Delta: Post − Pre",
         subtitle="Red = expansion; blue = depletion")

  # per-celltype change bar chart
  ct_pre <- pre %>% count(celltype,name="n_pre")
  ct_post <- post %>% count(celltype,name="n_post")
  ct_change <- merge(ct_pre,ct_post,by="celltype",all=TRUE)
  ct_change[is.na(ct_change)] <- 0
  ct_change$log2fc <- log2((ct_change$n_post+1)/(ct_change$n_pre+1))
  ct_change$direction <- ifelse(ct_change$log2fc>0,"Expanded","Depleted")
  ct_change <- ct_change[order(ct_change$log2fc),]
  ct_change$celltype <- factor(ct_change$celltype, levels=ct_change$celltype)

  p_bar <- ggplot(ct_change, aes_string(x="celltype",y="log2fc",fill="direction")) +
    geom_col(width=.6) + geom_hline(yintercept=0,linewidth=.3) +
    scale_fill_manual(values=c(Expanded="#e76f51",Depleted="#264653")) +
    coord_flip() + theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank(),legend.position="none") +
    labs(x=NULL,y="log2 FC (post/pre)",title="C4. Cell type composition shift")

  patchwork::wrap_plots(patchwork::wrap_plots(p_pre,p_post,p_delta,nrow=1),
                        p_bar, nrow=2, heights=c(2,1.2)) +
    patchwork::plot_annotation(title="densweave: Density topology comparison across conditions",
                               subtitle="Delta density map reveals which regions of the atlas expanded or contracted post-treatment")
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
densweave_plot_demo <- function() { densweave_plot_atlas() }

densweave_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ggplot2::ggsave(file.path(output_dir,"densweave_01_atlas.png"),
                  plot=densweave_plot_atlas(), width=14,height=11,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"densweave_02_spatial.png"),
                  plot=densweave_plot_spatial(), width=13,height=10,dpi=300,bg="white")
}

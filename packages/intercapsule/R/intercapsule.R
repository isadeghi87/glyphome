## intercapsule: Structured Interaction Capsule Visualization
##
## Visual grammar: "Interaction capsules" — compound glyphs where each cell-cell
## or gene-gene interaction is a structured object encoding:
##   Sender → Receiver → Channel/Context → Confidence → Evidence → Sensitivity
##
## Unlike standard network edges, capsules are 6-dimensional visual objects that
## can be inspected, filtered, and compared across conditions.
##
## Figures generated:
##   A. Spatial ligand-receptor interaction map: capsule glyphs between cell types
##      with sender/receiver roles, channels, and confidence
##   B. Host-microbiome interaction capsules: structured interactions between
##      immune cells and microbial taxa
##   C. Interaction comparison: same capsule set across two conditions (tumor vs normal)

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Palettes ─────────────────────────────────────────────────────────────────
intercapsule_ct_pal <- function() {
  c("CD8 T"="#264653","CD4 T"="#2a9d8f","Treg"="#e9c46a","Macrophage"="#e76f51",
    "cDC"="#577590","Tumor"="#d62828","Fibroblast"="#bc6c25","B cell"="#6a994e",
    "NK"="#7b2cbf","Endothelial"="#f4a261")
}
intercapsule_channel_pal <- function() {
  c("PD1-PDL1"="#264653","CD28-B7"="#2a9d8f","CXCL9-CXCR3"="#e9c46a",
    "CCL5-CCR5"="#e76f51","FAS-FASL"="#d62828","TNF-TNFR"="#577590",
    "IL10-IL10R"="#6d597a","TGFB1-TGFBR"="#bc6c25","VEGFA-FLT1"="#f4a261",
    "WNT5A-FZD5"="#7b2cbf","CD40-CD40L"="#8ecae6","CSF1-CSF1R"="#6a994e",
    "NOTCH1-DLL4"="#fb8500","HGF-MET"="#90be6d")
}

# ── Data A: Spatial LR interactions ──────────────────────────────────────────
intercapsule_spatial_data <- function(seed = 42) {
  set.seed(seed)
  senders <- c("Tumor","Tumor","Macrophage","cDC","Fibroblast","CD8 T","CD4 T","NK",
               "Tumor","Macrophage","Fibroblast","Endothelial","Treg","Tumor","B cell")
  receivers <- c("CD8 T","Treg","CD4 T","CD8 T","Tumor","Tumor","Macrophage","Tumor",
                 "Fibroblast","Fibroblast","Endothelial","Tumor","CD8 T","NK","cDC")
  channels <- c("PD1-PDL1","TGFB1-TGFBR","CD40-CD40L","CXCL9-CXCR3","HGF-MET",
                "FAS-FASL","TNF-TNFR","FAS-FASL","VEGFA-FLT1","CSF1-CSF1R",
                "NOTCH1-DLL4","VEGFA-FLT1","IL10-IL10R","CCL5-CCR5","CD40-CD40L")
  n <- length(senders)
  data.frame(
    sender=senders, receiver=receivers, channel=channels,
    confidence=pmin(1,pmax(.1, stats::rbeta(n,5,3))),
    evidence=sample(c("Co-expression","Proximity","Literature","All three"), n, replace=TRUE,
                    prob=c(.3,.25,.15,.30)),
    sensitivity=pmin(1,pmax(0, stats::rbeta(n,3,5))),
    effect_size=stats::rnorm(n, 0, 1.5),
    context=sample(c("Tumor margin","Immune hub","Stroma","Core"), n, replace=TRUE),
    stringsAsFactors=FALSE)
}

# ── Data B: Host-microbiome interactions ─────────────────────────────────────
intercapsule_microbiome_data <- function(seed = 7) {
  set.seed(seed)
  hosts <- c("CD8 T","CD4 T","Treg","Macrophage","cDC","NK","B cell","CD8 T",
             "Macrophage","CD4 T","Treg","NK","cDC","B cell","Macrophage")
  microbes <- c("F.nucleatum","B.fragilis","C.difficile","E.coli","L.rhamnosus",
                "A.muciniphila","B.longum","F.nucleatum","C.difficile","L.rhamnosus",
                "B.fragilis","E.coli","A.muciniphila","B.longum","F.nucleatum")
  mechanism <- c("TLR4","PSA","TcdB","LPS","SLP","Mucin","SCFA","TLR4","TcdB",
                 "SLP","PSA","LPS","Mucin","SCFA","TLR4")
  n <- length(hosts)
  data.frame(
    sender=microbes, receiver=hosts, channel=mechanism,
    confidence=pmin(1,pmax(.1, stats::rbeta(n,4,4))),
    evidence=sample(c("16S","Metagenomics","Culture","Multi-omic"), n, replace=TRUE),
    sensitivity=pmin(1,pmax(0, stats::rbeta(n,2,4))),
    effect=c(1,-1,-1,1,1,1,1,1,-1,1,-1,1,1,1,1),  # pro/anti-inflammatory
    stringsAsFactors=FALSE)
}

# ── Data C: Condition comparison ─────────────────────────────────────────────
intercapsule_condition_data <- function(seed = 99) {
  set.seed(seed)
  base <- intercapsule_spatial_data(seed)
  d1 <- base; d1$condition <- "Tumor"
  d2 <- base; d2$condition <- "Normal"
  # in normal tissue: immunosuppressive interactions weaker
  suppress_ch <- c("PD1-PDL1","TGFB1-TGFBR","IL10-IL10R")
  d2$confidence[d2$channel %in% suppress_ch] <- d2$confidence[d2$channel %in% suppress_ch] * .3
  d2$effect_size[d2$channel %in% suppress_ch] <- d2$effect_size[d2$channel %in% suppress_ch] * .4
  # activating interactions stronger in normal
  d2$confidence[!d2$channel %in% suppress_ch] <- pmin(1, d2$confidence[!d2$channel %in% suppress_ch] * 1.3)
  rbind(d1, d2)
}

# ── Helper: capsule glyph layout ─────────────────────────────────────────────
# Position senders on left column, receivers on right column; arcs between them
capsule_layout <- function(d) {
  senders_u <- unique(d$sender); receivers_u <- unique(d$receiver)
  ns <- length(senders_u); nr <- length(receivers_u)
  sx <- rep(0, ns); sy <- seq(ns, 1, length.out=ns)
  rx <- rep(4, nr); ry <- seq(nr, 1, length.out=nr)
  d$sx <- sx[match(d$sender, senders_u)]; d$sy <- sy[match(d$sender, senders_u)]
  d$rx <- rx[match(d$receiver, receivers_u)]; d$ry <- ry[match(d$receiver, receivers_u)]
  d$mx <- (d$sx + d$rx) / 2; d$my <- (d$sy + d$ry) / 2
  list(d=d, senders=data.frame(name=senders_u, x=sx, y=sy, stringsAsFactors=FALSE),
       receivers=data.frame(name=receivers_u, x=rx, y=ry, stringsAsFactors=FALSE))
}

# ── Plot A: Spatial LR capsules ──────────────────────────────────────────────
intercapsule_plot_spatial <- function() {
  d <- intercapsule_spatial_data()
  ct_pal <- intercapsule_ct_pal(); ch_pal <- intercapsule_channel_pal()
  lay <- capsule_layout(d); d <- lay$d; sn <- lay$senders; rc <- lay$receivers

  sn$fill <- unname(ct_pal[sn$name]); rc$fill <- unname(ct_pal[rc$name])
  d$ch_fill <- unname(ch_pal[d$channel])

  p <- ggplot() +
    # arcs: confidence = width, channel = color
    geom_curve(data=d, aes_string(x="sx",y="sy",xend="rx",yend="ry",
                                   linewidth="confidence",color="ch_fill",alpha="1-sensitivity"),
               curvature=.2, show.legend=FALSE) +
    scale_color_identity() + scale_alpha_identity() +
    scale_linewidth_continuous(range=c(.3, 3)) +
    # evidence marks at midpoint
    geom_point(data=d, aes_string(x="mx",y="my",shape="evidence"),
               size=2, fill="white", color="#343a40", stroke=.5) +
    scale_shape_manual(values=c("Co-expression"=21,"Proximity"=22,"Literature"=23,"All three"=24)) +
    # sender nodes
    geom_point(data=sn, aes_string(x="x",y="y"), fill=sn$fill, shape=21, size=6, stroke=.8,color="#2b2d42") +
    geom_text(data=sn, aes_string(x="x",y="y",label="name"), hjust=1.3, size=3) +
    # receiver nodes
    geom_point(data=rc, aes_string(x="x",y="y"), fill=rc$fill, shape=21, size=6, stroke=.8,color="#2b2d42") +
    geom_text(data=rc, aes_string(x="x",y="y",label="name"), hjust=-.3, size=3) +
    annotate("text",x=0,y=max(sn$y)+.8,label="SENDERS",fontface="bold",size=3.5,hjust=.5) +
    annotate("text",x=4,y=max(rc$y)+.8,label="RECEIVERS",fontface="bold",size=3.5,hjust=.5) +
    coord_cartesian(xlim=c(-2.5, 6.5)) +
    theme_void(base_size=10) +
    theme(legend.position="bottom") +
    labs(shape="Evidence type",
         title="A. Spatial ligand-receptor interaction capsules",
         subtitle="Arc width=confidence; opacity=1−sensitivity; midpoint shape=evidence type")

  p + patchwork::plot_annotation(
    title="intercapsule: 6-dimensional interaction capsules",
    subtitle="Sender → Receiver → Channel → Confidence → Evidence → Sensitivity in one glyph")
}

# ── Plot B: Host-microbiome ──────────────────────────────────────────────────
intercapsule_plot_microbiome <- function() {
  d <- intercapsule_microbiome_data()
  ct_pal <- intercapsule_ct_pal()
  lay <- capsule_layout(d); d <- lay$d; sn <- lay$senders; rc <- lay$receivers
  d$eff_col <- ifelse(d$effect > 0, "#2a9d8f", "#e76f51")
  rc$fill <- unname(ct_pal[rc$name])

  p <- ggplot() +
    geom_curve(data=d, aes_string(x="sx",y="sy",xend="rx",yend="ry",
                                   linewidth="confidence",color="eff_col",alpha="1-sensitivity"),
               curvature=.15, show.legend=FALSE) +
    scale_color_identity() + scale_alpha_identity() +
    scale_linewidth_continuous(range=c(.3,2.5)) +
    geom_point(data=d, aes_string(x="mx",y="my",shape="evidence"),
               size=2, fill="white", color="#343a40", stroke=.5) +
    scale_shape_manual(values=c("16S"=21,"Metagenomics"=22,"Culture"=23,"Multi-omic"=24)) +
    geom_point(data=sn, aes_string(x="x",y="y"), fill="#bc6c25", shape=23, size=5, stroke=.8, color="#2b2d42") +
    geom_text(data=sn, aes_string(x="x",y="y",label="name"), hjust=1.3, size=2.8, fontface="italic") +
    geom_point(data=rc, aes_string(x="x",y="y"), fill=rc$fill, shape=21, size=5, stroke=.8, color="#2b2d42") +
    geom_text(data=rc, aes_string(x="x",y="y",label="name"), hjust=-.3, size=2.8) +
    annotate("text",x=0,y=max(sn$y)+.8,label="MICROBES",fontface="bold",size=3.5) +
    annotate("text",x=4,y=max(rc$y)+.8,label="HOST IMMUNE",fontface="bold",size=3.5) +
    coord_cartesian(xlim=c(-2.5,6.5)) +
    theme_void(base_size=10) + theme(legend.position="bottom") +
    labs(shape="Evidence",
         title="B. Host-microbiome interaction capsules",
         subtitle="Green=pro-inflammatory, Orange=anti-inflammatory; diamond=microbial taxa")

  p + patchwork::plot_annotation(
    title="intercapsule: Microbiome-immune cross-kingdom interaction capsules",
    subtitle="Structured glyphs reveal mechanism, direction, and evidence quality")
}

# ── Plot C: Condition comparison ─────────────────────────────────────────────
intercapsule_plot_comparison <- function() {
  d <- intercapsule_condition_data()
  ch_pal <- intercapsule_channel_pal()

  # heatmap: confidence per interaction per condition
  d$interaction <- paste0(d$sender, " → ", d$receiver, "\n(", d$channel, ")")
  heat_d <- d %>% dplyr::select(interaction, condition, confidence)
  heat_w <- heat_d %>% pivot_wider(names_from=condition, values_from=confidence)
  heat_w$delta <- heat_w$Tumor - heat_w$Normal
  heat_w$interaction <- factor(heat_w$interaction, levels=heat_w$interaction[order(heat_w$delta)])

  p_heat <- ggplot(heat_d, aes_string(x="condition",y="interaction",fill="confidence")) +
    geom_tile(color="white",linewidth=.5) +
    scale_fill_gradient(low="#f8f9fa",high="#264653",limits=c(0,1),name="Confidence") +
    theme_minimal(base_size=9) +
    theme(panel.grid=element_blank(), axis.text.y=element_text(size=7)) +
    labs(x=NULL,y=NULL,title="C1. Capsule confidence: Tumor vs Normal")

  # delta bar chart
  heat_w_long <- heat_w[order(heat_w$delta),]
  heat_w_long$direction <- ifelse(heat_w_long$delta > 0, "Tumor-enriched", "Normal-enriched")
  p_delta <- ggplot(heat_w_long, aes_string(x="delta",y="interaction",fill="direction")) +
    geom_col(width=.6) +
    geom_vline(xintercept=0, linewidth=.4, color="#adb5bd") +
    scale_fill_manual(values=c("Tumor-enriched"="#e76f51","Normal-enriched"="#2a9d8f")) +
    theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank(),axis.text.y=element_blank()) +
    labs(x="Delta confidence (Tumor − Normal)",y=NULL,fill=NULL,title="C2. Condition delta")

  patchwork::wrap_plots(p_heat, p_delta, nrow=1, widths=c(2,1.5)) +
    patchwork::plot_annotation(title="intercapsule: Interaction capsule comparison across conditions",
                               subtitle="Immunosuppressive capsules amplified in tumor; activating capsules stronger in normal tissue")
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
intercapsule_plot_demo <- function() { intercapsule_plot_spatial() }

intercapsule_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ggplot2::ggsave(file.path(output_dir,"intercapsule_01_spatial_lr.png"),
                  plot=intercapsule_plot_spatial(), width=12,height=10,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"intercapsule_02_microbiome.png"),
                  plot=intercapsule_plot_microbiome(), width=12,height=10,dpi=300,bg="white")
}

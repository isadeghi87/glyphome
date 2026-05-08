## provtrace: Provenance-Native Bioinformatics Visualization
##
## Visual grammar: Each result is a "provenance card" showing not just WHAT the
## result is but HOW it was produced—pipeline variant, parameter sensitivity,
## reproducibility across methods, and agreement topology.
##
## Figures generated:
##   A. DE provenance: Sankey-like flow of genes across pipeline variants with
##      consensus heatmap and agreement topology ribbon
##   B. Clustering provenance: Alluvial cluster assignment across methods with
##      parameter sensitivity and stability scores
##   C. Method agreement network: Graph showing which methods co-identify results

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggraph)
library(igraph)

# ── Palettes ─────────────────────────────────────────────────────────────────
provtrace_method_pal <- function() {
  c("DESeq2"="#264653","edgeR"="#2a9d8f","limma-voom"="#e9c46a","Wilcoxon"="#f4a261",
    "MAST"="#e76f51","Seurat_SCT"="#577590","scran"="#bc6c25")
}
provtrace_stability_pal <- function() scales::col_numeric(c("#f8f9fa","#2a9d8f","#264653"),domain=c(0,1))

# ── Data: DE provenance ──────────────────────────────────────────────────────
provtrace_de_data <- function(seed = 42) {
  set.seed(seed)
  genes <- paste0("Gene_", sprintf("%03d", 1:80))
  methods <- c("DESeq2","edgeR","limma-voom","Wilcoxon","MAST")
  nm <- length(methods)
  # each gene has a "true" signal strength; methods with different sensitivity
  true_lfc <- stats::rnorm(80, 0, 2)
  sens <- c(1, .9, .85, .6, .95)  # relative sensitivity
  spec <- c(.92, .88, .90, .70, .85) # specificity
  rows <- list()
  for (j in seq_along(methods)) {
    detect <- abs(true_lfc)*sens[j] + stats::rnorm(80,0,.3) > 1.0
    detect <- detect & (stats::runif(80) < spec[j] | abs(true_lfc) > 1.5)
    lfc_obs <- true_lfc * sens[j] + stats::rnorm(80, 0, .25)
    pv <- 10^(-(abs(lfc_obs)*2.5 + stats::runif(80,-.5,.5)))
    pv[!detect] <- pv[!detect] * 10^stats::runif(sum(!detect), 1, 3)
    rows[[j]] <- data.frame(gene=genes, method=methods[j], logFC=lfc_obs,
                            pvalue=pmin(1,pv), detected=detect, stringsAsFactors=FALSE)
  }
  d <- do.call(rbind, rows)
  # add consensus count per gene
  cons <- aggregate(detected ~ gene, d, sum)
  names(cons)[2] <- "n_methods"
  merge(d, cons, by="gene")
}

# ── Data: Clustering provenance ──────────────────────────────────────────────
provtrace_cluster_data <- function(seed = 7) {
  set.seed(seed)
  methods <- c("Louvain","Leiden_0.5","Leiden_1.0","Leiden_1.5","SC3","HDBSCAN","Seurat_SCT")
  nc <- c(8, 6, 10, 14, 9, 7, 11)  # cluster counts per method
  cells <- paste0("Cell_", sprintf("%04d", 1:200))
  rows <- list()
  # base cluster assignment
  base <- sample(1:8, 200, replace=TRUE, prob=c(.15,.12,.18,.10,.10,.15,.12,.08))
  for (j in seq_along(methods)) {
    cl <- base
    # add method-specific noise
    flip_rate <- c(.00, .08, .05, .12, .10, .15, .07)[j]
    flip <- stats::runif(200) < flip_rate
    cl[flip] <- sample(1:nc[j], sum(flip), replace=TRUE)
    cl <- as.integer(factor(cl))  # renumber
    rows[[j]] <- data.frame(cell=cells, method=methods[j], cluster=paste0("C",cl), stringsAsFactors=FALSE)
  }
  do.call(rbind, rows)
}

# ── Plot A: DE provenance ────────────────────────────────────────────────────
provtrace_plot_de <- function() {
  d <- provtrace_de_data()

  # consensus heatmap: top 30 genes by max |logFC|
  top30 <- d %>% group_by(gene) %>% summarise(maxlfc = max(abs(logFC))) %>%
    arrange(desc(maxlfc)) %>% head(30) %>% pull(gene)
  dsub <- d[d$gene %in% top30, ]
  dsub$gene <- factor(dsub$gene, levels = rev(top30))
  dsub$sig <- ifelse(dsub$pvalue < .05, "sig", "ns")

  p_heat <- ggplot(dsub, aes(x = method, y = gene, fill = logFC)) +
    geom_tile(color = "white", linewidth = .3) +
    geom_point(data = dsub[dsub$sig == "ns", ], aes(x = method, y = gene),
               shape = 4, size = 1.5, color = "#6c757d", inherit.aes = FALSE) +
    scale_fill_gradient2(low = "#264653", mid = "#f8f9fa", high = "#e76f51", midpoint = 0, name = "logFC") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          panel.grid = element_blank(),
          legend.key.height = unit(.5, "cm"),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "#666666")) +
    labs(x = NULL, y = NULL, title = "A. Consensus DE heatmap (top 30 genes)",
         subtitle = "\u2715 = not significant (padj > 0.05) in that method")

  # agreement ribbon: number of methods detecting each gene
  gene_cons <- d %>% filter(gene %in% top30) %>% group_by(gene) %>%
    summarise(n_detect = sum(detected), mean_lfc = mean(logFC[detected]))
  gene_cons$gene <- factor(gene_cons$gene, levels = rev(top30))
  gene_cons$fill_col <- ifelse(gene_cons$n_detect == 5, "#264653",
                                ifelse(gene_cons$n_detect >= 4, "#2a9d8f",
                                       ifelse(gene_cons$n_detect >= 3, "#e9c46a", "#e76f51")))
  p_ribbon <- ggplot(gene_cons, aes(x = n_detect, y = gene)) +
    geom_col(aes(fill = fill_col), width = .7) +
    geom_text(aes(label = n_detect), size = 2.8, hjust = -.3, fontface = "bold") +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0, 6.2), breaks = 0:5) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = "# Methods", y = NULL, title = "Agreement")

  # logFC concordance with R² annotation
  w <- d %>% filter(method %in% c("DESeq2","MAST")) %>%
    dplyr::select(gene, method, logFC) %>% pivot_wider(names_from = method, values_from = logFC)
  r2 <- round(cor(w$DESeq2, w$MAST, use = "complete.obs")^2, 3)
  p_conc <- ggplot(w, aes(x = DESeq2, y = MAST)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#adb5bd") +
    geom_point(size = 1.8, alpha = .55, color = "#264653") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#e76f51", se = TRUE, alpha = .12, linewidth = .7) +
    annotate("text", x = min(w$DESeq2, na.rm = TRUE) + .5,
             y = max(w$MAST, na.rm = TRUE) - .3,
             label = paste0("R\u00b2 = ", r2), size = 4, color = "#264653", fontface = "bold") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(title = "B. logFC concordance: DESeq2 vs MAST", x = "DESeq2 logFC", y = "MAST logFC")

  patchwork::wrap_plots(
    patchwork::wrap_plots(p_heat, p_ribbon, nrow = 1, widths = c(4, 1)),
    p_conc, nrow = 2, heights = c(3, 1.5)
  ) + patchwork::plot_annotation(
    title = "provtrace: Differential expression provenance",
    subtitle = "Only genes reproducible across methods are trustworthy",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "#444444")
    )
  )
}

# ── Plot B: Clustering provenance ────────────────────────────────────────────
provtrace_plot_clustering <- function() {
  d <- provtrace_cluster_data()
  methods_ordered <- c("Louvain","Leiden_0.5","Leiden_1.0","Leiden_1.5","SC3","HDBSCAN","Seurat_SCT")

  # per-cell stability: fraction of method pairs that agree
  cell_stab <- d %>% group_by(cell) %>%
    summarise(n_unique = length(unique(cluster)),
              most_common = names(sort(table(cluster), decreasing = TRUE))[1],
              stability = max(table(cluster)) / length(cluster))

  # alluvial-like: method vs cluster as stacked bars
  freq <- d %>% group_by(method, cluster) %>% summarise(n = n(), .groups = "drop")
  freq$method <- factor(freq$method, levels = methods_ordered)
  n_cl <- length(unique(d$cluster))
  cl_pal <- scales::hue_pal()(max(n_cl, 20))[seq_len(n_cl)]
  names(cl_pal) <- sort(unique(d$cluster))
  p_alluvial <- ggplot(freq, aes(x = method, y = n, fill = cluster)) +
    geom_col(position = "stack", width = .7, color = "white", linewidth = .2) +
    scale_fill_manual(values = cl_pal) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none", panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "#666666")) +
    labs(x = NULL, y = "Cells", title = "A. Cluster composition across methods",
         subtitle = "Color = cluster label; methods reorder and split clusters differently")

  # stability density
  med_stab <- median(cell_stab$stability)
  p_stab <- ggplot(cell_stab, aes(x = stability)) +
    geom_histogram(fill = "#264653", color = "white", bins = 25, linewidth = .3) +
    geom_vline(xintercept = .5, linetype = "dashed", color = "#e76f51", linewidth = .7) +
    annotate("text", x = .52, y = Inf, vjust = 1.5, label = "Threshold = 0.5",
             size = 3, color = "#e76f51", fontface = "italic") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = "Stability score", y = "Cells", title = "B. Cell-level stability",
         subtitle = "Fraction of methods assigning most common cluster")

  # resolution sensitivity: number of clusters by method
  res_df <- d %>% group_by(method) %>% summarise(n_clusters = length(unique(cluster)))
  res_df$method <- factor(res_df$method, levels = methods_ordered)
  p_res <- ggplot(res_df, aes(x = method, y = n_clusters)) +
    geom_col(fill = "#2a9d8f", width = .6) +
    geom_text(aes(label = n_clusters), vjust = -.3, size = 3.5, fontface = "bold") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 11)) +
    labs(x = NULL, y = "# Clusters", title = "C. Resolution sensitivity")

  patchwork::wrap_plots(
    p_alluvial,
    patchwork::wrap_plots(p_stab, p_res, nrow = 1),
    nrow = 2, heights = c(2, 1.2)
  ) + patchwork::plot_annotation(
    title = "provtrace: Clustering provenance across methods and resolutions",
    subtitle = "Stability < 0.5 (dashed line): cells for which no majority assignment exists",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "#444444")
    )
  )
}

# ── Plot C: Method agreement network ──────────────────────────────────────────
provtrace_plot_agreement_net <- function() {
  d <- provtrace_de_data()
  methods <- unique(d$method); nm <- length(methods)
  det_mat <- d %>% dplyr::select(gene, method, detected) %>%
    pivot_wider(names_from = method, values_from = detected)
  mat <- as.matrix(det_mat[, -1])
  agree <- matrix(0, nm, nm, dimnames = list(methods, methods))
  for (i in 1:(nm - 1)) for (j in (i + 1):nm) {
    agree[i, j] <- sum(mat[, i] & mat[, j]) / sum(mat[, i] | mat[, j])
    agree[j, i] <- agree[i, j]
  }
  diag(agree) <- 1
  el <- which(upper.tri(agree), arr.ind = TRUE)
  edges <- data.frame(from = methods[el[, 1]], to = methods[el[, 2]],
                       jaccard = agree[el], stringsAsFactors = FALSE)
  edges <- edges[edges$jaccard > .25, ]

  g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = methods))
  V(g)$n_detected <- sapply(methods, function(m) sum(d$detected[d$method == m]))

  p_net <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = jaccard, alpha = jaccard), color = "#264653", show.legend = TRUE) +
    geom_node_point(aes(size = n_detected), fill = "#2a9d8f", color = "#264653", shape = 21, stroke = .8) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
    scale_edge_width(range = c(.5, 4), name = "Jaccard") +
    scale_edge_alpha(range = c(.3, .9), guide = "none") +
    scale_size_continuous(range = c(5, 14), name = "Genes detected") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, color = "#555555"),
          legend.position = "right") +
    labs(title = "Method agreement network",
         subtitle = "Edge width = Jaccard index of shared DE genes; node size = total detected")

  p_net + patchwork::plot_annotation(
    title = "provtrace: Pipeline method agreement topology",
    subtitle = "Methods that cluster together report similar DE gene sets",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "#444444")
    )
  )
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
provtrace_plot_demo <- function() { provtrace_plot_de() }

provtrace_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ggplot2::ggsave(file.path(output_dir,"provtrace_01_de_provenance.png"),
                  plot=provtrace_plot_de(), width=13,height=11,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"provtrace_03_agreement_net.png"),
                  plot=provtrace_plot_agreement_net(), width=10,height=9,dpi=300,bg="white")
}

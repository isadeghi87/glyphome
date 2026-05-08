## glyphome-common: shared theme, palette helpers, and figure-saving utilities
##
## All ten visualization modules in this workspace consume the same base
## theme and figure-saving conventions defined here. Sourcing this file gives
## downstream code a consistent typographic/visual language and avoids
## duplicating ggplot theme code in every package.
##
## Usage:
##   source("R/glyphome-common.R")
##   p + glyphome_theme()
##   glyphome_save(p, "fig.png", role = "primary")

# ── Base theme ───────────────────────────────────────────────────────────────

#' glyphome base theme
#'
#' A minimal, publication-grade ggplot2 theme used by every module in the
#' workspace. Scales sensibly for composite figures.
#'
#' @param base_size base font size in pt (default 10)
#' @param base_family font family (default "sans")
#' @return a ggplot2 theme object
glyphome_theme <- function(base_size = 10, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = ggplot2::element_text(color = "#444444",
                                                  size = base_size - 1,
                                                  margin = ggplot2::margin(b = 6)),
      plot.caption       = ggplot2::element_text(color = "#666666",
                                                  size = base_size - 2, hjust = 0),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_line(color = "#e9ecef", linewidth = 0.25),
      strip.text         = ggplot2::element_text(face = "bold"),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.key.size    = ggplot2::unit(0.7, "lines")
    )
}

#' Void theme variant for image-like panels (spatial maps, embeddings)
glyphome_theme_void <- function(base_size = 10) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = "#444444", size = base_size - 1)
    )
}

# ── Figure-size table by panel role ─────────────────────────────────────────
# Aligned with the figure-composition standards: hero panels get more area,
# auxiliary panels get less.
glyphome_size <- function(role = c("primary", "hero", "supporting", "auxiliary")) {
  role <- match.arg(role)
  switch(role,
    hero       = list(width = 16, height = 12),
    primary    = list(width = 13, height = 10),
    supporting = list(width = 10, height = 8),
    auxiliary  = list(width = 8,  height = 6)
  )
}

#' Save a publication-grade figure with role-aware sizing
#'
#' @param plot a ggplot or patchwork object
#' @param path destination path
#' @param role one of "hero", "primary", "supporting", "auxiliary"
#' @param dpi output resolution (default 300)
#' @param ... passed to ggplot2::ggsave
glyphome_save <- function(plot, path, role = "primary", dpi = 300, ...) {
  s <- glyphome_size(role)
  ggplot2::ggsave(path, plot = plot,
                  width = s$width, height = s$height,
                  dpi = dpi, bg = "white", ...)
  invisible(path)
}

# ── Module registry ──────────────────────────────────────────────────────────
# Single source of truth for which modules exist, what role they play, and
# which figures they generate. Used by the master figure generator and the
# tutorial book.
glyphome_modules <- function() {
  data.frame(
    module = c("uncertainr","provtrace","scalefusion","counterpath","intercapsule",
               "densweave","clinladder","survivora","cohortprism","riskweave"),
    family = c("uncertainty","provenance","scale","dynamics","interaction",
               "topology","clinical","clinical","clinical","clinical"),
    one_liner = c(
      "Uncertainty as a first-class visual element",
      "Provenance and pipeline dependency inside the figure",
      "Cross-scale semantic continuity across resolution levels",
      "Temporal and counterfactual trajectory fans",
      "Interaction-first capsule glyphs for multimodal biology",
      "Density-first topology quilts replacing scatter",
      "Decision-oriented evidence ladders for clinical genomics",
      "Prognostic landscape cartography (hazard as terrain)",
      "Spectral patient prisms with continuous + discrete + temporal",
      "Multi-objective treatment decision frontier"
    ),
    n_figures = c(3,2,2,3,2,2,3,3,3,3),
    stringsAsFactors = FALSE
  )
}

#' Source all module files into the calling environment
#'
#' @param base_dir workspace root (default: detected from getwd)
glyphome_load_all <- function(base_dir = getwd()) {
  pkgs <- glyphome_modules()$module
  for (pkg in pkgs) {
    src <- file.path(base_dir, "packages", pkg, "R", paste0(pkg, ".R"))
    if (file.exists(src)) source(src, chdir = TRUE)
  }
  invisible(pkgs)
}

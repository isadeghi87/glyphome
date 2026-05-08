#!/usr/bin/env Rscript
## Master figure generator for the glyphome suite.
##
## Loads shared utilities, then sources every module and runs its
## <module>_save_demo() entry point. PNGs are written to `figures/`.
##
## Usage:
##   Rscript scripts/generate_all_figures.R
##   Rscript scripts/generate_all_figures.R --module densweave intercapsule

suppressPackageStartupMessages({
  required <- c("ggplot2","dplyr","tidyr","patchwork","scales",
                "ggraph","igraph","MASS","ggnewscale","grid")
  for (p in required) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("Missing required package: ", p)
    library(p, character.only = TRUE)
  }
})

base_dir <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = FALSE),
  error = function(e) getwd())
if (!nzchar(base_dir) || !dir.exists(base_dir)) base_dir <- getwd()

source(file.path(base_dir, "R", "glyphome-common.R"))
source(file.path(base_dir, "R", "glyphome-render.R"))

fig_dir <- file.path(base_dir, "figures")

# Optional CLI: --module <name1> <name2> ...
args <- commandArgs(trailingOnly = TRUE)
modules <- NULL
if (length(args) >= 2 && args[1] == "--module") {
  modules <- args[-1]
}

cat("=== glyphome: rendering figure set ===\n")
cat(sprintf("Workspace: %s\n", base_dir))
cat(sprintf("Output:    %s\n", fig_dir))
if (!is.null(modules)) cat(sprintf("Modules:   %s\n", paste(modules, collapse = ", ")))
cat("\n")

ok <- glyphome_render_all(fig_dir, base_dir = base_dir, modules = modules)

n_png <- length(list.files(fig_dir, pattern = "\\.png$"))
cat(sprintf("\n=== Done. %d/%d modules rendered, %d PNGs in %s ===\n",
            length(ok), nrow(glyphome_modules()), n_png, fig_dir))

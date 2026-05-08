# AGENTS.md — visualization_package (glyphome v0.2)

## Project Overview
The workspace is the **glyphome** umbrella R package + ten stand-alone
visualization sub-packages for publication-grade bioinformatics figures.
Each sub-package targets a specific gap in the conventional figure vocabulary
(uncertainty, provenance, scale, dynamics, interaction, topology, clinical
decision support).

## Architecture
- `DESCRIPTION` / `NAMESPACE` / `R/` — glyphome umbrella package (shared theme,
  size table, `glyphome_render_all()` driver, module registry).
- `packages/<module>/` — ten stand-alone R packages, each with its own
  `DESCRIPTION`, auto-generated `NAMESPACE`, and a single monolithic
  `R/<module>.R` containing palettes, data generators, plot constructors,
  and a `<module>_save_demo()` entry point.
- `scripts/generate_all_figures.R` — CLI driver (supports `--module`).
- `figures/` — 26 curated PNGs (4 redundant figures removed in v0.2).
- `book/` — complete bookdown tutorial and reference manual (16 chapters).

## Figure Generation Skill Pack
Palettes and visualization patterns from this package have been abstracted into the
**Figure Generation Skill Pack** at: `/home/isadeghi/projects/figure_generation_skill/`

Key palettes are documented in:
`figure_generation_skill/instructions/color_and_accessibility.md`

## Packages

| Package | File | Purpose | Key Functions |
|---|---|---|---|
| densweave | `packages/densweave/R/densweave.R` | Atlas-scale UMAP density, spatial density | `densweave_ct_pal()`, `densweave_atlas_data()`, `densweave_spatial_data()` |
| intercapsule | `packages/intercapsule/R/intercapsule.R` | LR interaction capsule glyphs | `intercapsule_ct_pal()`, `intercapsule_channel_pal()` |
| scalefusion | `packages/scalefusion/R/scalefusion.R` | Multi-resolution scale ladders | `scalefusion_spatial_pal()`, `scalefusion_ct_pal()` |
| uncertainr | `packages/uncertainr/R/uncertainr.R` | Uncertainty-native UMAPs | `uncertainr_celltype_palette()`, `uncertainr_unc_scales()` |
| clinladder | `packages/clinladder/R/clinladder.R` | Clinical genomics review panels | `clinladder_*` |
| provtrace | `packages/provtrace/R/provtrace.R` | Provenance-aware DE visualization | `provtrace_*` |
| counterpath | `packages/counterpath/R/counterpath.R` | Temporal/counterfactual visualization | `counterpath_*` |
| survivora | `packages/survivora/R/survivora.R` | Prognostic landscape cartography (HR terrain + stability contours) | `survivora_plot_landscape()`, `survivora_plot_locality_km()`, `survivora_plot_repro_forest()` |
| cohortprism | `packages/cohortprism/R/cohortprism.R` | Spectral patient prisms — continuous signatures + mutations + clinical events | `cohortprism_plot_wall()`, `cohortprism_plot_single()`, `cohortprism_plot_projection()` |
| riskweave | `packages/riskweave/R/riskweave.R` | Multi-objective treatment decision frontier (efficacy × toxicity × confidence) | `riskweave_plot_frontier()`, `riskweave_plot_contributions()`, `riskweave_plot_sensitivity()` |

## Usage
Each R file is self-contained. Source directly:
```r
source("/home/isadeghi/projects/visualization_package/packages/densweave/R/densweave.R")
pal <- densweave_ct_pal()
```

## Palette Reference
All palettes are named character vectors (cell type name → hex color).
`densweave_ct_pal()` is the primary 30-color atlas palette recommended for UMAP plots.
`intercapsule_ct_pal()` is recommended for LR interaction visualizations.

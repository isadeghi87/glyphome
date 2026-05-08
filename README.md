# glyphome — advanced bioinformatics visualization suite

[![Version](https://img.shields.io/badge/version-0.2.0-2a9d8f)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/license-MIT-264653)](LICENSE)
[![Modules](https://img.shields.io/badge/modules-10-e76f51)](packages/)
[![Figures](https://img.shields.io/badge/curated_figures-26-e9c46a)](figures/)

`glyphome` is a meta-package that bundles **ten focused R visualization
modules**, each addressing a distinct gap in the conventional bioinformatics
figure vocabulary — from uncertainty as a first-class channel to
multi-objective treatment decision frontiers. The suite ships with a
complete [bookdown tutorial and reference manual](book/) plus a uniform
render-all pipeline.

---

## Installation

### Requirements

- R >= 4.1.0
- A C/C++ toolchain (`build-essential` on Debian/Ubuntu, Xcode CLT on macOS,
  Rtools on Windows) — needed by `igraph` and `MASS`.

### One-line install (umbrella + all 10 modules)

```r
install.packages("remotes")

remotes::install_github("isadeghi87/glyphome")               # umbrella

# Install all ten sub-packages from their subdirectories:
mods <- c("uncertainr","provtrace","scalefusion","counterpath",
          "intercapsule","densweave","clinladder","survivora",
          "cohortprism","riskweave")
for (m in mods)
  remotes::install_github("isadeghi87/glyphome",
                          subdir = paste0("packages/", m))
```

### Install a single module

If you only need one figure family (e.g. density-first topology), install
just that module — the sub-packages have no cross-dependency:

```r
remotes::install_github("isadeghi87/glyphome", subdir = "packages/densweave")
library(densweave)
densweave_plot_atlas()
```

### Development install from a local clone

```bash
git clone git@github.com:isadeghi87/glyphome.git
cd glyphome

# system deps once (Ubuntu/Debian example)
sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev \
                        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
                        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

```r
# from R, in the cloned directory
install.packages(c("remotes","devtools","bookdown"))

# install all CRAN dependencies first
remotes::install_deps(dependencies = TRUE)

# install the umbrella package
devtools::install(".")

# install all sub-packages
for (m in list.dirs("packages", recursive = FALSE))
  devtools::install(m)
```

### Conda / mamba environment (optional)

```bash
mamba create -n glyphome -c conda-forge \
  r-base r-ggplot2 r-patchwork r-dplyr r-tidyr r-scales \
  r-mass r-ggraph r-igraph r-ggnewscale r-remotes r-bookdown
mamba activate glyphome
Rscript -e 'remotes::install_github("isadeghi87/glyphome")'
```

## Quick tour

### Render the full curated figure set

```r
library(glyphome)
glyphome_load_all()                       # source every module
glyphome_render_all("figures")            # writes 26 PNGs
```

or from a shell:

```bash
Rscript scripts/generate_all_figures.R
Rscript scripts/generate_all_figures.R --module densweave intercapsule
```

### A single figure in five lines

```r
library(densweave)
library(glyphome)
p <- densweave_plot_atlas()
glyphome_save(p, "atlas.png", role = "hero")
```

### Use the shared theme on your own ggplots

```r
library(ggplot2); library(glyphome)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  glyphome_theme() +
  labs(title = "Title", subtitle = "Subtitle in muted grey")
```

### Build the tutorial locally

```r
install.packages("bookdown")
bookdown::render_book("book/")
browseURL("book/_book/index.html")
```

### Tutorial structure

The tutorial is written as a full package reference manual, not just a gallery.
It starts with installation and the common grammar, then moves through the ten
modules and ends with composition rules for publication-grade multi-panel
figures:

1. Installation, quick tour, and rendering workflow.
2. Shared grammar: themes, palettes, role-based sizing, and saving.
3. One chapter per module with concept, API, runnable examples, and design notes.
4. Composite-figure assembly and QC.
5. How to add an eleventh module.

## The ten modules

| Module | Family | What it does |
|---|---|---|
| [uncertainr](book/03-uncertainr.Rmd) | uncertainty | Decomposed uncertainty channels on every glyph |
| [provtrace](book/04-provtrace.Rmd)   | provenance  | Pipeline-fragility traces inside the figure |
| [scalefusion](book/05-scalefusion.Rmd) | scale     | Cross-scale semantic continuity |
| [counterpath](book/06-counterpath.Rmd) | dynamics  | Counterfactual trajectory fans |
| [intercapsule](book/07-intercapsule.Rmd) | interaction | Multi-channel capsule glyphs for L–R interactions |
| [densweave](book/08-densweave.Rmd)   | topology    | Density terrain replacing scatter |
| [clinladder](book/09-clinladder.Rmd) | clinical    | Decision-oriented evidence ladders |
| [survivora](book/10-survivora.Rmd)   | clinical    | Hazard as terrain, stability as contour |
| [cohortprism](book/11-cohortprism.Rmd) | clinical  | Patient prisms (continuous + discrete + temporal) |
| [riskweave](book/12-riskweave.Rmd)   | clinical    | Multi-objective treatment decision frontier |

## Architecture

```
glyphome/                              ← glyphome umbrella package
├── DESCRIPTION  NAMESPACE  LICENSE
├── R/                                 # shared theme, sizes, render-all
│   ├── glyphome-package.R
│   ├── glyphome-common.R
│   └── glyphome-render.R
├── packages/<module>/                 # 10 stand-alone R packages
│   ├── DESCRIPTION
│   ├── NAMESPACE                      # auto-generated from <module>_*
│   └── R/<module>.R
├── scripts/
│   └── generate_all_figures.R         # CLI driver
├── figures/                           # 26 curated PNGs
└── book/                              # bookdown tutorial/reference (16 chapters)
```

Every module follows the same contract:

- one or more `<module>_*_data()` data generators,
- one or more `<module>_plot_*()` plot constructors,
- a `<module>_save_demo(output_dir)` entry point.

This contract is what allows `glyphome_render_all()` to drive the entire
suite with a few lines of code, and what allows the bookdown reference to
stay in sync with the module set.

## Tutorial and reference manual

Build the book locally:

```r
bookdown::render_book("book/")
browseURL("book/_book/index.html")
```

The tutorial covers:

- **Foundations.** Architecture, shared visual grammar, role-based panel sizing.
- **Modules.** One chapter per module: concept, figures, API, worked example,
  design notes.
- **Composing publication figures.** Layout patterns, the QC checklist, the
  refinement loop, when to split a figure.
- **Extending glyphome.** How to author and register an eleventh module.

## What changed in v0.2

- Promoted to a true umbrella package (`glyphome`) with shared theme,
  size table, and `glyphome_render_all()` driver.
- Auto-generated `NAMESPACE` for every sub-package; standardised
  `DESCRIPTION` metadata.
- Removed four redundant figures (kept 26 of the strongest 30):
  `densweave_03_comparison`, `intercapsule_03_comparison`,
  `scalefusion_03_disagreement`, `provtrace_02_clustering`.
- Added the bookdown tutorial/reference manual under `book/` (16 chapters).
- New CLI flags for the master generator (`--module`).

## License

MIT — see [LICENSE](LICENSE).

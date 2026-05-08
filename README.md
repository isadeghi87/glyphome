# bioviz — advanced bioinformatics visualization suite

[![Version](https://img.shields.io/badge/version-0.2.0-2a9d8f)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/license-MIT-264653)](LICENSE)
[![Modules](https://img.shields.io/badge/modules-10-e76f51)](packages/)
[![Figures](https://img.shields.io/badge/curated_figures-26-e9c46a)](figures/)

`bioviz` is a meta-package that bundles **ten focused R visualization
modules**, each addressing a distinct gap in the conventional bioinformatics
figure vocabulary — from uncertainty as a first-class channel to
multi-objective treatment decision frontiers. The suite ships with a
[ComplexHeatmap-style reference book](book/) (`bookdown`) and a uniform
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

remotes::install_github("isadeghi87/bioviz")               # umbrella

# Install all ten sub-packages from their subdirectories:
mods <- c("uncertainr","provtrace","scalefusion","counterpath",
          "intercapsule","densweave","clinladder","survivora",
          "cohortprism","riskweave")
for (m in mods)
  remotes::install_github("isadeghi87/bioviz",
                          subdir = paste0("packages/", m))
```

### Install a single module

If you only need one figure family (e.g. density-first topology), install
just that module — the sub-packages have no cross-dependency:

```r
remotes::install_github("isadeghi87/bioviz", subdir = "packages/densweave")
library(densweave)
densweave_plot_atlas()
```

### Development install from a local clone

```bash
git clone git@github.com:isadeghi87/bioviz.git
cd bioviz

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
mamba create -n bioviz -c conda-forge \
  r-base r-ggplot2 r-patchwork r-dplyr r-tidyr r-scales \
  r-mass r-ggraph r-igraph r-ggnewscale r-remotes r-bookdown
mamba activate bioviz
Rscript -e 'remotes::install_github("isadeghi87/bioviz")'
```

## Quick tour

### Render the full curated figure set

```r
library(bioviz)
bioviz_load_all()                       # source every module
bioviz_render_all("figures")            # writes 26 PNGs
```

or from a shell:

```bash
Rscript scripts/generate_all_figures.R
Rscript scripts/generate_all_figures.R --module densweave intercapsule
```

### A single figure in five lines

```r
library(densweave)
library(bioviz)
p <- densweave_plot_atlas()
bioviz_save(p, "atlas.png", role = "hero")
```

### Use the shared theme on your own ggplots

```r
library(ggplot2); library(bioviz)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  bioviz_theme() +
  labs(title = "Title", subtitle = "Subtitle in muted grey")
```

### Build the reference book locally

```r
install.packages("bookdown")
bookdown::render_book("book/")
browseURL("book/_book/index.html")
```

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
visualization_package/                  ← bioviz umbrella package
├── DESCRIPTION  NAMESPACE  LICENSE
├── R/                                 # shared theme, sizes, render-all
│   ├── bioviz-package.R
│   ├── bioviz-common.R
│   └── bioviz-render.R
├── packages/<module>/                 # 10 stand-alone R packages
│   ├── DESCRIPTION
│   ├── NAMESPACE                      # auto-generated from <module>_*
│   └── R/<module>.R
├── scripts/
│   └── generate_all_figures.R         # CLI driver
├── figures/                           # 26 curated PNGs
└── book/                              # bookdown reference (15 chapters)
```

Every module follows the same contract:

- one or more `<module>_*_data()` data generators,
- one or more `<module>_plot_*()` plot constructors,
- a `<module>_save_demo(output_dir)` entry point.

This contract is what allows `bioviz_render_all()` to drive the entire
suite with a few lines of code, and what allows the bookdown reference to
stay in sync with the module set.

## Reference book

Build the book locally:

```r
bookdown::render_book("book/")
browseURL("book/_book/index.html")
```

The book covers:

- **Foundations.** Architecture, shared visual grammar, role-based panel sizing.
- **Modules.** One chapter per module: concept, figures, API, worked example,
  design notes.
- **Composing publication figures.** Layout patterns, the QC checklist, the
  refinement loop, when to split a figure.
- **Extending bioviz.** How to author and register an eleventh module.

## What changed in v0.2

- Promoted to a true umbrella package (`bioviz`) with shared theme,
  size table, and `bioviz_render_all()` driver.
- Auto-generated `NAMESPACE` for every sub-package; standardised
  `DESCRIPTION` metadata.
- Removed four redundant figures (kept 26 of the strongest 30):
  `densweave_03_comparison`, `intercapsule_03_comparison`,
  `scalefusion_03_disagreement`, `provtrace_02_clustering`.
- Added the bookdown reference manual under `book/` (15 chapters).
- New CLI flags for the master generator (`--module`).

## License

MIT — see [LICENSE](LICENSE).

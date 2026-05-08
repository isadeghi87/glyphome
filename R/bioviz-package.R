#' bioviz: Advanced Bioinformatics Visualization Suite
#'
#' bioviz is a meta-package that bundles ten focused visualization modules,
#' each addressing a distinct gap in the conventional bioinformatics figure
#' vocabulary. Modules are independently usable but share a common visual
#' language (theme, size table, save helper) defined in this package.
#'
#' @section Module families:
#'  \itemize{
#'    \item \strong{Uncertainty} — \code{uncertainr}
#'    \item \strong{Provenance} — \code{provtrace}
#'    \item \strong{Scale}      — \code{scalefusion}
#'    \item \strong{Dynamics}   — \code{counterpath}, \code{intercapsule}
#'    \item \strong{Topology}   — \code{densweave}
#'    \item \strong{Clinical}   — \code{clinladder}, \code{survivora},
#'                                \code{cohortprism}, \code{riskweave}
#'  }
#'
#' @section Quick start:
#' \preformatted{
#'   library(bioviz)
#'   bioviz_load_all()
#'   bioviz_render_all("figures")
#' }
#'
#' @section Tutorial:
#' A long-form tutorial in the style of the ComplexHeatmap reference book
#' is provided under \code{book/}. Build with
#' \code{bookdown::render_book("book/")}.
#'
#' @docType package
#' @name bioviz-package
#' @aliases bioviz
NULL

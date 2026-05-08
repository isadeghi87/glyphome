## glyphome: render-all driver
##
## Top-level entry point that loads every module and renders the curated
## figure set (one to three panels per module). Used by both the master
## generator script and the tutorial book.

#' Render the full glyphome figure set
#'
#' Loads every module in [glyphome_modules()] and runs its `<module>_save_demo()`
#' entry point, writing PNGs to `output_dir`.
#'
#' @param output_dir directory to write figures into; created if needed.
#' @param base_dir   workspace root containing `packages/`
#' @param modules    optional character vector to render only a subset
#' @return invisible character vector of module names that succeeded.
glyphome_render_all <- function(output_dir,
                              base_dir = getwd(),
                              modules = NULL) {
  reg <- glyphome_modules()
  if (!is.null(modules)) reg <- reg[reg$module %in% modules, , drop = FALSE]
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ok <- character(0)
  for (i in seq_len(nrow(reg))) {
    m <- reg$module[i]
    src <- file.path(base_dir, "packages", m, "R", paste0(m, ".R"))
    if (!file.exists(src)) {
      message(sprintf("[%s] source file missing: %s", m, src)); next
    }
    env <- new.env()
    tryCatch({
      sys.source(src, envir = env)
      fn <- get(paste0(m, "_save_demo"), envir = env)
      fn(output_dir)
      ok <- c(ok, m)
      message(sprintf("[%s] rendered", m))
    }, error = function(e) {
      message(sprintf("[%s] FAILED: %s", m, conditionMessage(e)))
    })
  }
  invisible(ok)
}

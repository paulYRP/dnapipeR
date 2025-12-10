#' Generate DNAm PDF Report from DNAm.Rmd - inst/scripts/.R
#'
#' @export
dnamReport <- function(
    output = "DNAm_Report.pdf",
    outputDir = "reports",

    qcDir = "figures/preprocessingMinfiEwasWater/enMix",
    preprocessingDir = "figures/preprocessingMinfiEwasWater/qc",
    postprocessingDir = "figures/preprocessingMinfiEwasWater/metrics",
    svaDir = "figures/svaEnmix/sva",
    glmDir = "figures/methylationGLM_T1",
    glmmDir = "figures/methylationGLMM_T1T2",
    reportTitle = "DNA methylation",
    author = "School of Biomedical Sciences",
    date = format(Sys.Date(), "%B %d, %Y")
) {

  # Helper: Convert to absolute path
  makeAbs <- function(path) {
    if (grepl("^([A-Za-z]:|/)", path)) return(path)
    return(file.path(getwd(), path))
  }

  # Locate Rmd + Runner Script
  rmd <- system.file("scripts", "DNAm.Rmd", package = "dnapipeR")
  script <- system.file("scripts", "DNAm.R", package = "dnapipeR")

  if (rmd == "" || script == "")
    stop("DNAm.Rmd or DNAm.R not found in package.")

  outDirAbs <- makeAbs(outputDir)
  dir.create(outDirAbs, recursive = TRUE, showWarnings = FALSE)

  figDir <- file.path(outputDir, "figures")
  dir.create(figDir, recursive = TRUE, showWarnings = FALSE)


  # Build argument list for external Rscript
  arg_list <- c(
    "--rmd",        shQuote(rmd),
    "--output",     shQuote(output),
    "--outputDir",  shQuote(outDirAbs),
    "--qcDir",      shQuote(makeAbs(qcDir)),
    "--preDir",     shQuote(makeAbs(preprocessingDir)),
    "--postDir",    shQuote(makeAbs(postprocessingDir)),
    "--svaDir",     shQuote(makeAbs(svaDir)),
    "--glmDir",     shQuote(makeAbs(glmDir)),
    "--glmmDir",    shQuote(makeAbs(glmmDir)),
    "--title",      shQuote(reportTitle),
    "--author",     shQuote(author),
    "--date",       shQuote(date),
    "--figDir", shQuote(figDir)
  )

  # Construct Rscript command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))

  # Run command
  message("Running dnamReport():")
  message(cmd)

  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}

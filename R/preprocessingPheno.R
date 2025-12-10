#' Run preprocessingPheno.R as an external script from inst/script/.R
#'
#'
#' @export
preprocessingPheno <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    sepType = "",
    betaPath = "rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    mPath = "rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    cnPath = "rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    SampleID = "Sample_Name",
    timeVar = "Timepoint",
    timepoints = "1,2",
    combineTimepoints = "1,2",
    outputPheno = "data/preprocessingPheno",
    outputRData = "rData/preprocessingPheno/metrics",
    outputRDataMerge = "rData/preprocessingPheno/mergeData",
    sexColumn = "Sex",
    outputLogs = "logs",
    outputDir = "data/preprocessingPheno"
) {
  
  # Locate script inside installed package
  script <- system.file("scripts", "preprocessingPheno.R", package = "dnapipeR")
  if (script == "")
    stop("Script preprocessingPheno.R not found in package.")
  
  # Build argument list (shQuote for safety)
  arg_list <- c(
    "--phenoFile", shQuote(phenoFile),
    "--sepType", shQuote(sepType),
    "--betaPath", shQuote(betaPath),
    "--mPath", shQuote(mPath),
    "--cnPath", shQuote(cnPath),
    "--SampleID", shQuote(SampleID),
    "--timeVar", shQuote(timeVar),
    "--timepoints", shQuote(timepoints),
    "--combineTimepoints", shQuote(combineTimepoints),
    "--outputPheno", shQuote(outputPheno),
    "--outputRData", shQuote(outputRData),
    "--outputRDataMerge", shQuote(outputRDataMerge),
    "--sexColumn", shQuote(sexColumn),
    "--outputLogs", shQuote(outputLogs),
    "--outputDir", shQuote(outputDir)
  )
  
  # Build full command for printing
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))
  
  # User-visible messages
  message("Running preprocessingPheno:")
  message(cmd)
  
  # Silent execution (Windows + Linux/HPC)
  invisible(
    if (.Platform$OS.type == "windows") {
      shell(cmd, intern = FALSE, translate = TRUE)
    } else {
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
  )
}

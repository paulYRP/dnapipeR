#' Run svaEnmix.R as an external script from inst/script/.R
#'
#'
#' @export
svaEnmix <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    rgsetData = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData",
    sepType = "",
    outputLogs = "logs",
    nSamples = NA,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    SentrixIDColumn = "Sentrix_ID",
    SentrixPositionColumn = "Sentrix_Position",
    ctrlSvaPercVar = 0.90,
    ctrlSvaFlag = 1,
    scriptLabel = "svaEnmix",
    tiffWidth = 2000,
    tiffHeight = 1000,
    tiffRes = 150
) {
  
  # Locate script inside installed package
  script <- system.file("scripts", "svaEnmix.R", package = "dnapipeR")
  if (script == "")
    stop("Script svaEnmix.R not found in package.")
  
  # Build argument list
  arg_list <- c(
    "--phenoFile", shQuote(phenoFile),
    "--rgsetData", shQuote(rgsetData),
    "--sepType", shQuote(sepType),
    "--outputLogs", shQuote(outputLogs),
    "--nSamples", nSamples,
    "--SampleID", shQuote(SampleID),
    "--arrayType", shQuote(arrayType),
    "--annotationVersion", shQuote(annotationVersion),
    "--SentrixIDColumn", shQuote(SentrixIDColumn),
    "--SentrixPositionColumn", shQuote(SentrixPositionColumn),
    "--ctrlSvaPercVar", ctrlSvaPercVar,
    "--ctrlSvaFlag", ctrlSvaFlag,
    "--scriptLabel", shQuote(scriptLabel),
    "--tiffWidth", tiffWidth,
    "--tiffHeight", tiffHeight,
    "--tiffRes", tiffRes
  )
  
  # Build system command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))
  
  message("Running svaEnmix:")
  message(cmd)
  
  # Platform specific execution
  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}

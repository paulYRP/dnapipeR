#' Run preprocessingMinfiEwasWater.R as an external script from inst/script/.R
#'
#'
#' @export
preprocessingMinfiEwasWater <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/pheno.csv",
    idatFolder = "data/preprocessingMinfiEwasWater/idats",
    outputLogs = "logs",
    nSamples = NA,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    scriptLabel = "preprocessingMinfiEwasWater",
    baseDataFolder = "rData",
    sepType = "",
    tiffWidth = 2000,
    tiffHeight = 1000,
    tiffRes = 150,
    qcCutoff = 10.5,
    detPtype = "m+u",
    detPThreshold = 0.05,
    funnormSeed = 123,
    normMethods = "adjustedfunnorm",
    sexColumn = "Sex",
    pvalThreshold = 0.01,
    chrToRemove = "chrX,chrY",
    snpsToRemove = "SBE,CpG",
    mafThreshold = 0.1,
    crossReactivePath = "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv",
    plotGroupVar = "Sex",
    lcRef = "salivaEPIC",
    phenoOrder = "Sample_Name;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position",
    lcPhenoDir = "data/preprocessingMinfiEwasWater"
) {

  # Resolve script inside the package
  script <- system.file("scripts", "preprocessingMinfiEwasWater.R", package = "dnapipeR")

  if (script == "") stop("Script preprocessingMinfiEwasWater.R not found in package.")

  # Build command-line arguments
  arg_list <- c(
    "--phenoFile", shQuote(phenoFile),
    "--sepType", shQuote(sepType),
    "--idatFolder", shQuote(idatFolder),
    "--outputLogs", shQuote(outputLogs),
    "--nSamples", nSamples,
    "--SampleID", shQuote(SampleID),
    "--arrayType", shQuote(arrayType),
    "--annotationVersion", shQuote(annotationVersion),
    "--scriptLabel", shQuote(scriptLabel),
    "--baseDataFolder", shQuote(baseDataFolder),
    "--tiffWidth", tiffWidth,
    "--tiffHeight", tiffHeight,
    "--tiffRes", tiffRes,
    "--qcCutoff", qcCutoff,
    "--detPtype", shQuote(detPtype),
    "--detPThreshold", detPThreshold,
    "--funnormSeed", funnormSeed,
    "--normMethods", shQuote(normMethods),
    "--sexColumn", shQuote(sexColumn),
    "--pvalThreshold", pvalThreshold,
    "--chrToRemove", shQuote(chrToRemove),
    "--snpsToRemove", shQuote(snpsToRemove),
    "--mafThreshold", mafThreshold,
    "--crossReactivePath", shQuote(crossReactivePath),
    "--plotGroupVar", shQuote(plotGroupVar),
    "--lcRef", shQuote(lcRef),
    "--phenoOrder", shQuote(phenoOrder),
    "--lcPhenoDir", shQuote(lcPhenoDir)
  )

  # Build system command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))

  message("Running preprocessingMinfiEwasWater:")
  message(cmd)

  # Platform specific execution
  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}

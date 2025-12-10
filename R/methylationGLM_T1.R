#' Run methylationGLM_T1.R as an external script from inst/script/.R
#'
#'
#' @export
methylationGLM_T1 <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1.RData",
    outputLogs = "logs",
    outputRData = "rData/methylationGLM_T1/models",
    outputPlots = "figures/methylationGLM_T1",
    phenotypes = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore",
    covariates = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells",
    factorVars = "Sex,Ethnicity,TraumaDefinition",
    cpgPrefix = "cg",
    cpgLimit = NA,
    nCores = 32,
    plotWidth = 2000,
    plotHeight = 1000,
    plotDPI = 150,
    interactionTerm = NULL,
    libPath = NULL,
    glmLibs = "glm2",
    prsMap = NULL,
    summaryPval = NA,
    summaryResidualSD = TRUE,
    saveSignificantCpGs = FALSE,
    significantCpGDir = "preliminaryResults/cpgs/methylationGLM_T1",
    significantCpGPval = 0.05,
    saveTxtSummaries = TRUE,
    chunkSize = 10000,
    summaryTxtDir = "preliminaryResults/summary/methylationGLM_T1/glm",
    fdrThreshold = 0.05,
    padjmethod = "fdr",
    annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group",
    annotatedGLMOut = "data/methylationGLM_T1"
) {

  # Locate script inside the installed package
  script <- system.file("scripts", "methylationGLM_T1.R", package = "dnapipeR")
  if (script == "")
    stop("Script methylationGLM_T1.R not found in package.")

  # Build command-line argument list
  arg_list <- c(
    "--inputPheno",           shQuote(inputPheno),
    "--outputLogs",           shQuote(outputLogs),
    "--outputRData",          shQuote(outputRData),
    "--outputPlots",          shQuote(outputPlots),
    "--phenotypes",           shQuote(phenotypes),
    "--covariates",           shQuote(covariates),
    "--factorVars",           shQuote(factorVars),
    "--cpgPrefix",            shQuote(cpgPrefix),
    "--cpgLimit",             cpgLimit,
    "--nCores",               nCores,
    "--plotWidth",            plotWidth,
    "--plotHeight",           plotHeight,
    "--plotDPI",              plotDPI,
    "--summaryPval",          summaryPval,
    "--summaryResidualSD",    summaryResidualSD,
    "--saveSignificantCpGs",  saveSignificantCpGs,
    "--significantCpGDir",    shQuote(significantCpGDir),
    "--significantCpGPval",   significantCpGPval,
    "--saveTxtSummaries",     saveTxtSummaries,
    "--chunkSize",            chunkSize,
    "--summaryTxtDir",        shQuote(summaryTxtDir),
    "--fdrThreshold",         fdrThreshold,
    "--padjmethod",           shQuote(padjmethod),
    "--annotationPackage",    shQuote(annotationPackage),
    "--annotationCols",       shQuote(annotationCols),
    "--annotatedGLMOut",      shQuote(annotatedGLMOut)
  )

  # ---- CONDITIONAL FLAGS (skip if NULL) ----
  if (!is.null(interactionTerm))
    arg_list <- c(arg_list, "--interactionTerm", shQuote(interactionTerm))

  if (!is.null(libPath))
    arg_list <- c(arg_list, "--libPath", shQuote(libPath))

  if (!is.null(prsMap))
    arg_list <- c(arg_list, "--prsMap", shQuote(prsMap))

  # Build complete system command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))

  message("Running methylationGLM_T1:")
  message(cmd)

  # Run depending on OS
  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}

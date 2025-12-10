#' Run methylationGLMM_T1T2.R as an external script from inst/script/.R
#'
#' @export
methylationGLMM_T1T2 <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData",
    outputLogs = "logs/",
    outputRData = "rData/methylationGLMM_T1T2/models",
    outputPlots = "figures/methylationGLMM_T1T2",
    personVar = "person",
    timeVar = "Timepoint",
    phenotypes = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore",
    covariates = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells",
    factorVars = "Sex,Ethnicity,TraumaDefinition,Timepoint",
    lmeLibs = "lme4,lmerTest",
    prsMap = NULL,
    libPath = NULL,
    cpgPrefix = "cg",
    cpgLimit = NA,
    nCores = 32,
    summaryPval = NA,
    plotWidth = 2000,
    plotHeight = 1000,
    plotDPI = 150,
    interactionTerm = NULL,
    saveSignificantInteractions = TRUE,
    significantInteractionDir = "preliminaryResults/cpgs/methylationGLMM_T1T2",
    significantInteractionPval = 0.05,
    saveTxtSummaries = TRUE,
    chunkSize = 10000,
    summaryTxtDir = "preliminaryResults/summary/methylationGLMM_T1T2/lmer",
    fdrThreshold = 0.05,
    padjmethod = "fdr",
    annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group",
    annotatedLMEOut = "data/methylationGLMM_T1T2"
) {
  
  # Locate script inside installed package
  script <- system.file("scripts", "methylationGLMM_T1T2.R", package = "dnapipeR")
  if (script == "")
    stop("Script methylationGLMM_T1T2.R not found in package.")
  
  # BUILD ARGUMENT LIST
  arg_list <- c(
    "--inputPheno",               shQuote(inputPheno),
    "--outputLogs",               shQuote(outputLogs),
    "--outputRData",              shQuote(outputRData),
    "--outputPlots",              shQuote(outputPlots),
    "--personVar",                shQuote(personVar),
    "--timeVar",                  shQuote(timeVar),
    "--phenotypes",               shQuote(phenotypes),
    "--covariates",               shQuote(covariates),
    "--factorVars",               shQuote(factorVars),
    "--lmeLibs",                  shQuote(lmeLibs),
    "--cpgPrefix",                shQuote(cpgPrefix),
    "--cpgLimit",                 cpgLimit,
    "--nCores",                   nCores,
    "--summaryPval",              summaryPval,
    "--plotWidth",                plotWidth,
    "--plotHeight",               plotHeight,
    "--plotDPI",                  plotDPI,
    "--saveSignificantInteractions", saveSignificantInteractions,
    "--significantInteractionDir",   shQuote(significantInteractionDir),
    "--significantInteractionPval",  significantInteractionPval,
    "--saveTxtSummaries",            saveTxtSummaries,
    "--chunkSize",               chunkSize,
    "--summaryTxtDir",           shQuote(summaryTxtDir),
    "--fdrThreshold",            fdrThreshold,
    "--padjmethod",              shQuote(padjmethod),
    "--annotationPackage",       shQuote(annotationPackage),
    "--annotationCols",          shQuote(annotationCols),
    "--annotatedLMEOut",         shQuote(annotatedLMEOut)
  )
  
  # CONDITIONAL FLAGS
  if (!is.null(interactionTerm))
    arg_list <- c(arg_list, "--interactionTerm", shQuote(interactionTerm))
  
  if (!is.null(libPath))
    arg_list <- c(arg_list, "--libPath", shQuote(libPath))
  
  if (!is.null(prsMap))
    arg_list <- c(arg_list, "--prsMap", shQuote(prsMap))
  
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))
  
  message("Running methylationGLMM_T1T2:")
  message(cmd)
  
  # Run depending on OS
  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}

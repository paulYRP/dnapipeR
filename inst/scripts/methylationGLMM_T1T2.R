#!/usr/bin/env Rscript
# ==============================================================================
# DNAm LME Analysis Script (Longitudinal: Timepoint 1 vs 2)
# Script Name: methylationGLMM_T1T2.R
# Description: Performs CpG-by-phenotype Linear Mixed-Effects (LME) models
#              across selected phenotypes and covariates using DNA methylation
#              beta values at T1 and T2, modeling within-person changes over time.
# ==============================================================================

# Usage Example (Full version):

# Rscript methylationGLMM_T1T2.R \
#   --inputPheno rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData \
#   --outputLogs logs/methylationGLMM_T1T2 \
#   --personVar person \
#   --timeVar Timepoint \
#   --phenotypes DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore \
#   --covariates Sex,Age,Ethnicity,TraumaDefinition,Leukocytes.EWAS,Epithelial.cells.EWAS \
#   --factorVars Sex,Ethnicity,TraumaDefinition \
#   --lmeLibs lmerTest \
#   --libPath ~/R/x86_64-pc-linux-gnu-library/4.4 \
#   --cpgPrefix cg \
#   --cpgLimit 1000 \
#   --nCores 32 \
#   --interactionTerm Timepoint \
#   --saveSignificantInteractions \
#   --significantInteractionPval 0.05 \
#   --significantInteractionDir preliminaryResults/cpgs/methylationGLMM_T1T2 \
#   --saveTxtSummaries \
#   --summaryTxtDir preliminaryResults/summary/methylationGLMM_T1T2/lmer \
#   --fdrThreshold 0.05 \
#   --padjmethod BH \
#   --annotationPackage IlluminaHumanMethylationEPICv2anno.20a1.hg38 \
#   --annotationCols Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group \
#   --annotatedLMEOut data/methylationGLMM_T1T2/annotatedLME.csv
# ==============================================================================

# Usage Example (Default values with key parameters):

# d \
#   --inputPheno rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData \
#   --outputRData rData/methylationGLMM_T1T2/models
# ==============================================================================

# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(optparse)
        library(dplyr)
        library(lme4)
        library(lmerTest)
        library(parallel)
        library(ggplot2)
        library(minfi)
        library(IlluminaHumanMethylationEPICv2manifest)
        library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
        library(ggrepel)
})
# ==============================================================================

# DNAm LME Analysis Script — Input Arguments

#   --inputPheno                  [FILE]   RData file containing longitudinal pheno + beta data (T1 and T2 combined)
#   --personVar                   [STR]    Column identifying repeated individuals across timepoints
#   --timeVar                     [STR]    Column indicating timepoints (e.g., "T1", "T2")
#   --phenotypes            [STR]    Comma-separated list of phenotype variables to model (e.g., "DASS_Depression,DASS_Anxiety")
#   --covariates            [STR]    Comma-separated list of covariate names (e.g., "Sex,Age,...")
#   --factorVars            [STR]    Covariates to convert to factor, comma-separated (e.g., "Sex,Ethnicity")
#   --lmeLibs                     [STR]    Comma-separated list of libraries for LME modeling (e.g., "lmerTest")
#   --libPath                     [STR]    R library path for parallel cluster evaluation
#   --cpgPrefix                   [STR]    Prefix pattern to match CpG probe IDs (default: "cg")
#   --cpgLimit                    [INT]    Maximum number of CpGs to include (default: 1000); NA for all
#   --nCores                      [INT]    Number of cores to parallelize LME fitting
#   --interactionTerm             [STR]    Variable name for interaction with phenotype (e.g., Timepoint)
#   --saveSignificantInteractions [LOGIC]  Whether to save CpGs with significant interaction terms
#   --significantInteractionDir   [DIR]    Directory to write significant CpG interaction results
#   --significantInteractionPval  [NUM]    P-value cutoff to define significant CpG interactions
#   --saveTxtSummaries            [LOGIC]  Whether to save LME summary tables as tab-delimited TXT
#   --summaryTxtDir               [DIR]    Directory to write LME summary text outputs
#   --fdrThreshold                [NUM]    FDR threshold for significance (default: 0.05)
#   --padjmethod                  [STR]    Method for multiple testing correction (default: "BH")
#   --annotationPackage           [STR]    Annotation object name (e.g., IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#   --annotationCols              [STR]    Comma-separated columns from annotation object to retain (e.g., "Name,chr,pos,...")
#   --outputLogs                  [DIR]    Directory to save log file
#   --annotatedLMEOut             [DIR]    Directory where annotated LME summary CSV will be saved
# ==============================================================================

# ----------- Command Line Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--inputPheno", default = "rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData", help = "Input RData file with longitudinal pheno + beta values [default: %default]"),
        make_option("--outputLogs", default = "logs/", help = "Directory to save log file [default: %default]"),
        make_option("--outputRData", default = "rData/methylationGLMM_T1T2/models", help = "Directory to save LME model outputs [default: %default]"),
        make_option("--outputPlots", default = "figures/methylationGLMM_T1T2", help = "Directory to save diagnostic plots [default: %default]"),
        make_option("--personVar", default = "person", help = "Column name identifying unique individuals [default: %default]"),
        make_option("--timeVar", default = "Timepoint", help = "Column name indicating timepoints (e.g., T1, T2) [default: %default]"),
        make_option("--phenotypes", default = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore", help = "Comma-separated phenotype scores [default: %default]"),
        make_option("--covariates", default = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells", help = "Comma-separated covariates [default: %default]"),
        make_option("--factorVars", default = "Sex,Ethnicity,TraumaDefinition,Timepoint", help = "Variables to convert to factor [default: %default]"),
        make_option("--lmeLibs", default = "lme4,lmerTest", help = "Comma-separated list of libraries for LME models [default: %default]"),
        make_option("--prsMap", default = NULL, help = "Optional: comma-separated mapping of phenotype to PRS covariates (e.g., phenotype:PRS)"),
        make_option("--libPath", default = NULL, help = "Library path for parallel nodes (NULL disables setting). Aqua: ~/R/x86_64-pc-linux-gnu-library/4.4"),
        make_option("--cpgPrefix", default = "cg", help = "Regex pattern to match CpG columns [default: %default]"),
        make_option("--cpgLimit", default = NA, type = "integer", help = "Limit number of CpGs to test [default: all]"),
        make_option("--nCores", type = "integer", default = 32, help = "Number of cores for parallel processing [default: %default]"),
        make_option("--summaryPval", type = "double", default = NA, help = "Optional p-value filter threshold for summary extraction [default: none]"),
        make_option("--plotWidth", default = 2000, type = "integer", help = "Plot width in pixels [default: %default]"),
        make_option("--plotHeight", default = 1000, type = "integer", help = "Plot height in pixels [default: %default]"),
        make_option("--plotDPI", default = 150, type = "integer", help = "Plot DPI (resolution) [default: %default]"),
        make_option("--interactionTerm", default = NULL, help = "Variable to interact with phenotype in LME model [default: None]"),
        make_option("--saveSignificantInteractions", action = "store_true", default = TRUE, help = "Enable saving significant interaction CpG results [default: %default]"),
        make_option("--significantInteractionDir", default = "preliminaryResults/cpgs/methylationGLMM_T1T2", help = "Directory to save significant interaction CpGs [default: %default]"),
        make_option("--significantInteractionPval", type = "double", default = 0.05, help = "P-value threshold for significant interactions [default: %default]"),
        make_option("--saveTxtSummaries", action = "store_true", default = TRUE, help = "Whether to save GLM summaries as TXT files [default: %default]"),
        make_option("--chunkSize", type="integer", default=10000, help="Number of CpGs to process per worker [default %default]"),
        make_option("--summaryTxtDir", default = "preliminaryResults/summary/methylationGLMM_T1T2/lmer", help = "Output directory to save summary text files [default: %default]", metavar = "DIR"),
        make_option("--fdrThreshold", type = "double", default = 0.05, help = "Padj-value threshold to determine significance [default: %default]"),
        make_option("--padjmethod", default = "fdr", help = "Method for multiple testing correction [default: %default]"),
        make_option("--annotationPackage", default = "IlluminaHumanMethylationEPICv2anno.20a1.hg38", help = "Annotation object name [default: %default]"),
        make_option("--annotationCols", default = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group", help = "Comma-separated annotation columns to retain [default: %default]"),
        make_option("--annotatedLMEOut", default = "data/methylationGLMM_T1T2", help = "Path to save final annotated LME results [default: %default]")

)))

# Set default values for optional arguments
if (is.character(opt$cpgLimit) && tolower(opt$cpgLimit) == "na") {
        opt$cpgLimit <- NA
}

if (is.character(opt$summaryPval) && tolower(opt$summaryPval) == "na") {
        opt$summaryPval <- NA
}

if (!is.null(opt$prsMap)) {
        opt$prsMapList <- setNames(
                sapply(strsplit(unlist(strsplit(opt$prsMap, ",")), ":"), `[`, 2),
                sapply(strsplit(unlist(strsplit(opt$prsMap, ",")), ":"), `[`, 1)
        )
} else {
        opt$prsMapList <- list()
}

opt$lmeLibList   <- strsplit(opt$lmeLibs, ",")[[1]]
opt$covariateList <- strsplit(opt$covariates, ",")[[1]]
opt$phenotypeList <- strsplit(opt$phenotypes, ",")[[1]]
opt$factorVarsList <- strsplit(opt$factorVars, ",")[[1]]

dir.create(opt$summaryTxtDir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$significantInteractionDir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$annotatedLMEOut, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputPlots, recursive = TRUE, showWarnings = FALSE)
# ==============================================================================

# ----------- Setup Logging -----------
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)
logFilePath <- file.path(opt$outputLogs, "log_methylationGLMM_T1T2.txt")
logCon <- file(logFilePath, open = "wt")
sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting DNAm LME Analysis (Timepoint 1 vs 3) ====\n")
cat("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input phenotype + beta file: ", opt$inputPheno, "\n")
cat("Person ID variable: ", opt$personVar, "\n")
cat("Timepoint variable: ", opt$timeVar, "\n")
cat("Phenotypes: ", opt$phenotypes, "\n")
cat("Covariates: ", opt$covariates, "\n")
cat("PRS mapping: ", if (is.null(opt$prsMap)) "None" else paste(unlist(opt$prsMap), collapse = ", "), "\n")
cat("Factor variables: ", opt$factorVars, "\n")
cat("LME libraries: ", opt$lmeLibs, "\n")
cat("Library path: ", ifelse(is.null(opt$libPath), "Default (R system)", opt$libPath), "\n")
cat("CpG prefix pattern: ", opt$cpgPrefix, "\n")
cat("CpG limit: ", ifelse(is.na(opt$cpgLimit), "All CpGs", opt$cpgLimit), "\n")
cat("Number of cores: ", opt$nCores, "\n")
cat("Interaction term: ", opt$interactionTerm, "\n")
cat("Save significant interactions: ", opt$saveSignificantInteractions, "\n")
cat("Significant interaction p-value cutoff: ", opt$significantInteractionPval, "\n")
cat("Directory to save significant interactions: ", opt$significantInteractionDir, "\n")
cat("Save summary TXT files: ", opt$saveTxtSummaries, "\n")
cat("Summary p-value threshold: ", opt$summaryPval, "\n")
cat("Padjmethod: ", opt$padjmethod, "\n")
cat("FDR threshold for significance: ", opt$fdrThreshold, "\n")
cat("Multiple testing correction method: ", opt$padjmethod, "\n")
cat("Save significant interaction CpGs: ", opt$saveSignificantInteractions, "\n")
cat("Significant interaction p-value cutoff: ", opt$significantInteractionPval, "\n")
cat("Directory to save significant interaction CpGs: ", opt$significantInteractionDir, "\n")
cat("Summary TXT output directory: ", opt$summaryTxtDir, "\n")
cat("Annotation package: ", opt$annotationPackage, "\n")
cat("Annotation columns: ", opt$annotationCols, "\n")
cat("FDR threshold for significance: ", opt$fdrThreshold, "\n")
cat("Output logs directory: ", opt$outputLogs, "\n")
cat("Directory to save annotated LME results: ", opt$annotatedLMEOut, "\n")
cat("=======================================================================\n\n")

# ----------- Load Data -----------
cat("Loading input phenotype + beta file...\n")
inpt <- load(opt$inputPheno)
phenoBT1T2 <- get(inpt)

# ----------- Check and Create person Column if Missing -----------
if (!(opt$personVar %in% colnames(phenoBT1T2))) {
  cat(paste0("Column '",
             opt$personVar, "' not found. Creating it from 'SID'...\n"))

  phenoBT1T2[[opt$personVar]] <- as.numeric(factor(gsub("[AB]$", "",
                                                        phenoBT1T2$SID)))

  cat("Example mapping of SID to person ID:\n")
  print(head(phenoBT1T2[order(phenoBT1T2[[opt$personVar]],
                              phenoBT1T2$SID),
                        c("SID", opt$personVar)], 20
  ))

  cat("Count of records per person ID:\n")
  print(table(phenoBT1T2[[opt$personVar]]))
} else {
  cat(paste0("Column '",
             opt$personVar, "' already exists. Skipping creation.\n"))
}

cat("Checking structure of merged longitudinal dataset...\n")
print(table(table(phenoBT1T2[[opt$personVar]])))
print(table(phenoBT1T2[[opt$timeVar]]))
print(unique(phenoBT1T2[[opt$timeVar]]))
print(dim(phenoBT1T2))
cat("=======================================================================\n")

# ----------- Summary Stats by Timepoint -----------
cat("Summary statistics for phenotype scores by timepoint:\n")
phenoBT1T2 %>%
        select(all_of(opt$timeVar), all_of(opt$phenotypeList)) %>%
        group_by(.data[[opt$timeVar]]) %>%
        summarise(across(everything(), list(
                mean = ~mean(., na.rm = TRUE),
                sd   = ~sd(., na.rm = TRUE),
                n    = ~sum(!is.na(.))
        ))) %>%
        print(width = Inf)
cat("=======================================================================\n")

# ----------- Run LME Function -----------
lme <- function(
                phenoScore,
                merge,
                personVar = opt$personVar,
                timeVar = opt$timeVar,
                covariates = opt$covariateList,
                factorVars = opt$factorVarsList,
                interactionTerm = opt$interactionTerm,
                cpgPrefix = opt$cpgPrefix,
                cpgLimit = opt$cpgLimit,
                nCore = opt$nCores,
                libPath = opt$libPath,
                lmeLibs = opt$lmeLibList
) {
        if (!phenoScore %in% colnames(merge)) {
                stop(paste("Phenotype", phenoScore, "not found"))
        }

        covariateNames <- c(timeVar, phenoScore, covariates)
        cpgCol <- grep(paste0("^", cpgPrefix), colnames(merge), value = TRUE)
        if (!is.na(cpgLimit)) {
                cpgCol <- head(cpgCol, as.numeric(cpgLimit))
        }

        cl <- makeCluster(nCore)
        clusterExport(cl, varlist = c("merge", "phenoScore", "personVar",
                                      "timeVar", "covariateNames", "factorVars",
                                      "interactionTerm", "libPath", "lmeLibs"),
                      envir = environment())

        clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                        .libPaths(libPath)
                }
                sapply(lmeLibs, function(pkg) {
                        if (!require(pkg, character.only = TRUE)) {
                                stop(paste("Failed to load package:", pkg))
                        }
                })
        })

        fit <- function(cpg) {
                tryCatch({
                        modelData <- merge[, c(personVar, covariateNames, cpg)]
                        colnames(modelData)[ncol(modelData)] <- "beta"
                        for (var in factorVars) {
                                if (var %in% colnames(modelData)) {
                                        modelData[[var]] <- as.factor(modelData[[var]])
                                }
                        }

                        if (!is.null(interactionTerm) && interactionTerm != "") {
                          intFormula <- paste(phenoScore, "*", interactionTerm)
                          fixedTerms <- setdiff(covariateNames, c(phenoScore,
                                                                  interactionTerm))
                          fixedPart <- paste(c(intFormula, fixedTerms), collapse = " + ")
                        } else {
                          fixedTerms <- setdiff(covariateNames, phenoScore)
                          fixedPart <- paste(c(phenoScore, fixedTerms), collapse = " + ")
                        }

                        form <- as.formula(paste("beta ~",
                                                 fixedPart,
                                                 "+ (1|",
                                                 personVar, ")"))

                        model <- lmer(form,
                                      data = modelData,
                                      na.action = na.exclude, REML = TRUE)

                        list(
                                coef = summary(model)$coefficients,
                                residuals = residuals(model),
                                fitted = fitted(model),
                                ranef = ranef(model),
                                fixef = fixef(model)
                        )
                }, error = function(e) NULL)
        }

        fitList <- parLapply(cl, cpgCol, fit)
        names(fitList) <- cpgCol
        stopCluster(cl)
        return(fitList)
}

# ----------- Run LMEs for All Phenotypes -----------
cat("Running LMEs for all phenotypes...\n")
for (pheno in opt$phenotypeList) {
        outputFile <- file.path(opt$outputRData, paste0(pheno, "LME.RData"))

        if (file.exists(outputFile)) {
                cat("Loading existing LME for:", pheno, "\n")
                load(outputFile)
                assign(paste0(pheno, "LME"), fitResult, envir = .GlobalEnv)
                next
        }

        cat("Running LME for:", pheno, "\n")

        prsVar <- if (pheno %in% names(opt$prsMapList)) opt$prsMapList[[pheno]] else NULL
        allCovariates <- if (!is.null(prsVar)) c(opt$covariateList,
                                                 prsVar) else opt$covariateList

        if (!is.null(opt$interactionTerm) && opt$interactionTerm != "") {
          fixedTerms <- setdiff(allCovariates, c(pheno, opt$interactionTerm))
          modelFormula <- paste("~",
                                paste(c(paste0(pheno, "*", opt$interactionTerm),
                                        fixedTerms),
                                      collapse = " + "),
                                "+ (1|", opt$personVar, ")")
        } else {
          fixedTerms <- setdiff(allCovariates, pheno)
          modelFormula <- paste("~",
                                paste(c(pheno, fixedTerms),
                                      collapse = " + "),
                                "+ (1|", opt$personVar, ")")
        }



        cat("Formula:", modelFormula, "\n")


        fitResult <- lme(
                phenoScore = pheno,
                merge = phenoBT1T2,
                personVar = opt$personVar,
                timeVar = opt$timeVar,
                covariates = allCovariates,
                factorVars = opt$factorVarsList,
                interactionTerm = opt$interactionTerm,
                cpgPrefix = opt$cpgPrefix,
                cpgLimit = opt$cpgLimit,
                nCore = opt$nCores,
                libPath = opt$libPath,
                lmeLibs = opt$lmeLibList
        )

        save(fitResult, file = file.path(opt$outputRData,
                                         paste0(pheno, "LME.RData")))

        assign(paste0(pheno, "LME"), fitResult, envir = .GlobalEnv)

}
cat("Finished running LMEs for all phenotypes.\n")
cat("=======================================================================\n")

# ----------- Extract LME Interaction Summary -----------
cpgsLME <- function(
                fitList,
                phenotype,
                interactionTerm = opt$interactionTerm,
                pValue = opt$summaryPval,
                nCore = opt$nCores,
                libPath = opt$libPath,
                lmeLibs = opt$lmeLibList,
                chunkSize = opt$chunkSize

) {
        cat("Extracting LME interaction effects for:", phenotype, "\n")

        if (is.null(interactionTerm) || interactionTerm == "") {
          cat("No interaction term provided — extracting main effects only.\n")
        } else {
          cat("Interaction term detected:", interactionTerm, "\n")
        }


        cpgNames <- names(fitList)
        if (is.null(chunkSize)) {
          chunkSize <- max(1000, floor(length(cpgNames) / (nCore * 4)))
        }
        cat("Total CpGs:", length(cpgNames), "| Using chunkSize:", chunkSize, "\n")

        splitIntoChunks <- function(x, size) {
          split(x, ceiling(seq_along(x) / size))
        }

        cpgNames <- names(fitList)
        cpgChunks <- splitIntoChunks(cpgNames, chunkSize)

        cl <- makeCluster(nCore)
        clusterExport(
                cl,
                varlist = c("fitList", "pValue", "interactionTerm",
                            "phenotype", "libPath", "lmeLibs"),
                envir = environment()
        )

        clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                        .libPaths(libPath)
                }
                sapply(lmeLibs, function(pkg) {
                        if (!require(pkg, character.only = TRUE)) {
                                stop(paste("Failed to load package:", pkg))
                        }
                })
        })

        results <- parLapplyLB(cl, cpgChunks, function(chunk) {
          outList <- vector("list", length(chunk))
          idx <- 1
          for (cpg in chunk) {
            modelOutput <- fitList[[cpg]]
            if (is.null(modelOutput) || is.null(modelOutput$coef)) next

            coefTable <- modelOutput$coef

            if (!is.null(interactionTerm) && interactionTerm != "") {
              pattern <- paste0("^", phenotype, ".*:", interactionTerm)
              matchedTerms <- grep(pattern, rownames(coefTable), value = TRUE)
            } else {
              pattern <- paste0("^", phenotype)
              matchedTerms <- grep(pattern, rownames(coefTable), value = TRUE)
            }

            if (length(matchedTerms) == 0) next

            tmp <- tryCatch({
              do.call(rbind, lapply(matchedTerms, function(term) {
                coefRow <- coefTable[term, ]
                data.frame(
                  CpG = cpg,
                  Interaction.Term = term,
                  Estimate = coefRow["Estimate"],
                  Std.Error = coefRow["Std. Error"],
                  t.value = coefRow["t value"],
                  P.value = coefRow["Pr(>|t|)"],
                  row.names = NULL
                )
              }))})
            if (!is.null(tmp)) {
              outList[[idx]] <- tmp
              idx <- idx + 1
              }
          }
          if (idx > 1) {
            return(data.table::rbindlist(outList[1:(idx-1)],
                                         use.names = TRUE, fill = TRUE))
          } else {
            return(NULL)
          }
        })

        stopCluster(cl)

        summary <- do.call(rbind, results)

        if (!is.na(pValue)) {
                summary <- subset(summary, `P.value` < pValue)
        }


        cat("Finished extracting interaction for:", phenotype, "\n")
        return(summary)
}

# ----------- Run and Save LME Summaries -----------
cat("Running LME summary extraction for all phenotypes...\n")

for (pheno in opt$phenotypeList) {
        summaryFile <- file.path(opt$outputRData, paste0(pheno,
                                                         "SummaryLME.RData"))

        if (file.exists(summaryFile)) {
                cat("Loading existing summary for:", pheno, "\n")
                load(summaryFile)
                assign(paste0(pheno, "SummaryLME"), fitResult,
                       envir = .GlobalEnv)
                next
        }

        cat("Processing interaction effects for:", pheno, "\n")

        fitObjectName <- paste0(pheno, "LME")
        if (!exists(fitObjectName)) {
                cat("Object not found:", fitObjectName, "- skipping\n")
                next
        }

        fitObject <- get(fitObjectName)

        fitResult <- cpgsLME(
                fitList = fitObject,
                phenotype = pheno,
                nCore = opt$nCores,
                pValue = opt$summaryPval,
                libPath = opt$libPath,
                lmeLibs = opt$lmeLibList,
                interactionTerm = opt$interactionTerm
        )

        save(fitResult, file = summaryFile)

        assign(paste0(pheno, "SummaryLME"), fitResult, envir = .GlobalEnv)

        cat("Saved:", paste0(pheno, "SummaryLME.RData"), "\n")
}

cat("Finished extracting all LME summaries.\n")
cat("=======================================================================\n")

# ----------- Save Significant LME Interaction Results -----------
saveSignificantInteractions <- function(
                resultList,
                resultName = deparse(substitute(resultList)),
                baseDir = opt$significantInteractionDir,
                pvalThreshold = opt$significantInteractionPval,
                interactionTerm = opt$interactionTerm
) {
        resultDir <- file.path(baseDir, resultName)
        if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

        if (is.null(interactionTerm) || interactionTerm == "") {
          cat("No interaction term detected — extracting main effects for",
              resultName, "\n")
          pattern <- paste0("^", resultName)
        } else {
          cat("Interaction term detected:", interactionTerm,
              "- extracting interaction effects for", resultName, "\n")
          pattern <- paste0("^", resultName, ".*:", interactionTerm)
        }

        for (i in seq_along(resultList)) {
          coefTable <- resultList[[i]]$coef
          cpgName <- names(resultList)[i]
          matchedRows <- grep(pattern, rownames(coefTable), value = FALSE)

          if (length(matchedRows) > 0) {
            termPvals <- coefTable[matchedRows, "Pr(>|t|)"]
            if (any(termPvals  < pvalThreshold, na.rm = TRUE)) {
              cpgDir <- file.path(resultDir, cpgName)
              if (!dir.exists(cpgDir)) dir.create(cpgDir)

              outputFile <- file.path(cpgDir, paste0(cpgName, ".txt"))
              write.table(coefTable, file = outputFile,
                          sep = "\t", quote = FALSE)
            }
          }
        }
}

cat("Saving significant LME interaction results...\n")

cat("Saving significant interaction CpGs to:", opt$significantInteractionDir, "\n")
for (pheno in opt$phenotypeList) {
        modelObjName <- paste0(pheno, "LME")
        if (exists(modelObjName)) {
                resultObj <- get(modelObjName)
                saveSignificantInteractions(resultList = resultObj, resultName = pheno)
        }
}
cat("Finished saving significant CpG interaction outputs.\n")
cat("=======================================================================\n")

# ----------- Save Summary Tables as TXT -----------
saveSummaryToTxt <- function(
                summaryDF,
                outputFile
) {
        if ("P.value" %in% colnames(summaryDF)) {
                summaryDF <- summaryDF[order(summaryDF$P.value), ]
        } else if ("Pr(>|t|)" %in% colnames(summaryDF)) {
                summaryDF <- summaryDF[order(summaryDF[["Pr(>|t|)"]]), ]
        }

        dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)

        write.table(
                summaryDF,
                file = outputFile,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
        )

        cat("Saved summary table to:", outputFile, "\n")
}

# ----------- Run and Save TXT Summaries -----------

if (opt$saveTxtSummaries) {
        cat("Saving summaries to TXT files...\n")

        phenoList <- strsplit(opt$phenotypes, ",")[[1]]

        for (pheno in phenoList) {

                outputFile <- file.path(
                        opt$summaryTxtDir,
                        paste0(pheno, "SummaryLME.txt")
                )

                summaryObjName <- paste0(pheno, "SummaryLME")

                if (exists(summaryObjName)) {
                        summaryObj <- get(summaryObjName)

                        saveSummaryToTxt(
                                summaryDF = summaryObj,
                                outputFile = outputFile
                        )
                }
        }

        cat("Finished writing all GLM summaries to TXT files.\n")
}
cat("=======================================================================\n")

# ----------- Diagnostic Plots for LME Results -----------
diagnosticPlotsLME <- function(summary,
                               phenoBeta,
                               variable,
                               fdrThreshold = opt$fdrThreshold,
                               padjmethod = opt$padjmethod,
                               cpgPrefix = opt$cpgPrefix,
                               outputDir = opt$outputPlots,
                               plotWidth = opt$plotWidth,
                               plotHeight = opt$plotHeight,
                               plotDPI = opt$plotDPI
                               ) {

        ## Multiple Testing Correction
        summary$FDR <- p.adjust(summary$P.value, method = padjmethod)

        ## Genomic Inflation Factor (λ)
        chisq <- qchisq(1 - summary$P.value, df = 1)
        lambda <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1), 3)
        message("Genomic inflation factor (λ): ", lambda)

        ## Q-Q Plot
        pvals <- summary$P.value
        pvals <- pvals[!is.na(pvals)]

        tiff(filename = file.path(outputDir,
                                  paste0("qqplot_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")

        qqplot(-log10(ppoints(length(pvals))), -log10(sort(pvals)),
               main = paste("Q-Q Plot of p-values for", variable,
                            "\nGenomic Inflation Factor λ =", lambda),
               xlab = "Expected -log10(p)",
               ylab = "Observed -log10(p)",
               pch = 16, col = "black")
        abline(0, 1, col = "red")
        dev.off()

        ## Mean Beta Values
        cpgCols <- grep(paste0("^", cpgPrefix),
                        colnames(phenoBeta), value = TRUE)
        betas <- phenoBeta[, cpgCols]
        meanBeta <- colMeans(betas, na.rm = TRUE)
        meanBetaLog <- log2(meanBeta)

        # Merge log2meanBeta with summary
        summary$log2meanBeta <- meanBetaLog[summary$CpG]
        summary <- summary[!is.na(summary$log2meanBeta), ]

        ## Plot: Residual Proxy (SD estimate) vs Mean Methylation
        tiff(filename = file.path(outputDir,
                                  paste0("residualSD_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")

        plot(summary$log2meanBeta,
             summary$Std.Error,
             pch = 20,
             col = rgb(0, 0, 0, 0.4),
             xlab = "log2(Average Beta)",
             ylab = "Standard Error",
             main = "SD vs Average Beta (proxy from Std. Error)")
        lines(lowess(summary$log2meanBeta, summary$Std.Error),
              col = "red", lwd = 2)
        dev.off()

        ## Plot: Significance vs Variability (colored by FDR)
        p <- ggplot(summary, aes(x = -log10(P.value),
                            y = Std.Error, color = FDR < fdrThreshold)) +
                geom_point(alpha = 0.6) +
                geom_text_repel(data = subset(summary,
                                              FDR < fdrThreshold),
                                aes(label = CpG)) +
                scale_color_manual(values = c("FALSE" = "grey50",
                                              "TRUE" = "firebrick")) +
                labs(title = paste("Standard Error vs Significance for", variable),
                     x = "-log10(p-value)",
                     y = "Standard Error",
                     color = paste("FDR <", fdrThreshold)) +
                theme_minimal()

        tiff(
                filename = file.path(outputDir,
                                     paste0("residualSignificance_",
                                            variable, ".tiff")),
                width = plotWidth,
                height = plotHeight,
                res = plotDPI, type = "cairo"
        )
        print(p)
        dev.off()
}
# ==============================================================================

# ----------- Generate Diagnostic Plots for All Phenotypes -----------
cat("Generating diagnostic plots for all phenotypes...\n")

# Parse phenotype list
phenotypeNames <- strsplit(opt$phenotypes, ",")[[1]]

for (pheno in phenotypeNames) {

        summaryName <- paste0(pheno, "SummaryLME")

        if (exists(summaryName)) {
                cat("Generating diagnostic plots for:", pheno, "\n")

                diagnosticPlotsLME(
                        summary = get(summaryName),
                        phenoBeta = phenoBT1T2,
                        variable = pheno,
                        cpgPrefix = opt$cpgPrefix,
                        fdrThreshold = opt$fdrThreshold,
                        padjmethod = opt$padjmethod,
                        outputDir = opt$outputPlots,
                        plotWidth = opt$plotWidth,
                        plotHeight = opt$plotHeight,
                        plotDPI = opt$plotDPI
                )
        } else {
                cat("Warning: Summary object not found for phenotype:", pheno, "\n")
        }
}

cat("Finished generating all diagnostic plots.\n")
cat("Plots saved to:", opt$outputPlots, "\n")
cat("=======================================================================\n")

# ----------- Load Annotation Object -----------
cat("Loading annotation object:", opt$annotationPackage, "\n")
annotationObject <- getAnnotation(get(opt$annotationPackage))

cat("Annotation loaded with", nrow(annotationObject), "probes\n")
cat("Available columns:\n")
print(colnames(annotationObject))
cat("=======================================================================\n")

# ----------- Function: Annotate LME Summaries -----------
annotateLME <- function(
                summaryList,
                annotationObject,
                annotationCols = strsplit(opt$annotationCols, ",")[[1]]
) {
        modelNames <- names(summaryList)

        cat("Merging LME summaries...\n")

        cleanedSummaries <- list()

        for (modelName in modelNames) {
                df <- summaryList[[modelName]]

                if (!all(c("CpG", "Interaction.Term", "P.value") %in% colnames(df))) {
                        warning(paste("Skipping", modelName, ": required columns not found"))
                        next
                }

                dfSplit <- split(df, df$Interaction.Term)

                modelTables <- lapply(names(dfSplit), function(term) {
                        subDf <- dfSplit[[term]]

                        interactionSuffix <- gsub(paste0("^", modelName, "\\."), "", term)
                        pCol <- paste0(modelName, "_", interactionSuffix, "_P.Value")

                        subDf[[pCol]] <- subDf$P.value
			          subDf <- as.data.frame(subDf)
                        subDf <- subDf[, c("CpG", pCol)]
                        return(subDf)
                })

                mergedModel <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE),
                                      modelTables)
                cleanedSummaries[[modelName]] <- mergedModel
        }

        mergedSummary <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE),
                                cleanedSummaries)

        annDF <- as.data.frame(annotationObject)

	annDF$CpG <- rownames(annDF)

        if (!is.null(annotationCols)) {
                annDF <- annDF[, unique(c("CpG", annotationCols)), drop = FALSE]
        }

        annotatedResults <- merge(mergedSummary, annDF, by = "CpG", all.x = TRUE)
        colnames(annotatedResults)[1] <- "IlmnID"

        cat("Annotation completed. Annotated CpGs:", nrow(annotatedResults), "\n")
        return(annotatedResults)
}

# ----------- Execute Annotation and Save CSV -----------
cat("Running annotation of LME summary results...\n")

# Build summaryList dynamically
phenotypeNames <- strsplit(opt$phenotypes, ",")[[1]]
summaryList <- setNames(
        lapply(phenotypeNames, function(pheno) {
                get(paste0(pheno, "SummaryLME"))
        }),
        phenotypeNames
)

cat("\n Phenotypes to annotate:\n")
print(names(summaryList))

for (nm in names(summaryList)) {
  cat("First few rows of", nm, "summary:\n")
  print(head(summaryList[[nm]]))
  cat("Column names in", nm, "summary:\n")
  print(colnames(summaryList[[nm]]))
}

# Annotate
annotatedLME <- annotateLME(
        summaryList = summaryList,
        annotationObject = annotationObject,
        annotationCols = strsplit(opt$annotationCols, ",")[[1]]
)

# Save as CSV inside the output directory
annotatedLMEPath <- file.path(opt$annotatedLMEOut, "annotatedLME.csv")

cat("Saving annotated LME results to:", annotatedLMEPath, "\n")

write.csv(
        annotatedLME,
        file = annotatedLMEPath,
        row.names = FALSE
)

cat("Finished writing annotated LME results.\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)

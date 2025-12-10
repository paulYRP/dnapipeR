# !/usr/bin/env Rscript
# ==============================================================================
# DNAm GLM Analysis Script (Timepoint 1)
# Script Name: methylationGLM_T1.R
# Description: Performs CpG-by-phenotype GLMs across selected phenotypes
#              and covariates using DNA methylation beta values at T1.
# ==============================================================================

# Usage Example (Full version):

# Rscript methylationGLM_T1.R \
#   --inputPheno rData/preprocessingPheno/mergeData/phenoBetaT1.RData \
#   --outputLogs logs/methylationGLM_T1 \
#   --outputRData rData/methylationGLM_T1/models \
#   --outputPlots figures/methylationGLM_T1 \
#   --phenotypes PCL_SUM,PCL5_B,PCL5_C,PCL5_D,PCL5_E,PTGIX_SUM,DASS_D,DASS_S,DASS_A,SSS8_SUM \
#   --covariates Sex,Age,Ethnicity,TraumaDefinition,Leukocytes.EWAS,Epithelial.cells.EWAS,PTSD_PRS \
#   --factorVars Sex,Ethnicity,TraumaDefinition \
#   --cpgPrefix cg \
#   --cpgLimit 100000 \
#   --nCores 32 \
#   --plotWidth 2000 \
#   --plotHeight 1000 \
#   --plotDPI 150 \
#   --libPath ~/R/x86_64-pc-linux-gnu-library/4.4 \
#   --glmLibs glm2 \
#   --prsMap DASS_D:MDD_PRS,DASS_S:PTSD_PRS,DASS_A:GAD_PRS,SSS8_SUM:MDD_PRS \
#   --summaryPval 0.01 \
#   --summaryResidualSD \
#   --saveSignificantCpGs \
#   --significantCpGPval 0.05 \
#   --saveTxtSummaries \
#   --summaryTxtDir preliminaryResults/summary/methylationGLM_T1/glm \
#   --fdrThreshold 0.05 \
#   --padjmethod BH \
#   --annotationPackage IlluminaHumanMethylationEPICv2anno.20a1.hg38 \
#   --annotationCols Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group \
#   --annotatedGLMOut data/methylationGLM_T1/annotatedGLM.csv
# ==============================================================================

# Usage Example (Default values with key parameters):

# Rscript methylationGLM_T1.R \
#   --inputPheno rData/preprocessingPheno/mergeData/phenoBetaT1.RData \
#   --outputRData rData/methylationGLM_T1/models
# ==============================================================================

# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(optparse)
        library(parallel)
        library(glm2)
        library(ggplot2)
        library(minfi)
        library(IlluminaHumanMethylationEPICv2manifest)
        library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
        library(ggrepel)
})
# ==============================================================================

# DNAm GLM Analysis Script — Input Arguments

#   --inputPheno            [FILE]   RData file containing phenotype + beta data (e.g., phenoBetaT1.RData)
#   --outputLogs            [DIR]    Directory to save log files and execution metadata
#   --outputRData           [DIR]    Directory to save fitted GLM objects and summaries
#   --outputPlots           [DIR]    Folder to store diagnostic and output plots (.tiff)
#   --phenotypes            [STR]    Comma-separated list of phenotype variables to model (e.g., "DASS_Depression,DASS_Anxiety")
#   --covariates            [STR]    Comma-separated list of covariate names (e.g., "Sex,Age,...")
#   --factorVars            [STR]    Covariates to convert to factor, comma-separated (e.g., "Sex,Ethnicity")
#   --cpgPrefix             [STR]    Prefix pattern to match CpG probe IDs (e.g., "cg|ch")
#   --cpgLimit              [INT]    Maximum number of CpGs to include in GLMs (e.g., 100000); NA for all
#   --nCores                [INT]    Number of parallel cores to use for fitting GLMs
#   --plotWidth             [INT]    Width of plots in pixels (default: 2000)
#   --plotHeight            [INT]    Height of plots in pixels (default: 1000)
#   --plotDPI               [INT]    Plot resolution in dots per inch (DPI)
#   --libPath               [STR]    Path to R libraries used in parallel workers
#   --glmLibs               [STR]    Comma-separated list of libraries to load in GLM evaluation (e.g., "glm2")
#   --prsMap                [STR]    Optional: Comma-separated mapping of phenotype to PRS covariates (e.g., "phenotype:PRS")
#   --summaryPval           [NUM]    P-value threshold to retain significant GLM summary rows
#   --summaryResidualSD     [LOGIC]  Whether to include residual SD in the GLM summary output
#   --saveSignificantCpGs   [LOGIC]  Enable saving CpG coefficient tables for significant CpGs only
#   --significantCpGDir     [DIR]    Directory to write significant CpG coefficient tables
#   --significantCpGPval    [NUM]    P-value cutoff to define CpGs as significant
#   --saveTxtSummaries      [LOGIC]  Whether to save GLM summary tables as tab-delimited text files
#   --summaryTxtDir         [DIR]    Directory to write GLM summary text files
#   --fdrThreshold          [NUM]    FDR threshold for significant CpGs (default: 0.05)
#   --padjmethod            [STR]    Method for multiple testing correction (e.g., "BH", "bonferroni")
#   --annotationPackage     [STR]    Annotation object name (e.g., IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#   --annotationCols        [STR]    Comma-separated columns from annotation object to include (e.g., "Name,chr,pos,...")
#   --annotatedGLMOut       [FILE]   Output CSV file to save final annotated GLM table
# ==============================================================================

# ----------- Command Line Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--inputPheno", default = "rData/preprocessingPheno/mergeData/phenoBetaT1.RData", help = "Input RData file with pheno+beta [default: %default]"),
        make_option("--outputLogs", default = "logs/", help = "Directory to save log file [default: %default]"),
        make_option("--outputRData", default = "rData/methylationGLM_T1/models", help = "Directory to save GLM result objects [default: %default]"),
        make_option("--outputPlots", default = "figures/methylationGLM_T1", help = "Directory for output plots [default: %default]"),
        make_option("--phenotypes", default = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore", help = "Comma-separated phenotype scores [default: %default]"),
        make_option("--covariates", default = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells", help = "Comma-separated covariates [default: %default]"),
        make_option("--factorVars", default = "Sex,Ethnicity,TraumaDefinition", help = "Variables to convert to factor [default: %default]"),
        make_option("--cpgPrefix", default = "cg", help = "Regex pattern to match CpG columns cg|ch [default: %default]"),
        make_option("--cpgLimit", default = NA, type = "integer", help = "Limit number of CpGs to test [default: all]"),
        make_option("--nCores", default = 32, type = "integer", help = "Number of CPU cores to use [default: %default]"),
        make_option("--plotWidth", default = 2000, type = "integer", help = "Plot width in pixels [default: %default]"),
        make_option("--plotHeight", default = 1000, type = "integer", help = "Plot height in pixels [default: %default]"),
        make_option("--plotDPI", default = 150, type = "integer", help = "Plot DPI (resolution) [default: %default]"),
        make_option("--interactionTerm", default = NULL, help = "Optional interaction with phenotype [default: %default]"),
        make_option("--libPath", default = NULL, help = "Default library path used in clusterEvalQ. Aqua: ~/R/x86_64-pc-linux-gnu-library/4.4"),
        make_option("--glmLibs", default = "glm2", help = "Comma-separated list of libraries to load in each worker node"),
        make_option("--prsMap", default = NULL, help = "Optional: comma-separated mapping of phenotype to PRS covariates (e.g., phenotype:PRS)"),
        make_option("--summaryPval", type = "double", default = NA, help = "Optional p-value filter threshold for summary extraction [default: none]"),
        make_option("--summaryResidualSD", action = "store_true", default = TRUE, help = "Include residual standard deviation in summary [default: %default]"),
        make_option("--saveSignificantCpGs", action = "store_true", default = FALSE, help = "Enable saving significant CpG results for each GLM model [default: %default]"),
        make_option("--significantCpGDir", default = "preliminaryResults/cpgs/methylationGLM_T1", help = "Directory to store significant CpG tables [default: %default]"),
        make_option("--significantCpGPval", type = "double", default = 0.05, help = "P-value threshold to determine significance [default: %default]"),
        make_option("--saveTxtSummaries", action = "store_true", default = TRUE, help = "Whether to save GLM summaries as TXT files [default: %default]"),
        make_option("--chunkSize", type="integer", default=10000, help="Number of CpGs to process per worker [default %default]"),
        make_option("--summaryTxtDir", default = "preliminaryResults/summary/methylationGLM_T1/glm", help = "Output directory to save summary text files [default: %default]", metavar = "DIR"),
        make_option("--fdrThreshold", type = "double", default = 0.05, help = "FDR threshold for significant CpGs [default: %default]"),
        make_option("--padjmethod", default = "fdr", help = "Method for multiple testing correction [default: %default]"),
        make_option("--annotationPackage", default = "IlluminaHumanMethylationEPICv2anno.20a1.hg38", help = "Annotation object name [default: %default]"),
        make_option("--annotationCols", default = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group", help = "Comma-separated annotation columns to retain [default: %default]"),
        make_option("--annotatedGLMOut", default = "data/methylationGLM_T1", help = "Path to save final annotated GLM results [default: %default]")


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


dir.create(opt$annotatedGLMOut, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$summaryTxtDir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$significantCpGDir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputPlots, recursive = TRUE, showWarnings = FALSE)

opt$glmLibList   <- strsplit(opt$glmLibs, ",")[[1]]
opt$covariateList <- strsplit(opt$covariates, ",")[[1]]
opt$phenotypeList <- strsplit(opt$phenotypes, ",")[[1]]
opt$factorVarsList <- strsplit(opt$factorVars, ",")[[1]]
# ==============================================================================

# ----------- Setup Logging -----------
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)
logFilePath <- file.path(opt$outputLogs, "log_methylationGLM_T1.txt")
logCon <- file(logFilePath, open = "wt")
sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting DNAm GLM Analysis (Timepoint 1) ====\n")
cat("Start time: ", format(Sys.time()), "\n")
cat("Input phenotype + beta file: ", opt$inputPheno, "\n")
cat("Output RData folder: ", opt$outputRData, "\n")
cat("Output logs folder: ", opt$outputLogs, "\n")
cat("Output plots folder: ", opt$outputPlots, "\n")
cat("Phenotypes: ", opt$phenotypes, "\n")
cat("Covariates: ", opt$covariates, "\n")
cat("Factor variables: ", opt$factorVars, "\n")
cat("CpG column prefix: ", opt$cpgPrefix, "\n")
cat("CpG limit: ", ifelse(is.na(opt$cpgLimit), "All", opt$cpgLimit), "\n")
cat("Number of cores: ", opt$nCores, "\n")
cat("Interaction term: ", opt$interactionTerm, "\n")
cat("Library path: ", opt$libPath, "\n")
cat("GLM libraries: ", opt$glmLibs, "\n")
cat("PRS mapping: ", ifelse(is.null(opt$prsMap), "None", paste(opt$prsMap, collapse = ", ")), "\n")
cat("Include Residual SD in summary: ", opt$summaryResidualSD, "\n")
cat("Summary p-value filter: ", ifelse(is.na(opt$summaryPval), "None", opt$summaryPval), "\n")
cat("Save TXT summaries: ", opt$saveTxtSummaries, "\n")
cat("Save summary tables: ", opt$saveTxtSummaries, "\n")
cat("Summary output folder: ", opt$summaryTxtDir, "\n")
cat("FDR threshold: ", opt$fdrThreshold, "\n")
cat("P-value adjustment method: ", opt$padjmethod, "\n")
cat("Save significant CpGs: ", opt$saveSignificantCpGs, "\n")
cat("Significant CpG output folder: ", opt$significantCpGDir, "\n")
cat("Significance p-value cutoff: ", opt$significantCpGPval, "\n")
cat("Annotation package: ", opt$annotationPackage, "\n")
cat("Annotation columns: ", opt$annotationCols, "\n")
cat("Annotated output CSV: ", opt$annotatedGLMOut, "\n")
cat("=======================================================================\n")

# ----------- Load Data -----------
inpt <- load(opt$inputPheno)
phenoBT1 <- get(inpt)

cat("Loaded phenotype + beta data from:", opt$inputPheno, "\n")
phenotypes <- strsplit(opt$phenotypes, ",")[[1]]
covariates <- strsplit(opt$covariates, ",")[[1]]
factorVars <- strsplit(opt$factorVars, ",")[[1]]

# ----------- Analysis Info -----------
cat("phenoBT1 dimensions:", dim(phenoBT1), "\n")
cgCount <- sum(grepl("^cg", colnames(phenoBT1)))
chCount <- sum(grepl("^ch", colnames(phenoBT1)))
cat("Number of CpG columns (cg):", cgCount, "\n")
cat("Number of CpG columns (ch):", chCount, "\n")

cat("Checking structure of the dataset...\n")
if (!is.null(opt$interactionTerm) && opt$interactionTerm != "") {
  if (opt$interactionTerm %in% names(phenoBT1)) {
    print(table(phenoBT1[[opt$interactionTerm]], useNA = "ifany"))
  } else {
    cat("Warning: interactionTerm not found in phenoBT1 columns.\n")
  }
} else {
  cat("No interaction term specified; skipping table summary.\n")
}
cat("=======================================================================\n")

# ----------- Convert Factors -----------
for (var in factorVars) {
        if (var %in% colnames(phenoBT1)) {
                phenoBT1[[var]] <- as.factor(phenoBT1[[var]])
        }
}

# ----------- Missing Summary & Distribution -----------
cat("Missing summary:\n")
print(sapply(phenoBT1[, c(phenotypes, covariates)], function(x) sum(is.na(x))))
cat("=======================================================================\n")

cat("Summary statistics:\n")
print(summary(phenoBT1[, c(phenotypes, covariates)]))

for (var in phenotypes) {

        if (is.numeric(phenoBT1[[var]])) {
          p <- ggplot(phenoBT1, aes_string(x = var)) +
            geom_histogram(bins = 30, fill = "steelblue", color = "white") +
            labs(title = paste("Distribution of", var),
                 x = var, y = "Frequency") +
            theme_minimal()

        } else {
          p <- ggplot(phenoBT1, aes_string(x = var)) +
            geom_bar(fill = "steelblue") +
            labs(title = paste("Distribution of", var), x = var, y = "Count") +
            theme_minimal()
        }

        tiff(filename = file.path(opt$outputPlots, paste0("hist_", var, ".tiff")),
             width = opt$plotWidth,
             height = opt$plotHeight,
             res = opt$plotDPI, type = "cairo")

        print(p)
        dev.off()
}

catVars <- intersect(factorVars, colnames(phenoBT1))
for (var in catVars) {
        p <- ggplot(phenoBT1, aes_string(x = var)) +
                geom_bar(fill = "darkorange") +
                labs(title = paste("Distribution of",
                                   var), x = var, y = "Count") +
                theme_minimal()
        tiff(filename = file.path(opt$outputPlots, paste0("bar_", var, ".tiff")),
             width = opt$plotWidth,
             height = opt$plotHeight,
             res = opt$plotDPI, type = "cairo")
        print(p)
        dev.off()
}

numVars <- setdiff(covariates, factorVars)
numVars <- intersect(numVars, colnames(phenoBT1))

for (var in numVars) {
        p <- ggplot(phenoBT1, aes_string(x = var)) +
                geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
                labs(title = paste("Distribution of", var),
                     x = var, y = "Frequency") +
                theme_minimal()

        tiff(filename = file.path(opt$outputPlots, paste0("hist_", var, ".tiff")),
             width = opt$plotWidth,
             height = opt$plotHeight,
             res = opt$plotDPI, type = "cairo")
        print(p)
        dev.off()
}

cat("Plots saved to:", opt$outputPlots, "\n")
cat("=======================================================================\n")

# ----------- GLM Function -----------
glm <- function(
                phenoScore,
                merge,
                covariates = opt$covariateList,
                factorVars = opt$factorVarsList,
                cpgPrefix = opt$cpgPrefix,
                cpgLimit = opt$cpgLimit,
                nCore = opt$nCores,
                libPath = opt$libPath,
                glmLibList = opt$glmLibList,
                interactionTerm = NULL

) {
        if (!phenoScore %in% colnames(merge)) {
                stop(paste("Phenotype", phenoScore, "not found in the dataset"))
        }

        if (!is.null(interactionTerm) && interactionTerm != "" && interactionTerm %in% colnames(merge)) {
          interactionPart <- paste(phenoScore, "*", interactionTerm)
          fixedTerms <- setdiff(c(covariates), c(phenoScore, interactionTerm))
          fullFormula <- paste("~", paste(c(interactionPart, fixedTerms), collapse = " + "))
        } else {
          fullFormula <- paste("~", paste(c(phenoScore, covariates), collapse = " + "))
        }

        cpgCol <- grep(paste0("^", cpgPrefix), colnames(merge), value = TRUE)
        if (!is.na(cpgLimit)) {
                cpgCol <- head(cpgCol, as.numeric(cpgLimit))
        }
        cl <- makeCluster(nCore)
        clusterExport(
                cl,
                varlist = c("phenoScore", "interactionTerm",
                            "covariates",
                            "factorVars",
                            "libPath",
                            "glmLibList", "fullFormula"),
                envir = environment()
        )

        clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                        .libPaths(libPath)
                }
                sapply(glmLibList, function(pkg) {
                        if (!require(pkg, character.only = TRUE)) {
                                stop(paste("Failed to load package:", pkg))
                        }
                })
        })

        fit <- function(cpg) {
                tryCatch({
                        vars <- unique(c(phenoScore,
                                                 covariates,
                                                 cpg,
                                                 if (!is.null(interactionTerm) && interactionTerm != "") interactionTerm))

                        model <- merge[, vars, drop = FALSE]
                        colnames(model)[ncol(model)] <- "beta"

                        for (var in factorVars) {
                                if (var %in% colnames(model)) {
                                        model[[var]] <- as.factor(model[[var]])
                                }
                        }

                        fit <- glm2(
                                formula = as.formula(paste("beta", fullFormula)),
                                data = model,
                                family = gaussian(),
                                na.action = na.exclude
                        )

                        return(list(
                                coef = summary(fit)$coefficients,
                                residuals = residuals(fit),
                                fitted = fitted(fit)
                        ))
                }, error = function(e) NULL)
        }

        fitList <- parLapply(cl, cpgCol, fit)
        names(fitList) <- cpgCol
        stopCluster(cl)

        return(fitList)
}

# ----------- Run GLMs -----------
cat("Running GLMs for all phenotypes...\n")
for (pheno in opt$phenotypeList) {
        outputFile <- file.path(opt$outputRData, paste0(pheno, "GLM.RData"))

        if (file.exists(outputFile)) {
                cat("Loading existing GLM for:", pheno, "\n")
                load(outputFile)
                assign(paste0(pheno, "GLM"), fitResult, envir = .GlobalEnv)
                next
        }

        cat("Running GLM for:", pheno, "\n")

        # Include PRS if defined for this phenotype
        prsVar <- if (pheno %in% names(opt$prsMapList)) opt$prsMapList[[pheno]] else NULL
        allCovariates <- if (!is.null(prsVar)) c(opt$covariateList, prsVar) else opt$covariateList

        # Interaction
        modelFormula <- if (!is.null(opt$interactionTerm) && opt$interactionTerm != "") {
          paste("~", paste(c(paste0(pheno, "*", opt$interactionTerm), allCovariates), collapse = " + "))
        } else {
          paste("~", paste(c(pheno, allCovariates), collapse = " + "))
        }
        cat("Formula:", modelFormula, "\n")

        fitResult <- glm(
                phenoScore = pheno,
                merge = phenoBT1,
                covariates = allCovariates,
                factorVars = opt$factorVarsList,
                cpgPrefix = opt$cpgPrefix,
                cpgLimit = opt$cpgLimit,
                nCore = opt$nCores,
                libPath = opt$libPath,
                glmLibList = opt$glmLibList,
                interactionTerm = opt$interactionTerm
        )

        save(fitResult, file = outputFile)

        assign(paste0(pheno, "GLM"), fitResult, envir = .GlobalEnv)

}
cat("Finished running GLMs for all phenotypes.\n")
cat("=======================================================================\n")

# ----------- Extract GLM Summary in Parallel -----------
cpgsGLM <- function(
                fitList,
                variable,
                interactionTerm = NULL,
                includeResidualSD = opt$summaryResidualSD,
                pValue = opt$summaryPval,
                nCore = opt$nCores,
                libPath = opt$libPath,
                glmLibList = opt$glmLibList,
                chunkSize = opt$chunkSize

) {
  cat("Starting extraction for variable:", variable, "\n")
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
    varlist = c("fitList",
                "variable",
                "interactionTerm",
                "includeResidualSD",
                "pValue",
                "libPath",
                "glmLibList"),
    envir = environment()
  )

  clusterEvalQ(cl, {
    if (!is.null(libPath)) {
      .libPaths(libPath)
    }
    sapply(glmLibList, function(pkg) {
      if (!require(pkg, character.only = TRUE)) {
        stop(paste("Failed to load package:", pkg))
      }
    })
  })

  results <- parLapplyLB(cl, cpgChunks, function(chunk) {
    outList <- vector("list", length(chunk))
    idx <- 1
    for (cpg in chunk) {
      modelObj <- fitList[[cpg]]
      if (is.null(modelObj) || is.null(modelObj$coef)) next

      coefTable <- modelObj$coef
      if (!is.null(interactionTerm) && interactionTerm != "") {
        pattern <- paste0("^",
                          variable, ".*:", interactionTerm)
        matchedRows <- grep(pattern, rownames(coefTable), value = TRUE)

        if (length(matchedRows) == 0) {
          matchedRows <- grep(paste0("^", variable),
                              rownames(coefTable), value = TRUE)
        }
      } else {
        matchedRows <- grep(paste0("^", variable),
                            rownames(coefTable), value = TRUE)
      }

      if (length(matchedRows) == 0) next

      tmp <- coefTable[matchedRows, , drop = FALSE]
      tmp <- as.data.frame(tmp)
      tmp$CpG <- cpg
      tmp$Coefficient <- rownames(tmp)

      if (includeResidualSD && !is.null(modelObj$residuals)) {
        tmp$ResidualSD <- sd(modelObj$residuals, na.rm = TRUE)
      }

      if (!is.null(tmp)) {
        outList[[idx]] <- tmp
        idx <- idx + 1
      }
    }
    do.call(rbind, outList)
  })

  stopCluster(cl)

  summary <- do.call(rbind, results)

  if (is.null(summary) || nrow(summary) == 0) {
    warning("No CpG-level results extracted for variable:", variable)
    return(NULL)
  }

  cols <- c("CpG", "Coefficient", "Estimate", "Std. Error", "t value",
            "Pr(>|t|)")
  if (includeResidualSD) cols <- c(cols, "ResidualSD")
  summary <- summary[, cols, drop = FALSE]

  if (!is.na(pValue)) {
    summary <- subset(summary, `Pr(>|t|)` < pValue)
  }

  rownames(summary) <- NULL

  cat("Completed GLM summary extraction for:", variable, "\n")
  return(summary)
}

cat("=======================================================================\n")
phenotypeList <- strsplit(opt$phenotypes, ",")[[1]]

# ----------- Run and Save GLM Summaries -----------
for (pheno in phenotypeList) {
        summaryFile <- file.path(opt$outputRData, paste0(pheno,
                                                         "SummaryGLM.RData"))

        if (file.exists(summaryFile)) {
                cat("Loading existing summary for:", pheno, "\n")
                load(summaryFile)
                assign(paste0(pheno, "SummaryGLM"), fitResult,
                       envir = .GlobalEnv)
                next
        }

        cat("Running GLM summary extraction for:", pheno, "\n")

        fitObjectName <- paste0(pheno, "GLM")
        fitObject <- get(fitObjectName)

        fitResult <- cpgsGLM(
                fitList = fitObject,
                variable = pheno,
                includeResidualSD = opt$summaryResidualSD,
                pValue = opt$summaryPval,
                nCore = opt$nCores,
                libPath = opt$libPath,
                glmLibList = opt$glmLibList
        )

        save(fitResult, file = summaryFile)

        assign(paste0(pheno, "SummaryGLM"), fitResult, envir = .GlobalEnv)

        cat("Saved:", paste0(pheno, "SummaryGLM.RData"), "\n")
}
cat("=======================================================================\n")

# ----------- Save Significant CpGs -----------
saveSignificantCpGs <- function(
                resultList,
                resultName,
                baseDir = opt$significantCpGDir,
                pvalThreshold = opt$significantCpGPval,
                interactionTerm = NULL
) {
  resultDir <- file.path(baseDir, resultName)
  if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

  for (i in seq_along(resultList)) {
    coefTable <- resultList[[i]]$coef
    cpgName <- names(resultList)[i]

    if (!is.null(interactionTerm) && interactionTerm != "") {
      pattern <- paste0("^", resultName, ".*:", interactionTerm)
      matchedRows <- grep(pattern, rownames(coefTable), value = TRUE)

      if (length(matchedRows) == 0) {
        matchedRows <- grep(paste0("^", resultName),
                            rownames(coefTable), value = TRUE)
        if (length(matchedRows) > 0) {
          cat("Interaction term", interactionTerm, "dropped for CpG", cpgName,
              "- extracting main effect for", resultName, "\n")
        }
      }
    } else {
      matchedRows <- grep(paste0("^", resultName),
                          rownames(coefTable), value = TRUE)
    }

    if (length(matchedRows) > 0) {
      matchedPvals <- coefTable[matchedRows, "Pr(>|t|)"]

      if (any(matchedPvals < pvalThreshold, na.rm = TRUE)) {
        cpgDir <- file.path(resultDir, cpgName)
        if (!dir.exists(cpgDir)) dir.create(cpgDir)
        outputFile <- file.path(cpgDir, paste0(cpgName, ".txt"))
        write.table(coefTable, file = outputFile, sep = "\t", quote = FALSE)
      }
    }
  }
}

# ---------- Save Significant CpGs to Directory -----------
if (opt$saveSignificantCpGs) {
        cat("Saving significant CpGs to:", opt$significantCpGDir, "\n")

        for (pheno in strsplit(opt$phenotypes, ",")[[1]]) {
                modelObjName <- paste0(pheno, "GLM")
                if (exists(modelObjName)) {
                        resultObj <- get(modelObjName)
                        saveSignificantCpGs(
                                resultList = resultObj,
                                resultName = pheno
                        )
                }
        }

        cat("Finished saving significant CpG model outputs.\n")
}
cat("=======================================================================\n")

# ----------- Save GLM Summary Tables -----------
saveSummaryToTxt <- function(
                summaryDF,
                outputFile,
                sortByPval = TRUE
) {
        if (sortByPval) {
                if ("P.value" %in% colnames(summaryDF)) {
                        summaryDF <- summaryDF[order(summaryDF$P.value), ]
                } else if ("Pr(>|t|)" %in% colnames(summaryDF)) {
                        summaryDF <- summaryDF[order(summaryDF[["Pr(>|t|)"]]), ]
                }
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

# ----------- Save All Summaries from Phenotype List -----------

if (opt$saveTxtSummaries) {
        cat("Saving summaries to TXT files...\n")

        phenoList <- strsplit(opt$phenotypes, ",")[[1]]

        for (pheno in phenoList) {

                outputFile <- file.path(
                        opt$summaryTxtDir,
                        paste0(pheno, "SummaryGLM.txt")
                )

                summaryObjName <- paste0(pheno, "SummaryGLM")

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

# ----------- Diagnostic Plot Function -----------
diagnosticPlots <- function(
                summary,
                betaMatrix,
                variable,
                cpgPrefix = opt$cpgPrefix,
                padjmethod = opt$padjmethod,
                fdrThreshold = opt$fdrThreshold,
                outputDir = opt$outputPlots,
                plotWidth = opt$plotWidth,
                plotHeight = opt$plotHeight,
                plotDPI = opt$plotDPI
) {
        cat("Generating diagnostic plots for:", variable, "\n")

        dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

        summary$FDR <- p.adjust(summary$`Pr(>|t|)`, method = padjmethod)

        chisq <- qchisq(1 - summary$`Pr(>|t|)`, df = 1)
        lambda <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1), 3)
        message("Genomic inflation factor (λ): ", lambda)

        # -- Q-Q Plot of P-values --
        pvals <- summary$`Pr(>|t|)`
        pvals <- pvals[!is.na(pvals)]
        tiff(filename = file.path(outputDir,
                                  paste0("qqplot_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")
        qqplot(
                -log10(ppoints(length(pvals))),
                -log10(sort(pvals)),
                main = paste("Q-Q Plot of p-values for", variable,
                             "\nGenomic Inflation Factor λ =", lambda),
                xlab = "Expected -log10(p)",
                ylab = "Observed -log10(p)",
                pch = 16,
                col = "black"
        )
        abline(0, 1, col = "red")
        dev.off()

        # -- Residual SD vs log2(mean Beta) --
        cpgCols <- grep(paste0("^", cpgPrefix),
                        colnames(betaMatrix), value = TRUE)
        betaMat <- betaMatrix[, cpgCols, drop = FALSE]
        meanBeta <- colMeans(betaMat, na.rm = TRUE)
        meanBetaLog <- log2(meanBeta)
        summary$log2meanBeta <- meanBetaLog[summary$CpG]
        summary <- summary[!is.na(summary$log2meanBeta), ]

        tiff(filename = file.path(outputDir,
                                  paste0("residualSD_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")
        plot(summary$log2meanBeta,
             summary$ResidualSD,
             pch = 20,
             col = rgb(0, 0, 0, 0.4),
             xlab = "log2(Average Beta)",
             ylab = "Residual Standard Deviation",
             main = "plotSA-style: Residual SD vs Average Beta")
        lines(lowess(summary$log2meanBeta, summary$ResidualSD),
              col = "red", lwd = 2)
        dev.off()

        # -- Residual SD vs Significance --
        p <- ggplot(summary,
                    aes(x = -log10(`Pr(>|t|)`),
                        y = ResidualSD,
                        color = FDR < fdrThreshold)) +
                geom_point(alpha = 0.6) +
                geom_text_repel(data = subset(summary,
                                              FDR < fdrThreshold),
                                aes(label = CpG)) +
                scale_color_manual(values = c("FALSE" = "grey50",
                                              "TRUE" = "firebrick")) +
                labs(
                        title = paste("Residual SD vs Significance for",
                                      variable),
                        x = "-log10(p-value)",
                        y = "Residual SD",
                        color = paste("FDR <", fdrThreshold)
                ) +
                theme_minimal()

        # Save as TIFF using tiff() + print() + dev.off()
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

# ----------- Generate Diagnostic Plots -----------
cat("Generating diagnostic plots for all phenotypes...\n")
for (pheno in strsplit(opt$phenotypes, ",")[[1]]) {
        summaryName <- paste0(pheno, "SummaryGLM")
        if (exists(summaryName)) {
                diagnosticPlots(
                        summary = get(summaryName),
                        betaMatrix = phenoBT1,
                        variable = pheno,
                        fdrThreshold = opt$fdrThreshold,
                        cpgPrefix = opt$cpgPrefix,
                        padjmethod = opt$padjmethod,
                        outputDir = opt$outputPlots,
                        plotWidth = opt$plotWidth,
                        plotHeight = opt$plotHeight,
                        plotDPI = opt$plotDPI
                )
        }
}
cat("Finished generating all diagnostic plots.\n")
cat("Plots saved to:", opt$outputPlots, "\n")
cat("=======================================================================\n")

# ----------- Load Annotation Data -----------
cat("Loading annotation object:", opt$annotationPackage, "\n")
annotationObject <- getAnnotation(get(opt$annotationPackage))

cat("Annotation loaded with", nrow(annotationObject), "probes\n")
cat("Annotation columns available:\n")
print(colnames(annotationObject))
cat("=======================================================================\n")

# ----------- Apply Annotation -----------
annotateGLM <- function(
                summaryList,
                annotationObject,
                annotationCols = strsplit(opt$annotationCols, ",")[[1]]
) {
        modelNames <- names(summaryList)

        cat("Merging GLM summaries...\n")
        cleanedSummaries <- lapply(modelNames, function(modelName) {
          df <- summaryList[[modelName]]

          coefNames <- unique(df$Coefficient)
          dfList <- lapply(coefNames, function(coefName) {
            subdf <- df[df$Coefficient == coefName, ]
            subdf[[paste0(coefName, "P.Value")]] <- subdf$`Pr(>|t|)`
            subdf <- subdf[, c("CpG", paste0(coefName, "P.Value"))]
            return(subdf)
          })

        # Merge multiple coefficient data.frames by CpG
        mergedCoef <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE),
                             dfList)
        return(mergedCoef)

        })

        mergedSummary <- Reduce(function(x, y) merge(x,
                                                     y,
                                                     by = "CpG",
                                                     all = TRUE),
                                cleanedSummaries)

        annDF <- as.data.frame(annotationObject)
        annDF$CpG <- rownames(annDF)

        cat("Merging annotation with GLM summary...\n")
        annotatedResults <- merge(mergedSummary,
                                  annDF, by = "CpG", all.x = TRUE)

        finalCols <- c("CpG",
                       unlist(lapply(cleanedSummaries, function(df) colnames(df)[-1])),
                       annotationCols)

        annotatedResults <- annotatedResults[, finalCols]
        colnames(annotatedResults)[1] <- "IlmnID"

        cat("Annotation completed. Annotated CpGs:",
            nrow(annotatedResults), "\n")
        return(annotatedResults)
}

# ----------- Execute Annotation and Save Results -----------

cat("Running annotation of GLM summary results...\n")

# Split phenotype names and fetch each corresponding summary object
phenotypeNames <- strsplit(opt$phenotypes, ",")[[1]]
summaryList <- setNames(
        lapply(phenotypeNames, function(pheno) {
                get(paste0(pheno, "SummaryGLM"))
        }),
        phenotypeNames
)

# Convert annotationCols from comma-separated string to vector
annotationColsVec <- strsplit(opt$annotationCols, ",")[[1]]

# Run the annotation function
annotatedGLM <- annotateGLM(
        summaryList = summaryList,
        annotationObject = getAnnotation(get(opt$annotationPackage)),
        annotationCols = annotationColsVec
)

# Save annotated results
cat("Saving annotated GLM summary to:", opt$annotatedGLMOut, "\n")
write.csv(
        annotatedGLM,
        file = file.path(opt$annotatedGLMOut, "annotatedGLM.csv"),
        row.names = FALSE)
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)

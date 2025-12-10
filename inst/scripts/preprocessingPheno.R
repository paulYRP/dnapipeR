# !/usr/bin/env Rscript
# ==============================================================================
# DNAm Phenotype Preprocessing Script
# Script Name: preprocessingPheno.R
# Description: Merges phenotype data with EWAS QC, subsets by timepoints,
#              and prepares beta/m matrices for downstream modeling.
# ==============================================================================
# Usage Example (Full version):
# ==============================================================================
# Rscript preprocessingPheno.R \
#   --phenoFile data/preprocessingMinfiEwasWater/pheno.csv \
#   --betaPath rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
#   --mPath rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
#   --cnPath rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
#   --timepoints 1,2,3 \
#   --combineTimepoints 1,3 \
#   --outputPheno data/preprocessingPheno/ \
#   --outputRData rData/preprocessingPheno \
#   --outputLogs logs/preprocessingPheno

# ==============================================================================
# Usage Example (Default values with key parameters):
# ==============================================================================
# Rscript preprocessingPheno.R \
#   --phenoFile data/preprocessingMinfiEwasWater/pheno.csv \
#   --phenoEWAS data/preprocessingEwastools/pheno_ewasQC.csv

# ==============================================================================
# DNAm Preprocessing Script — Input Arguments
# ==============================================================================
#   --phenoFile              [FILE]   Path to primary phenotype CSV/TSV file
#   --phenoEWAS              [FILE]   Path to EWAS-based QC phenotype file
#   --betaPath               [FILE]   RData file containing beta matrix
#   --mPath                  [FILE]   RData file containing M-value matrix
#   --cnPath                 [FILE]   RData file containing CN matrix
#   --dropColumnsPhenoEWAS   [STR]    Columns to drop from phenoEWAS (comma-separated)
#   --colsToRenamePhenoEWAS  [STR]    Columns in phenoEWAS to append '.EWAS' suffix (comma-separated)
#   --mergeKey               [STR]    Column name to merge pheno and phenoEWAS (default: "SID")
#   --factorVars             [STR]    Variables to convert to factors (comma-separated)
#   --factorPrefixes         [STR]    Prefixes for factor levels (must match factorVars order)
#   --timepoints             [STR]    Comma-separated list of timepoints to split/subset (default: "1,2,3")
#   --combineTimepoints      [STR]    Timepoints to combine for longitudinal analysis (default: "1,3")
#   --outputPheno            [DIR]   Output CSV file for merged phenotype (default: "data/preprocessingPheno/")
#   --outputRData            [DIR]    Output folder for timepoint-subset and merged RData (default: "rData/preprocessingPheno")
#   --outputLogs             [DIR]    Folder to save logs (default: "logs/preprocessingPheno")
# ==============================================================================


# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(optparse)
        library(dplyr)
        library(tidyverse)
        library(readr)
        library(dplyr)
        library(stringr)
        library(purrr)
})

# ----------- Define Input Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--phenoFile", default = "data/preprocessingMinfiEwasWater/phenoLC.csv", help = "Input phenotype CSV file [default: %default]", metavar = "FILE"),
        make_option("--sepType", default = "" ,  type = "character", help = "Separator for phenotype file (e.g., ',' or '\\t')[default: %default] means NULL"),
        make_option("--betaPath", default = "rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData", help = "Path to Beta matrix RData [default: %default]", metavar = "FILE"),
        make_option("--mPath", default = "rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData", help = "Path to M-values RData [default: %default]", metavar = "FILE"),
        make_option("--cnPath", default = "rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData", help = "Path to CN matrix RData [default: %default]", metavar = "FILE"),
        make_option("--SampleID", type = "character", default = "SampleName", help = "Column name to use as sample IDs"),
        make_option("--timeVar", type = "character", default = "Timepoint", help = "Column name that identifies timepoints (e.g., 'Timepoint', 'Visit', 'Wave')"),
        make_option("--timepoints", default = "1,2", help = "Timepoints to subset for beta and M values", metavar = "T1,T2,T3"),
        make_option("--combineTimepoints", default = "1,2", help = "Timepoints to combine for longitudinal analysis", metavar = "T1,Tn"),
        make_option("--outputPheno", default = "data/preprocessingPheno/", help = "Path to save final phenotype CSVs", metavar = "DIR"),
        make_option("--outputRData", default = "rData/preprocessingPheno/metrics", help = "Directory to save processed RData objects metrics", metavar = "DIR"),
        make_option("--outputRDataMerge", default = "rData/preprocessingPheno/mergeData", help = "Directory to save processed RData objects mergedata", metavar = "DIR"),
        make_option("--sexColumn", type = "character", default = "Sex", help = "Column in sample data with sample sex (e.g., 'Sex', coded F/M)"),
        make_option("--outputLogs", default = "logs/", help = "Directory for all log output [default: %default]", metavar = "DIR"),
        make_option("--outputDir", default = "data/preprocessingPheno", help = "Directory to save Beta CSV and ZIP outputs [default: %default]", metavar = "DIR")
)))

dir.create(opt$outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputRDataMerge, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputPheno, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputDir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================

# ----------- Logging Setup -----------
logFilePath <- file.path(opt$outputLogs, "log_preprocessingPheno.txt")
logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
#===============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting Phenotype Preprocessing ====\n")
cat("Start Time:               ", format(Sys.time()), "\n")
cat("Log file path:            ", logFilePath, "\n\n")
cat("Phenotype file:           ", opt$phenoFile, "\n")
cat("Beta path:                ", opt$betaPath, "\n")
cat("M-values path:            ", opt$mPath, "\n")
cat("CN path:                  ", opt$cnPath, "\n\n")

cat("Identifier column:        ", opt$SampleID, "\n")
cat("Timepoint column:        ", opt$timeVar, "\n")
cat("Timepoints (if present):  ", opt$timepoints, "\n")
cat("Combine timepoints:       ", opt$combineTimepoints, "\n\n")
cat("Sex column:               ", opt$sexColumn, "\n")

cat("Output phenotype dir:     ", opt$outputPheno, "\n")
cat("RData metrics dir:        ", opt$outputRData, "\n")
cat("RData merge dir:          ", opt$outputRDataMerge, "\n\n")
cat("=======================================================================\n")

# ----------- Load Data -----------
load(opt$betaPath)
load(opt$mPath)
load(opt$cnPath)

cat("Beta dimensions: ", dim(beta), "\n")
cat("M dimensions: ", dim(m), "\n")
cat("CN dimensions: ", dim(cn), "\n")

# ----------- Read Phenotype File -----------
if (opt$sepType == "\\t") {
  sepChar <- "\t"
} else if (opt$sepType == "") {
  sepChar <- NULL
} else {
  sepChar <- opt$sepType
}

# Now read the phenotype file
if (!is.null(sepChar)) {
  pheno <- read.csv(opt$phenoFile, sep = sepChar)
} else {
  pheno <- read.csv(opt$phenoFile)
}

cat("Phenotype file loaded with",
    nrow(pheno), "samples and", ncol(pheno), "columns.\n")
cat("Preview of phenoLC:\n")
print(head(pheno[, 1:5]))
cat("=======================================================================\n")

# ----------- Subsetting Timepoints & Data Splitting -----------
timepoints <- as.numeric(strsplit(opt$timepoints, ",")[[1]])
cat("Subsetting to timepoints:", paste(timepoints, collapse = ", "), "\n")

# Print available timepoints in the phenotype
cat("Available values in", opt$timeVar, "column:\n")
print(table(pheno[[opt$timeVar]], useNA = "ifany"))

for (tp in timepoints) {
  # Subset phenotype by Timepoint
  phenoSub <- subset(pheno, pheno[[opt$timeVar]] == tp)
  assign(paste0("phenoT", tp), phenoSub)

  # Subset matrices using SID (SampleID) from phenoSub
  sids <- as.character(phenoSub[[opt$SampleID]])

  assign(paste0("betaT", tp), beta[, sids])
  assign(paste0("mT", tp),    m[,    sids])
  assign(paste0("cnT", tp),   cn[,   sids])
}

# Save each subset
for (tp in timepoints) {
        write.csv(get(paste0("phenoT", tp)), file = file.path(opt$outputPheno,
                                                              paste0("phenoT",
                                                                     tp, ".csv")),
                  row.names = FALSE)
        save(list = paste0("betaT", tp), file = file.path(opt$outputRData,
                                                          paste0("betaT",
                                                                 tp, ".RData")))
        save(list = paste0("mT", tp), file = file.path(opt$outputRData,
                                                       paste0("mT",
                                                              tp, ".RData")))
}

# ----------- Merge Combined Timepoints for Longitudinal Analysis -----------
cat("Combining timepoints:", opt$combineTimepoints, "\n")
combineTPs <- as.numeric(strsplit(opt$combineTimepoints, ",")[[1]])

combinedPhenoList <- lapply(combineTPs, function(tp) get(paste0("phenoT", tp)))
phenoCombined <- do.call(rbind, combinedPhenoList)

combineSuffix <- paste0("T", paste(combineTPs, collapse = "T"))
write.csv(phenoCombined,
          file = file.path(opt$outputPheno,
                                          paste0("pheno",
                                                 combineSuffix, ".csv")),
          row.names = FALSE)

cat("Saved combined phenotype file for T1T2 at:", opt$outputPheno, "\n")
cat("=======================================================================\n")

# ----------- Merge Beta Matrix with Phenotype ----------
mergeBeta <- function(phenoFrame, betaMatrix, id = opt$SampleID) {
        rownames(phenoFrame) <- phenoFrame[[id]]
        matched <- intersect(rownames(phenoFrame), colnames(betaMatrix))
        phenoFrame <- phenoFrame[matched, ]
        betaMatrix <- betaMatrix[, matched]

        betaTranp <- as.data.frame(t(betaMatrix))
        mergedData <- cbind(phenoFrame, betaTranp)
        return(mergedData)
}

# Perform merge for each timepoint

mergedList <- list()
for (tp in timepoints) {
        cat("Processing merge for timepoint:", tp, "\n")

        phenoObj <- paste0("phenoT", tp)
        betaObj <- paste0("betaT", tp)

        if (!exists(phenoObj) || !exists(betaObj)) {
                cat("Warning: One or both objects not found for T", tp, "\n", sep = "")
                next
        }

        phenoTemp <- get(phenoObj)
        betaTemp <- get(betaObj)

        cat("  - pheno rows:", nrow(phenoTemp), "\n")
        cat("  - beta cols:", ncol(betaTemp), "\n")

        mergedTemp <- tryCatch({
                mergeBeta(phenoTemp, betaTemp)
        }, error = function(e) {
                cat("[ERROR] mergeBeta failed for timepoint", tp, ":\n", conditionMessage(e), "\n")
                return(NULL)
        })

        if (!is.null(mergedTemp)) {
                mergedList[[as.character(tp)]] <- mergedTemp
                save(mergedTemp, file = file.path(opt$outputRDataMerge,
                                                  paste0("phenoBetaT", tp, ".RData")))
                cat("Saved merged object for T", tp, "\n", sep = "")
        } else {
                cat("Skipping save for T", tp, " due to error\n")
        }
}


# ----------- Combine merged phenotype + beta matrix -----------
combined <- do.call(rbind, mergedList[as.character(combineTPs)])
save(combined,
     file = file.path(opt$outputRDataMerge,
                      paste0("phenoBeta", combineSuffix, ".RData")))

cat("Combined data saved for timepoints:",
    paste(combineTPs, collapse = ", "), "\n")

cat("Merged data saved to:", opt$outputRDataMerge, "\n")
cat("=======================================================================\n")
# ----------- Preprocessing Betas for Horvath Calculator -----------

betaCSV <- as.data.frame(beta)
betaCSV <- tibble::rownames_to_column(betaCSV, var = "ProbeID")

# Inspect changes
dim(betaCSV)
print(head(betaCSV)[1:5, 1:5])

betaCSVPath <- file.path(opt$outputDir, "beta.csv")
write.csv(betaCSV, file = betaCSVPath, row.names = FALSE)
cat("Beta CSV file for ClockFundation saved to:", betaCSVPath, "\n")

zipFile <- file.path(opt$outputDir, "beta.zip")

zip(zipfile = zipFile, files = betaCSVPath, flags = "-j")

cat("Beta ZIP file for ClockFundation saved to:", zipFile, "\n")

# ----------- Preprocessing CSV for Horvath Calculator -----------

# Rename the column "SampleName" to "id"
pheno <- pheno %>%
  rename(id = opt$SampleID)

# Recode Sex only if values are not already "Male"/"Female"
uniqueSex <- unique(pheno[[opt$sexColumn]])

if (!all(uniqueSex %in% c("Male", "Female"))) {
  cat("Re-encoding Sex: 0 = Female, 1 = Male\n")
  pheno[[opt$sexColumn]] <- ifelse(pheno[[opt$sexColumn]] == 0, "Female", "Male")
} else {
  cat("Sex column already contains 'Male' and 'Female'. Skipping recoding.\n")
}

phenoCSVPath <- file.path(opt$outputDir, "phenoCF.csv")

write.csv(pheno, file = phenoCSVPath, row.names = FALSE)
cat("Sample file for ClockFundation saved to:", phenoCSVPath, "\n")

cat("=======================================================================\n")

# ----------- Close Logging -----------
cat("\nSession Info:\n")
print(sessionInfo())
# =============================================================================
sink(type = "message")
sink()
close(logCon)

#!/usr/bin/env Rscript
# ==============================================================================
# DNAm Surrogate Variable Estimation Script (Enmix-based)
# Script Name: svaEnmix.R
# Description: Derives surrogate variables from control probes using ewastools::ctrlsva
#              and performs exploratory analysis (Sentrix ID/Position effects).
# ==============================================================================

# Usage Example (Full version):
# Rscript svaEnmix.R \
#   --phenoFile data/preprocessing/pheno.csv \
#   --idatFolder data/idats/ \
#   --outputLogs logs/svaEnmix/ \
#   --nSamples 50 \
#   --SampleID SampleID \
#   --scriptLabel svaEnmix \
#   --sepType "," \
#   --SentrixIDColumn SentrixID \
#   --SentrixPositionColumn SentrixPosition \
#   --ctrlSvaPercVar 0.90 \
#   --ctrlSvaFlag 1 \
#   --qcTiffPath figures/svaEnmix/sva/ \
#   --tiffWidth 2000 --tiffHeight 1000 --tiffRes 150
# ==============================================================================

# Usage Example (Minimal):
# Rscript svaEnmix.R \
#   --phenoFile data/pheno.csv \
#   --idatFolder data/idats/ \
#   --SampleID SampleID
# ==============================================================================

# ----------- Libraries -----------
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICv2manifest)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  library(Gviz)
  library(ENmix)
  library(ewastools)
  library(optparse)
  library(ggplot2)
  library(ggpubr)
  library(MASS)
})

# ==============================================================================
# DNAm Surrogate Variable Estimation — Input Arguments

#   --phenoFile              [FILE]   Path to phenotype/sample sheet (.csv or .tsv)
#   --idatFolder             [DIR]    Path to raw .idat files
#   --outputLogs             [DIR]    Directory to save logs and execution metadata
#   --nSamples               [INT]    Number of samples to use (for testing); NA = all
#   --SampleID               [STR]    Column name to assign sample names (e.g., SampleID)

#   --scriptLabel            [STR]    Script label for tagging folders/outputs (e.g., "svaEnmix")

#   --sepType                [STR]    Field separator in phenotype file ("," or "\\t")
#   --SentrixIDColumn        [STR]    Column name for Sentrix IDs (e.g., "SentrixID")
#   --SentrixPositionColumn  [STR]    Column name for Sentrix positions (e.g., "SentrixPosition")

#   --ctrlSvaPercVar         [NUM]    Percentage of variance explained by SVs (default = 0.90)
#   --ctrlSvaFlag            [INT]    Flag option for ctrlsva() (default = 1)

#   --qcTiffPath             [DIR]    Path for output TIFF figures
#   --tiffWidth              [INT]    Width of TIFF images (e.g., 2000)
#   --tiffHeight             [INT]    Height of TIFF images (e.g., 1000)
#   --tiffRes                [INT]    TIFF resolution in DPI (e.g., 150)
# ==============================================================================

# ----------- Command Line Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
  make_option("--phenoFile", default = "data/preprocessingMinfiEwasWater/phenoLC.csv", help = "Path to phenotype CSV file", metavar = "FILE"),
  make_option("--rgsetData", type = "character", default = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData", help = "Path to the RGSet RData object [default = %default]"),
  make_option("--sepType", default = "" ,  type = "character", help = "Separator for phenotype file (e.g., ',' or '\\t')"),
  make_option("--outputLogs", default = "logs/", help = "Directory for all output", metavar = "DIR"),
  make_option("--nSamples", type = "integer", default = NA, help = "Limit to first N samples [default: all]"),
  make_option("--SampleID", type = "character", default = "SampleName", help = "Column name to use as sample IDs"),
  make_option("--arrayType", default = "IlluminaHumanMethylationEPICv2", help = "Array platform name"),
  make_option("--annotationVersion", default = "20a1.hg38", help = "Annotation version"),
  make_option("--SentrixIDColumn", type = "character", default = "SentrixID", help = "Column name for Sentrix ID used in coloring SVA plots"),
  make_option("--SentrixPositionColumn", type = "character", default = "SentrixPosition", help = "Column name for Sentrix Position used in shaping SVA plots"),
  make_option("--ctrlSvaPercVar", type = "double", default = 0.90, help = "Percentage variance explained for ctrlSVA [default: 0.90]"),
  make_option("--ctrlSvaFlag", type = "integer", default = 1, help = "ctrlSVA flag parameter [default: 1]"),
  make_option("--scriptLabel", default = "svaEnmix", help = "Label for output folders/logs"),
  make_option("--tiffWidth", type = "integer", default = 2000),
  make_option("--tiffHeight", type = "integer", default = 1000),
  make_option("--tiffRes", type = "integer", default = 150)

)))

# ==============================================================================

# ----------- Logging Setup -----------
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)

logFilePath <- file.path(opt$outputLogs,"log_svaEnmix.txt")
logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting SVA Estimation with Enmix ====\n")
cat("Start time: ", format(Sys.time()), "\n\n")
cat("Log file path: ", logFilePath, "\n\n")
cat("Pheno file: ", opt$phenoFile, "\n")
cat("Separator type: ", ifelse(is.null(opt$sepType), "default, No separator", opt$sepType), "\n")
cat("IDAT folder: ", opt$idatFolder, "\n")
cat("Log directory: ", opt$outputLogs, "\n")
cat("Sample limit: ", ifelse(is.na(opt$nSamples), "All", opt$nSamples), "\n")
cat("SampleID column: ", opt$SampleID, "\n")
cat("Sentrix ID column: ", opt$SentrixIDColumn, "\n")
cat("Sentrix Position column: ", opt$SentrixPositionColumn, "\n")
cat("Script label: ", opt$scriptLabel, "\n")
cat("ctrlSva percvar: ", opt$ctrlSvaPercVar, "\n")
cat("ctrlSva flag: ", opt$ctrlSvaFlag, "\n")
cat("TIFF dimensions (WxH): ", opt$tiffWidth, "x", opt$tiffHeight, " at", opt$tiffRes, "dpi\n")
# =============================================================================

# ----------- Directory Setup for Figures -----------
dir.create(file.path("figures", opt$scriptLabel), showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path("data", opt$scriptLabel), showWarnings = FALSE,
           recursive = TRUE)
cat("=======================================================================\n")

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
  targets <- read.csv(opt$phenoFile, sep = sepChar)
} else {
  targets <- read.csv(opt$phenoFile)
}

if (!is.na(opt$nSamples) && opt$nSamples < nrow(targets)) {
  targets <- targets[1:opt$nSamples, ]
  cat("Subsetting to", opt$nSamples, "samples for testing.\n")
} else {
  cat("Using all", nrow(targets), "samples.\n")
}

cat("Phenotype file loaded with",
    nrow(targets), "samples and", ncol(targets), "columns.\n")
cat("Preview of targets:\n")
print(head(targets[, 1:6]))
cat("=======================================================================\n")

# ----------- Load IDAT Files into RGSet -----------
load(opt$rgsetData)

# Assign custom sample names
sampleNames(RGSet) <- targets[[opt$SampleID]]
cat("RGSet loaded with", length(sampleNames(RGSet)), "samples.\n")
cat("=======================================================================\n")

# ----------- Estimate Surrogate Variables from Control Probes -----------
sva <- ctrlsva(
  rgSet = RGSet,
  percvar = opt$ctrlSvaPercVar,
  flag = opt$ctrlSvaFlag
)
cat("Surrogate variables matrix (first few rows):\n")
print(head(sva))

# ----------- Plot SVA Colored by SentrixID -----------
sentrixID <- as.factor(pData(RGSet)[[opt$SentrixIDColumn]])

# Create TIFF output
svaSentrixPath <- file.path("figures",
                        opt$scriptLabel, "sva_SentrixID.tiff")
tiff(file = svaSentrixPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")

plot(sva[, 1], sva[, 2],
     col = rainbow(length(levels(sentrixID)))[sentrixID],
     pch = 16,
     xlab = "Surrogate Variable 1 (PC1)",
     ylab = "Surrogate Variable 2 (PC2)",
     main = "Surrogate Variables Colored by Chip (SentrixID)")
legend("topright", legend = levels(sentrixID),
       col = rainbow(length(levels(sentrixID))),
       pch = 16, title = "SentrixID", cex = 0.6)

dev.off()

cat("SVA Sentrix plot saved to: ", svaSentrixPath, "\n")

# ----------- Plot SVA Colored by SentrixPosition -----------
sentrixPos <- as.factor(pData(RGSet)[[opt$SentrixPositionColumn]])

svaPositionpath <- file.path("figures", opt$scriptLabel, "sva_SentrixPosition.tiff")

tiff(file = svaPositionpath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")

plot(sva[, 1], sva[, 2],
     col = rainbow(length(levels(sentrixPos)))[sentrixPos],
     pch = 16,
     xlab = "Surrogate Variable 1 (PC1)",
     ylab = "Surrogate Variable 2 (PC2)",
     main = "Surrogate Variables Colored by Sentrix Position")
legend("topright", legend = levels(sentrixPos),
       col = rainbow(length(levels(sentrixPos))),
       pch = 16, title = "SentrixPosition", cex = 0.6)

dev.off()

cat("SVA Position plot saved to: ", svaPositionpath, "\n")
cat("=======================================================================\n")

# ----------- Linear Models for Surrogate Variables (ANOVA) -----------
K <- ncol(sva)
cat("Number of surrogate variables (K):", K, "\n")

# Print class and levels of SentrixID
cat("SentrixID class:", class(sentrixID), "\n")
cat("SentrixID unique levels:", length(unique(sentrixID)), "\n")
print(table(sentrixID))

# Print class and levels of SentrixPosition
cat("SentrixPosition class:", class(sentrixPos), "\n")
cat("SentrixPosition unique levels:", length(unique(sentrixPos)), "\n")
print(table(sentrixPos))

# Print example row of SVA matrix
cat("First row of SVA matrix:\n")
print(sva[1, ])

# Confirm if sample names align
cat("Sample names in SVA matrix:",
    paste(rownames(sva)[1:5], collapse = ", "), "\n")
cat("Sample names in pData(RGSet):",
    paste(rownames(pData(RGSet))[1:5], collapse = ", "), "\n")

# Fit linear models for each surrogate variable
lmsvaFull <- lapply(1:K, function(i)
  lm(sva[, i] ~ SentrixID + SentrixPosition,
     data.frame("SentrixID" = sentrixID,
                "SentrixPosition" = sentrixPos))
)

lmsvaRed <- vector("list", K)

# ----------- Save summaries of full models -----------
capture.output(summary(lmsvaFull[[1]]),
               file = file.path("data", opt$scriptLabel, "summary_full_sva1.txt"))

if (K >= 2) {
  capture.output(summary(lmsvaFull[[2]]),
                 file = file.path("data", opt$scriptLabel, "summary_full_sva2.txt"))
}

# Perform backward elimination and write ANOVA output
for(i in 1:K){
  lmtmp = lmsvaFull[[i]]
  while(1){
    dttmp = dropterm(lmtmp, test = "F")
    if(max(dttmp$`Pr(F)`, na.rm = TRUE) > (0.05))
      ttmp = rownames(dttmp)[which.max(dttmp$`Pr(F)`)]
    else break
    lmtmp = update(lmtmp, paste(".~. - ", ttmp) )
    capture.output(dttmp,
                   file = file.path("data",
                                    opt$scriptLabel, paste0("dropterm_step_sva", i, ".txt")),
                   append = TRUE)
    capture.output(summary(lmtmp),
                   file = file.path("data",
                                    opt$scriptLabel, paste0("dropterm_model_sva", i, ".txt")),
                   append = TRUE)
  }

  lmsvaRed[[i]] = lmtmp
}

# ----------- Save ANOVA summaries for full and reduced models -----------
for (i in 1:K) {
  capture.output(anova(lmsvaFull[[i]]),
                 file = file.path("data",
                                  opt$scriptLabel, paste0("anova_full_sva", i, ".txt")))

  capture.output(anova(lmsvaRed[[i]]),
                 file = file.path("data",
                                  opt$scriptLabel, paste0("anova_reduced_sva", i, ".txt")))
}
cat("=======================================================================\n")

# ----------- Plot Matrix of Surrogate Variables Colored by SentrixID and Shape by Position -----------

# Prepare output TIFF file

svaSentrixPositionPath <- file.path("figures",
                                    opt$scriptLabel, "sva_SentrixIDPosition.tiff")

tiff(file = svaSentrixPositionPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes,
     type = "cairo")

# Prepare plotting layout
par(mfrow = c(K, K), family = "Times", las = 1)

# Extract and map IDs/positions
colorMap <- rainbow(length(levels(sentrixID)))
pchMap <- 1:length(levels(sentrixPos))

# Plot matrix
for (i in 1:K) {
  for (j in 1:K) {
    plot(sva[, j], sva[, i],
         col = colorMap[sentrixID],
         pch = pchMap[sentrixPos],
         xlab = paste("SV", j),
         ylab = paste("SV", i),
         main = "Effects of Sentrix ID (color) & Sentrix Position (shape)")

    # Legend only in top-left panel
    if (i == 1 && j == 1) {
      legend("topright",
             legend = levels(sentrixID),
             col = colorMap,
             pch = 15,
             title = "SentrixID",
             cex = 0.6)
      legend("bottomright",
             legend = levels(sentrixPos),
             pch = pchMap,
             title = "SentrixPosition",
             cex = 0.6)
    }
  }
}

# Close plotting device
dev.off()

cat("SVA Sentrix/Position plot saved to: ", svaSentrixPositionPath, "\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)

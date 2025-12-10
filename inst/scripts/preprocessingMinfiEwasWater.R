#!/usr/bin/env Rscript
# ==============================================================================
# DNAm Preprocessing Script (Minfi-based)
# Script Name: preprocessingMinfiEwasWater.R
# Description: Preprocesses DNAm data using Minfi package. 
# ==============================================================================

# Usage Example (Full version):

# Rscript preprocessingMinfiEwasWater.R \
#   --phenoFile data/preprocessingMinfiEwasWater/pheno.csv \
#   --idatFolder data/idats/ \
#   --outputLogs logs/preprocessingMinfiEwasWater/ \
#   --nSamples 100 \
#   --idColumns SID,Timepoint \
#   --arrayType IlluminaHumanMethylationEPICv2 \
#   --annotationVersion 20a1.hg38 \
#   --scriptLabel preprocessingMinfiEwasWater \
#   --baseDataFolder rData \
#   --qcTiffPath figures/preprocessingMinfiEwasWater/quality_control_MSet.tiff \
#   --tiffWidth 2000 \
#   --tiffHeight 1000 \
#   --tiffRes 150 \
#   --qcCutoff 10.5 \
#   --detPtype m+u \
#   --densityTiffPath figures/preprocessingMinfiEwasWater/densityBeta_MSet.tiff \
#   --pdfReportPath reports/qc_report_RGSet.pdf \
#   --funnormSeed 123 \
#   --normMethods "funnorm;quantile" \
#   --pvalThreshold 0.01 \
#   --chrToRemove chrX,chrY \
#   --snpsToRemove SBE,CpG \
#   --mafThreshold 0.5 \
#   --crossReactivePath data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv \
#   --plotGroupVar Ethnicity \
#   --betaMPlotPath figures/preprocessingMinfiEwasWater/densityBetaM_MSetF_Flt_Rxy_Ds_Rc.tiff
# ==============================================================================

# Usage Example (Default values with key parameters):

# Rscript preprocessingMinfiEwasWater.R \
#   --phenoFile data/preprocessingMinfiEwasWater/pheno.csv \
#   --idatFolder data/idats/ \
#   --nSamples 10 \ ## Testing with 10, NA for all
#   --pvalThreshold 0.01 \
#   --mafThreshold 0.5 \
#   --crossReactivePath data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv \
# ==============================================================================

# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(RColorBrewer)
        library(data.table)
        library(wateRmelon)
        library(minfi)
        library(IlluminaHumanMethylationEPICv2manifest)
        library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
        library(Gviz)
        library(ENmix)
        library(BeadSorted.Saliva.EPIC)
        library(ewastools)
        library(optparse)
        library(ggplot2)
        library(ggpubr)
        library(MASS)
})
# ==============================================================================

# DNAm Preprocessing Script — Input Arguments

#   --phenoFile            [FILE]   Path to phenotype/sample sheet (.csv or .tsv)
#   --idatFolder           [DIR]    Path to raw .idat files
#   --outputLogs           [DIR]    Directory to save logs and execution metadata
#   --nSamples             [INT]    Number of samples to use (for testing); NA means all
#   --idColumns            [STR]    Comma-separated column names to build sample IDs (e.g., "SID,Timepoint")
#   --arrayType            [STR]    Array type for annotation() (e.g., "IlluminaHumanMethylationEPICv2")
#   --annotationVersion    [STR]    Annotation string for array (e.g., "20a1.hg38")
#   --scriptLabel          [STR]    Script label for tagging folders/outputs (e.g., "preprocessingMinfiEwasWater")
#   --baseDataFolder       [DIR]    Root folder to store .RData outputs (default: "rData")
#
#   --qcTiffPath           [FILE]   Output TIFF file path for sample QC plot
#   --tiffWidth            [INT]    Width of TIFF images (e.g., 2000)
#   --tiffHeight           [INT]    Height of TIFF images (e.g., 1000)
#   --tiffRes              [INT]    Resolution (dpi) for TIFF images (e.g., 150)
#
#   --qcCutoff             [NUM]    Quality control cutoff value for bad sample detection
#   --detPtype             [STR]    Detection p-value type used ("m+u", "both", etc.)
#   --densityTiffPath      [FILE]   Output TIFF for beta density plot from MSet
#   --pdfReportPath        [FILE]   Path to save PDF report from qcReport()
#
#   --funnormSeed          [INT]    Random seed for normalization
#   --normMethods          [STR]    Normalization method(s), separated by ";" (e.g., "funnorm;swan")
#
#   --pvalThreshold        [NUM]    Detection p-value cutoff for probe filtering (e.g., 0.01)
#   --chrToRemove          [STR]    Comma-separated list of chromosomes to remove (e.g., "chrX,chrY")
#   --snpsToRemove         [STR]    SNP categories to remove (e.g., "SBE,CpG")
#   --mafThreshold         [NUM]    Minor allele frequency threshold (e.g., 0.5)
#
#   --crossReactivePath    [FILE]   Path to file with list of cross-reactive probes to remove
#   --plotGroupVar         [STR]    Phenotype variable used for coloring density plots (e.g., "Ethnicity")
#   --betaMPlotPath        [FILE]   Output TIFF for Beta & M-value density plots after filtering
# ==============================================================================

# ----------- Command Line Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--phenoFile", default = "data/preprocessingMinfiEwasWater/pheno.csv", help = "Path to phenotype CSV file", metavar = "FILE"),
        make_option("--sepType", default = "" ,  type = "character", help = "Separator for phenotype file (e.g., ',' or '\\t')"),
        make_option("--idatFolder", default = "data/preprocessingMinfiEwasWater/idats/", help = "Folder with IDAT files", metavar = "DIR"),
        make_option("--outputLogs", default = "logs/", help = "Directory for all output", metavar = "DIR"),
        make_option("--nSamples", type = "integer", default = NA, help = "Limit to first N samples [default: all]"),
        make_option("--SampleID", type = "character", default = "SampleName", help = "Column name to use as sample IDs"),
        make_option("--arrayType", default = "IlluminaHumanMethylationEPICv2", help = "Array platform name"),
        make_option("--annotationVersion", default = "20a1.hg38", help = "Annotation version"),
        make_option("--scriptLabel", default = "preprocessingMinfiEwasWater", help = "Label for output folders/logs"),
        make_option("--baseDataFolder", default = "rData", help = "Base folder for RData output"),
        make_option("--tiffWidth", type = "integer", default = 2000),
        make_option("--tiffHeight", type = "integer", default = 1000),
        make_option("--tiffRes", type = "integer", default = 150),
        make_option("--qcCutoff", type = "double", default = 10.5),
        make_option("--detPtype", default = "m+u", help = "Detection P-value type"),
        make_option("--detPThreshold", type = "double", default = 0.05, help = "Threshold for mean detection p-value to remove bad samples [default %default]"),
        make_option("--funnormSeed", type = "integer", default = 123, help = "Seed for normalization"),
        make_option("--normMethods", default = "adjustedfunnorm", help = "Normalization methods separated by ; (e.g., funnorm;swan)"),
        make_option("--sexColumn", type = "character", default = "Sex", help = "Column in phenotype (pData) with sample sex (e.g., 'Sex', coded F/M)"),
        make_option("--pvalThreshold", type = "double", default = 0.01),
        make_option("--chrToRemove", default = "chrX,chrY", help = "Chromosomes to remove"),
        make_option("--snpsToRemove", default = "SBE,CpG", help = "SNP positions to filter"),
        make_option("--mafThreshold", type = "double", default = 0.1),
        make_option("--crossReactivePath", default = "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv", type = "character", help = "Path to cross-reactive probe file"),
        make_option("--plotGroupVar", default = "Sex", help = "Grouping variable for density plots"),
        make_option("--lcRef", default = "salivaEPIC", help = "Reference for estimateLC (e.g., salivaEPIC, saliva, Reinius+Lin ...)"),
        make_option("--phenoOrder", default = "SampleName;Timepoint;Sex;PredSex;Basename;SentrixID;SentrixPosition", help = "Semicolon-separated leading column order; others appended"),
        make_option("--lcPhenoDir", default = "data/preprocessingMinfiEwasWater", help = "Output directory for the phenoLC.csv file [default: %default]")
        
        
)))

# Split comma/semicolon lists
opt$chrToRemoveList <- strsplit(opt$chrToRemove, ",")[[1]]
opt$snpList         <- strsplit(opt$snpsToRemove, ",")[[1]]
opt$normMethodList  <- strsplit(opt$normMethods, ";")[[1]]
# ==============================================================================

# ----------- Logging Setup -----------
logFilePath <- file.path(opt$outputLogs, "log_preprocessingMinfiEwasWater.txt")
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)

logCon <- file(logFilePath, open = "wt")  

sink(logCon, split = TRUE)                     
sink(logCon, type = "message")                
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting", opt$scriptLabel, "====\n")
cat("Start Time:               ", format(Sys.time()), "\n")
cat("Log file path:            ", logFilePath, "\n\n")
cat("Phenotype file:           ", opt$phenoFile, "\n")
cat("Separator type: ", ifelse(is.null(opt$sepType), "default (',')", opt$sepType), "\n")
cat("IDAT folder:              ", opt$idatFolder, "\n")
cat("nSamples limit:           ", ifelse(is.na(opt$nSamples), "all", opt$nSamples), "\n")
cat("SampleID column:          ", opt$SampleID, "\n")
cat("Array type:               ", opt$arrayType, "\n")
cat("Annotation version:       ", opt$annotationVersion, "\n")
cat("Base RData folder:        ", opt$baseDataFolder, "\n")
cat("TIFF size (w x h @ dpi):  ", opt$tiffWidth, " x ", opt$tiffHeight, " @ ", opt$tiffRes, "\n")
cat("QC cutoff (median):       ", opt$qcCutoff, "\n")
cat("Detection P-value type:   ", opt$detPtype, "\n\n")
cat("Detection p-value threshold:", opt$detPThreshold, "\n")
cat("Normalization methods:    ", paste(opt$normMethodList, collapse = ", "), "\n")
cat("Funnorm seed:             ", opt$funnormSeed, "\n")
cat("Sex column:               ", opt$sexColumn, "\n")
cat("Plot grouping variable:   ", opt$plotGroupVar, "\n\n")
cat("Probe filtering:\n")
cat("  P-value threshold:      ", opt$pvalThreshold, "\n")
cat("  Chromosomes to remove:  ", opt$chrToRemove, "\n")
cat("  SNP positions filter:   ", opt$snpsToRemove, "\n")
cat("  MAF threshold:          ", opt$mafThreshold, "\n")
cat("  Cross-reactive file:    ", ifelse(is.null(opt$crossReactivePath), 
                                         "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv", 
                                         opt$crossReactivePath), "\n\n")
cat("Cell composition (estimateLC):\n")
cat("  Reference:              ", opt$lcRef, "\n")
cat("  Leading pheno order:    ", opt$phenoOrder, "\n")
# =============================================================================

# ----------- Prepare Subfolders for metrics -----------
objectDir  <- file.path(opt$baseDataFolder, opt$scriptLabel, "objects")
normDir    <- file.path(opt$baseDataFolder, opt$scriptLabel, "normObjects")
metricsDir <- file.path(opt$baseDataFolder, opt$scriptLabel, "metrics")
filterDir  <- file.path(opt$baseDataFolder, opt$scriptLabel, "filterObjects")

dir.create(objectDir, recursive = TRUE, showWarnings = FALSE)
dir.create(normDir, recursive = TRUE, showWarnings = FALSE)
dir.create(metricsDir, recursive = TRUE, showWarnings = FALSE)
dir.create(filterDir, recursive = TRUE, showWarnings = FALSE)

# ----------- Prepare Subfolders for rData/qc -----------

qcDir <- file.path("rData", opt$scriptLabel, "qc")

if (!dir.exists(qcDir)) {
  dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)
}

# ----------- Prepare Subfolders for figures/metrics -----------

metricsFigDir <- file.path("figures", opt$scriptLabel, "metrics")

if (!dir.exists(metricsFigDir)) {
  dir.create(metricsFigDir, recursive = TRUE, showWarnings = FALSE)
}
# ----------- Prepare Subfolders for figures/qc -----------

qcFigDir <- file.path("figures", opt$scriptLabel, "qc")

if (!dir.exists(qcFigDir)) {
  dir.create(qcFigDir, recursive = TRUE, showWarnings = FALSE)
}
# ----------- Prepare Subfolders for figures/enmix -----------
# Target folder
enmixDir <- file.path("figures", opt$scriptLabel, "enMix")

# Create directory if missing
if (!dir.exists(enmixDir)) {
  dir.create(enmixDir, recursive = TRUE, showWarnings = FALSE)
}
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
print(head(targets[, 1:5])) 
cat("=======================================================================\n")

# ----------- Load IDAT Files into RGSet -----------
RGSet <- read.metharray.exp(
        base = opt$idatFolder,
        targets = targets,
        extended = FALSE,
        recursive = FALSE,
        verbose = FALSE
)

# Assign custom sample names
sampleNames(RGSet) <- targets[[opt$SampleID]]
cat("RGSet loaded with", length(sampleNames(RGSet)), "samples.\n")
cat("=======================================================================\n")

owd <- getwd(); on.exit(setwd(owd), add = TRUE)
setwd(enmixDir)

op <- options(bitmapType = "cairo")
on.exit(options(op), add = TRUE)

# Generate ENmix control plots (JPGs will be created in enmixDir)
plotCtrl(RGSet)

setwd(owd)

cat("Generated ENmix control JPGs in:", enmixDir, "\n")

cat("=======================================================================\n")

# ----------- Apply Annotation -----------
annotation(RGSet) <- c(
        array = opt$arrayType,
        annotation = opt$annotationVersion
)
cat("Applied annotation: ", paste(annotation(RGSet), collapse = ", "), "\n")
cat("Manifest used:\n")
show(getManifest(RGSet))
cat("=======================================================================\n")

# ----------- Save RGSet -----------
RGSetPath <- file.path(objectDir, "RGSet.RData")
save(RGSet, file = RGSetPath)
cat("RGSet saved to: ", RGSetPath, "\n")
cat("=======================================================================\n")

# ----------- Preprocess Raw (create MSet) -----------
cat("Running preprocessRaw() to generate MSet...\n")
MSet <- preprocessRaw(RGSet)
cat("MSet created with", ncol(MSet), "samples and", nrow(MSet), "probes.\n")
cat("=======================================================================\n")

# Save MSet object
MSetPath <- file.path(objectDir, "MSet.RData")
save(MSet, file = MSetPath)
cat("MSet saved to:", MSetPath, "\n")
cat("=======================================================================\n")

# Display methylated and unmethylated intensity
cat("Preview of methylated intensities:\n")
print(head(getMeth(MSet)[, 1:3]))
cat("=======================================================================\n")
cat("Preview of unmethylated intensities:\n")
print(head(getUnmeth(MSet)[, 1:3]))
cat("=======================================================================\n")

# ----------- Ratio Conversion and Genome Mapping -----------
cat("Converting MSet to RatioSet and GSet...\n")
RatioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
cat("RatioSet created.\n")
print(RatioSet)
cat("=======================================================================\n")
GSet <- mapToGenome(RatioSet)
cat("GSet created.\n")
print(GSet)
cat("=======================================================================\n")

# Save RatioSet and GSet
RatioSetPath <- file.path(objectDir, "RatioSet.RData")
GSetPath <- file.path(objectDir, "GSet.RData")
save(RatioSet, file = RatioSetPath)
save(GSet, file = GSetPath)
cat("=======================================================================\n")

# ----------- Extract Methylation Metrics -----------
cat("Extracting methylation raw level metrics from GSet, these metrics are not saved...
    \n")

beta <- getBeta(GSet)
cat("Preview of beta values:\n")
print(head(beta[, 1:5]))
cat("=======================================================================\n")

m <- getM(GSet)
cat("Preview of M-values:\n")
print(head(m[, 1:5]))
cat("=======================================================================\n")

cn <- getCN(GSet)
cat("Preview of copy number values:\n")
print(head(cn[, 1:5]))

cat("=======================================================================\n")

# ----------- Quality Control Plot (from MSet) -----------
cat("Running QC plotting from MSet object...\n")
qc <- getQC(MSet)

qcPath <- file.path("figures", opt$scriptLabel, "qc", "quality_control(MSet).tiff")
tiff(file = qcPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")

plotQC(qc, badSampleCutoff = opt$qcCutoff)
dev.off()

cat("QC plot saved to: ", qcPath, "\n")
cat("=======================================================================\n")

# ----------- Calculate Detection P-values -----------
cat("Calculating detection p-values...\n")

detP <- minfi::detectionP(RGSet, type = opt$detPtype)
cat("Detection p-values calculated using type: ", opt$detPtype, "\n")

cat("Preview of detection p-values:\n")
print(head(detP[, 1:5]))

detPpath <- file.path(qcDir, "detP_RGSet.RData")
save(detP, file = detPpath)
cat("Detection RData p-values saved to: ", detPpath, "\n")

detPlotPath <- file.path("figures", opt$scriptLabel, "qc", "detection_pvalues(RGSet).tiff")

tiff(file = detPlotPath, 
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")
barplot(colMeans(detP), 
        las=3, 
        cex.names=0.8, 
        ylab="Mean detection p-values")
abline(h=0.05,col="red", lwd = 2, lty = 2) 
dev.off()

cat("Detection plot p-values saved to: ", detPlotPath, "\n")
cat("=======================================================================\n")

# ----------- Remove samples based on detection P-values -----------
cat("Calculate the mean detection p-values across all samples...\n")
meanDetP <- colMeans(detP)

# === Identify Failed Samples ===
failedSamples <- names(meanDetP[meanDetP > opt$detPThreshold])
nFailed <- length(failedSamples)
nBefore <- ncol(RGSet)

cat("Number of failed samples:", nFailed, "\n")
if (nFailed > 0) {
  cat("Failed sample IDs:\n")
  cat(paste(failedSamples, collapse = ", "), "\n")
}

# === Remove Failed Samples from RGSet ===
RGSet <- RGSet[, !(colnames(RGSet) %in% failedSamples)]
nAfter <- ncol(RGSet)

cat("Samples before filtering:", nBefore, "\n")
cat("Samples after filtering:", nAfter, "\n")

# ----------- Save RGSet -----------
RGSetPath <- file.path(objectDir, "RGSet.RData")
save(RGSet, file = RGSetPath)
cat("RGSet saved after removing the failed samples to: ", RGSetPath, "\n")
cat("=======================================================================\n")

# ----------- Density Plot of Beta Values from MSet -----------
cat("Generating density plot of Beta values...\n")

phenoData <- pData(MSet)

# Ensure output directory exists
denBetaPath <- file.path("figures", 
                         opt$scriptLabel, "qc", "densityBeta(MSet).tiff")

tiff(filename = denBetaPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")

densityPlot(MSet,
            sampGroups = phenoData[[opt$plotGroupVar]],
            pal = brewer.pal(8, "Dark2"),
            main = paste("Density Plot of Beta Values by", opt$plotGroupVar),
            add = TRUE,
            legend = TRUE)

dev.off()

cat("Density plot saved to: ", denBetaPath, "\n")
cat("=======================================================================\n")

cat("Predicting sex based on Beta values...\n")
pSex <- getSex(GSet) 
head(pSex)

# -------------- Plot Sex predictions --------------
pSexPath <- file.path("figures", 
                         opt$scriptLabel, "qc", "sexPrediction(GSet).tiff")

tiff(filename = pSexPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = 70, type = "cairo")

plot(x = pSex$xMed, 
     y = pSex$yMed,
     type = "n", 
     xlab = "X chr, median total intensity (log2)", 
     ylab = "Y chr, median total intensity (log2)")
text(x = pSex$xMed, y = pSex$yMed, labels = targets[[opt$SampleID]], 
     col = ifelse(pSex$predictedSex == "M", "deepskyblue", "deeppink3"))
legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
dev.off()

cat("Predicted Sex plot saved to: ", pSexPath, "\n")
cat("=======================================================================\n")

# Create clinical sex plot
pSexD <- as.data.frame(pSex) 
pSexD <- merge(pSexD, targets, by.x="row.names", by.y = opt$SampleID)
head(pSexD[, 1:4])

if (is.character(targets$Sex) || is.factor(targets$Sex)) {
  targets$Sex <- ifelse(targets$Sex %in% c("F", "Female", "f", "female", "FEMALE"), 0, 1)
}

# -------------- Plot Sex predictions --------------
pSexClPath <- file.path("figures", 
                      opt$scriptLabel, "qc", "sexClinical(GSet).tiff")

tiff(filename = pSexClPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = 70, type = "cairo")

plot(x = pSexD$xMed, y = pSexD$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
text(x = pSexD$xMed, y = pSexD$yMed, labels = pSexD$Row.names, 
     col = ifelse(pSexD$Sex == "1", "deepskyblue", "deeppink3"))
legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
dev.off()

cat("Clinical Sex plot saved to: ", pSexClPath, "\n")
cat("=======================================================================\n")

# Bind the predicted sex to the targets file and identify any mismatches 
targets$PredSex <- pSex$predictedSex
# Convert F = 0 and M = 1 in the column predSex
targets$PredSex <- ifelse(targets$PredSex == "F", 0, 1)

# === Remove Failed Samples from targets ===
targets <- targets[!(targets[[opt$SampleID]] %in% failedSamples), ]

# Add PredSex to pData  
pData(RGSet)$PredSex <- targets$PredSex

cat("Mistmaches found")
print(targets[targets$Sex != targets$PredSex, 1:3])
cat("=======================================================================\n")

cat("Running normalization methods using Minfi and WateRmelon: ", 
    paste(opt$normMethodList, collapse = ", "), "\n")

sexVec <- NULL
if (!is.null(opt$sexColumn) && opt$sexColumn %in% colnames(pData(RGSet))) {
  sexVec <- pData(RGSet)[, opt$sexColumn]
} else {
  cat("Note: sexColumn not found in pData(RGSet). 
      Fallback to NULL; funnorm/adjustedfunnorm will run without sex covariate.\n")
}

normPaths <- c(); firstMethod <- TRUE  
for (method in opt$normMethodList) {
        cat("  → Applying normalization:", method, "\n")
        set.seed(opt$funnormSeed)
        
        normObj <- switch(
                method,
                "adjustedfunnorm" = adjustedFunnorm(RGSet, sex = sexVec),
                "funnorm"         = preprocessFunnorm(RGSet, sex = sexVec),
                "illumina"        = preprocessIllumina(RGSet),
                "quantile"        = preprocessQuantile(RGSet, sex = sexVec),
                "swan"            = preprocessSWAN(RGSet),
                stop(paste("Unknown normalization method:", method))
        )
        if (method %in% c("funnorm","adjustedfunnorm") && is.null(sexVec)) {
          cat("Requested method uses sex, but sex not provided; 
              proceeded with sex = NULL.\n")
        }
        
        if (firstMethod) {
                MSetF <- normObj
                firstMethod <- FALSE
        }
        
        normPath <- file.path(normDir, paste0("norm_", method, "_RGSet.RData"))
        save(normObj, file = normPath)
        normPaths <- c(normPaths, normPath)
        cat("Saved normalized object: ", normPath, "\n")
}

# -------------- Plot Row vs Normalise data --------------
rawNormlPath <- file.path("figures", 
                        opt$scriptLabel, "qc", "sexComparison_RawNorm(MSetF).tiff")

tiff(filename = rawNormlPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")
     
par(mfrow=c(1,2))
densityPlot(RGSet, 
            sampGroups=targets[[opt$sexColumn]],
            main="Raw", 
            legend=FALSE)
legend("top", 
       legend = levels(factor(targets[[opt$sexColumn]])), 
       text.col=brewer.pal(8,"Dark2"))

densityPlot(getBeta(MSetF), 
            sampGroups=targets[[opt$sexColumn]],
            main="Normalized", 
            legend=FALSE)
legend("top", 
       legend = levels(factor(targets[[opt$sexColumn]])), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

cat("Plot Raw vs Normalisation data saved to: ", rawNormlPath, "\n")
cat("=======================================================================\n")

# ----------- Probe Filtering Based on Detection P-values -----------
cat("Filtering probes with detection p-values: ", 
    opt$pvalThreshold, "...\n")

# Recompute detection p-values
detP <- minfi::detectionP(RGSet)

# Align detection p-values with normalized probes
detP <- detP[match(featureNames(MSetF), rownames(detP)), ]

# Identify probes retained across all samples
keep <- rowSums(detP < opt$pvalThreshold) == ncol(MSetF)
cat("Probes retained: ", sum(keep), "/", length(keep), "\n")

MSetF_Flt <- MSetF[keep, ]
MSetFfltPath <- file.path(filterDir, "removProbes_MSetF_Flt.RData")
save(MSetF_Flt, file = MSetFfltPath)
cat("Filtered object saved to: ", MSetFfltPath, "\n")
cat("=======================================================================\n")

# ----------- Filter Probes on Sex Chromosomes -----------
cat("Removing probes on chromosomes: ", paste(opt$chrToRemoveList, 
                                              collapse = ", "), "\n")
# Identify probes to remove
ann <- getAnnotation(RGSet)
removeProbes <- ann$Name[ann$chr %in% opt$chrToRemoveList]
keepChr <- !(featureNames(MSetF_Flt) %in% removeProbes)

MSetF_Flt_Rxy <- MSetF_Flt[keepChr, ]

cat("Remaining probes after removing selected chromosomes:\n")
print(table(getAnnotation(MSetF_Flt_Rxy)$chr))

rxyPath <- file.path(filterDir, "removChrXY_MSetF_Flt_Rxy.RData")
save(MSetF_Flt_Rxy, file = rxyPath)
cat("Sex chromosome-filtered object saved to: ", rxyPath, "\n")
cat("=======================================================================\n")

# ----------- Remove Probes with SNPs -----------
cat("Removing probes with SNPs at: ", paste(opt$snpList, collapse = ", "), 
    " with MAF >=", opt$mafThreshold, "\n")

# Apply SNP probe filtering
MSetF_Flt_Rxy_Ds <- dropLociWithSnps(
        MSetF_Flt_Rxy,
        snps = opt$snpList,
        maf = opt$mafThreshold
)
cat("Remaining probes after SNP filtering: ", nrow(MSetF_Flt_Rxy_Ds), "\n")

snpPath <- file.path(filterDir, paste0("removSNPs_MAF", opt$mafThreshold, 
                                       "_MSetF_Flt_Rxy_Ds.RData"))
save(MSetF_Flt_Rxy_Ds, file = snpPath)
cat("SNP-filtered object saved to: ", snpPath, "\n")
cat("=======================================================================\n")

# ----------- Remove Cross-Reactive Probes -----------
cat("Loading cross-reactive probe list from:\n", opt$crossReactivePath, "\n")

xReactiveProbes <- read.csv(opt$crossReactivePath, stringsAsFactors = FALSE)

# Filter out cross-reactive probes
keepCr <- !(featureNames(MSetF_Flt_Rxy_Ds) %in% xReactiveProbes$ProbeID)
cat("Probes retained after cross-reactive filter: ", sum(keepCr), "\n")

MSetF_Flt_Rxy_Ds_Rc <- MSetF_Flt_Rxy_Ds[keepCr, ]
rcPath <- file.path(filterDir, "removCrossReactive_MSetF_Flt_Rxy_Ds_Rc.RData")
save(MSetF_Flt_Rxy_Ds_Rc, file = rcPath)
cat("Cross-reactive-filtered object saved to: ", rcPath, "\n")
cat("=======================================================================\n")

# ----------- Final DNAm Matrices from Filtered Data -----------
cat("Extracting final DNAm matrices (M, Beta, CN)...\n")

# M-values
m <- getM(MSetF_Flt_Rxy_Ds_Rc)
mOut <- file.path(metricsDir, "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(m, file = mOut)
cat("M-values saved to: ", mOut, "\n")
print(head(m[, 1:5]))

# Beta-values
beta <- getBeta(MSetF_Flt_Rxy_Ds_Rc)
betaOut <- file.path(metricsDir, "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(beta, file = betaOut)
cat("Beta-values saved to: ", betaOut, "\n")
print(head(beta[, 1:5]))

# CN-values
cn <- getCN(MSetF_Flt_Rxy_Ds_Rc)
cnOut <- file.path(metricsDir, "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(cn, file = cnOut)
cat("CN-values saved to: ", cnOut, "\n")
print(head(cn[, 1:5]))
cat("=======================================================================\n")

# ----- Examine higher dimensions to look at other sources of variation -----

groupFactor <- factor(targets[[opt$plotGroupVar]])
groupSex <- factor(targets[[opt$sexColumn]])

mdsPath <- file.path("figures", 
                          opt$scriptLabel, 
                          "metrics", 
                          "examineMDS_PostFilteringCrossRect(MSetF_Flt_Rxy_Ds_Rc).tiff")

tiff(filename = mdsPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")

pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
plotMDS(getM(MSetF_Flt_Rxy_Ds_Rc),
        main="Timepoint",
        top=1000, gene.selection="common", 
        col=pal[groupFactor], dim=c(1,2))
legend("right", legend=levels(groupFactor), 
       text.col = brewer.pal(8,"Dark2"),
       cex=0.7, bg="white")
plotMDS(getM(MSetF_Flt_Rxy_Ds_Rc), 
        main="Sex",
        top=1000, gene.selection="common", 
        col=pal[groupSex], dim=c(2,3))
legend("topright", legend=levels(groupSex), 
       text.col = brewer.pal(8,"Dark2"),
       cex=0.7, bg="white")

dev.off()

cat("Plot examineMDS_PostFilteringCrossRect saved to: ", mdsPath, "\n")
cat("=======================================================================\n")

# ----------- Plot Density of Final Beta & M Values by Group Variable -----------
cat("Plotting final density plots for grouping variable: ", 
    opt$plotGroupVar, "\n")

betaMPlotPath <- file.path("figures", 
                     opt$scriptLabel, 
                     "metrics", 
                     "densityBeta&M(MSetF_Flt_Rxy_Ds_Rc).tiff")

# Create TIFF output
tiff(betaMPlotPath, 
     width = opt$tiffWidth, 
     height = opt$tiffHeight, 
     res = opt$tiffRes, type = "cairo")
par(mfrow = c(1, 2))

# Beta plot
densityPlot(beta,
            sampGroups = groupFactor,
            main = "Beta values",
            legend = FALSE,
            xlab = "Beta values")
legend("top", legend = levels(groupFactor), text.col = brewer.pal(8,"Dark2"))

# M-value plot
densityPlot(m,
            sampGroups = groupFactor,
            main = "M-values",
            legend = FALSE,
            xlab = "M values")
legend("topleft", legend = levels(groupFactor), text.col = brewer.pal(8,"Dark2"))

dev.off()
cat("Density plots saved to: ", betaMPlotPath, "\n")

cat("=======================================================================\n")

# ----------- Cell Type Estimation using Ewastool -----------

cat("Estimating cell composition with ref:", opt$lcRef, "\n")

lc <- ewastools::estimateLC(beta, 
                            ref = opt$lcRef, constrained = TRUE)
phenoLC <- cbind(targets, lc)

leadCols <- strsplit(opt$phenoOrder, ";", fixed = TRUE)[[1]]
leadCols <- leadCols[leadCols %in% colnames(phenoLC)]
phenoLC <- dplyr::select(phenoLC, dplyr::all_of(leadCols), dplyr::everything())

if (!dir.exists(opt$lcPhenoDir)) dir.create(opt$lcPhenoDir, recursive = TRUE)
lcPhenoOut <- file.path(opt$lcPhenoDir, "phenoLC.csv")
write.csv(phenoLC, 
          file = lcPhenoOut, 
          row.names = FALSE)
cat("Saved phenoLC:", lcPhenoOut, "\n")

cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")  
sink()                  
close(logCon)           


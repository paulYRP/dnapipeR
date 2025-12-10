#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(rmarkdown)
  library(optparse)
}))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--rmd", type="character"),
  make_option("--output", type="character"),
  make_option("--outputDir", type="character"),
  make_option("--qcDir", type="character"),
  make_option("--preDir", type="character"),
  make_option("--postDir", type="character"),
  make_option("--svaDir", type="character"),
  make_option("--glmDir", type="character"),
  make_option("--glmmDir", type="character"),
  make_option("--title", type="character"),
  make_option("--author", type="character"),
  make_option("--date", type="character"),
  make_option("--figDir", type="character")
)))

# Fix Windows path handling for Pandoc/LaTeX:
normalize_path <- function(x) gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))

opt$rmd        <- normalize_path(opt$rmd)
opt$outputDir  <- normalize_path(opt$outputDir)
opt$qcDir      <- normalize_path(opt$qcDir)
opt$preDir     <- normalize_path(opt$preDir)
opt$postDir    <- normalize_path(opt$postDir)
opt$svaDir     <- normalize_path(opt$svaDir)
opt$glmDir     <- normalize_path(opt$glmDir)
opt$glmmDir    <- normalize_path(opt$glmmDir)
opt$figDir    <- normalize_path(opt$figDir)

render(
  input = opt$rmd,
  output_file = opt$output,
  output_dir = opt$outputDir,
  params = list(
    reportTitle = opt$title,
    author = opt$author,
    date = opt$date,
    qcDir = opt$qcDir,
    preprocessingDir = opt$preDir,
    postprocessingDir = opt$postDir,
    svaDir = opt$svaDir,
    glmDir = opt$glmDir,
    glmmDir = opt$glmmDir,
    figDir = opt$figDir
  ),
  knit_root_dir = dirname(opt$rmd)
)

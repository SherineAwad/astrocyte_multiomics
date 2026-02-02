#!/usr/bin/env Rscript
library(ArchR)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript script.R <config.csv> <project_name> <genome> <threads>")
}

config_file <- args[1]
project_name <- args[2]
genome <- args[3]
threads <- as.integer(args[4])
output_dir <- project_name

addArchRGenome(genome)
addArchRThreads(threads = threads)

config <- read.csv(config_file)
atacFiles <- setNames(config$atac_file, config$sample_name)

for (i in 1:nrow(config)) {
  rna_file <- config$rna_file[i]
  if (!file.exists(rna_file)) {
    stop(paste("RNA file does not exist:", rna_file))
  }
}

ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

rna_data <- list()
for(i in 1:nrow(config)) {
  sample <- config$sample_name[i]
  proj <- ArchRProject(paste0(sample, ".arrow"), outputDirectory = sample, copyArrows = FALSE)
  
  seRNA <- tryCatch({
    import10xFeatureMatrix(input = config$rna_file[i], names = sample)
  }, error = function(e) {
    stop(paste("Failed to import RNA file:", config$rna_file[i], "\nError:", e$message))
  })
  
  k <- seqnames(rowRanges(seRNA)) %in% seqnames(proj@genomeAnnotation$chromSizes)
  seRNA <- seRNA[k,]
  seRNA <- seRNA[, colnames(seRNA) %in% rownames(proj@cellColData)]
  rna_data[[sample]] <- seRNA
}

all_counts <- do.call(cbind, lapply(rna_data, assay))
seRNA_all <- SummarizedExperiment(
  assays = list(counts = all_counts),
  rowRanges = rowRanges(rna_data[[1]])
)

proj_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = output_dir, copyArrows = TRUE)
proj_ALL <- addTileMatrix(proj_ALL, binarize = FALSE, force = TRUE)
proj_ALL <- addGeneExpressionMatrix(input = proj_ALL, seRNA = seRNA_all)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveArchRProject(
  ArchRProj = proj_ALL,
  outputDirectory = output_dir,
  load = FALSE
)

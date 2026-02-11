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

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 1,  # Lower for now
  minFrags = 100,  # Lower for now
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Create project FIRST
proj_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = output_dir, copyArrows = TRUE)

# Add RNA data ALL AT ONCE with proper alignment
seRNA_list <- list()
for(i in 1:nrow(config)) {
  sample <- config$sample_name[i]
  seRNA <- import10xFeatureMatrix(input = config$rna_file[i], names = sample)
  seRNA_list[[sample]] <- seRNA
}

# Combine RNA data
seRNA_combined <- do.call(cbind, seRNA_list)

# Add to project with forced cell matching
proj_ALL <- addGeneExpressionMatrix(
  input = proj_ALL, 
  seRNA = seRNA_combined,
  force = TRUE
)

# Save
saveArchRProject(
  ArchRProj = proj_ALL,
  outputDirectory = output_dir,
  load = FALSE
)

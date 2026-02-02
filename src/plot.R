#!/usr/bin/env Rscript

library(ArchR)
library(ggplot2)
library(argparse)

# Parser setup
parser <- ArgumentParser(description = 'Create UMAP plots from existing ArchR project')
parser$add_argument('--project_name', type = 'character', required = TRUE)
parser$add_argument('--suffix', type = 'character', default = '')
args <- parser$parse_args()

# Suppress Rplots.pdf creation
pdf(NULL)

# Load existing project from current directory
proj_ALL <- loadArchRProject(path = args$project_name, force = FALSE, showLogo = TRUE)

cat(sprintf("Loaded project with %d cells\n", nCells(proj_ALL)))

# Create UMAP plots
p1 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "Clusters_RNA",
  embedding = "UMAP_RNA",
  plotAs = "points",
  size = 0.5
)

p2 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "Clusters_ATAC",
  embedding = "UMAP_ATAC",
  plotAs = "points",
  size = 0.5
)

p3 <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "Clusters_Combined",
  embedding = "UMAP_Combined",
  plotAs = "points",
  size = 0.5
)

# Save as PDF instead of PNG
base_name <- basename(args$project_name)
pdf(paste0(base_name, args$suffix, "_RNA_UMAP.pdf"), width = 12, height = 8)
print(p1)
dev.off()

pdf(paste0(base_name, args$suffix, "_ATAC_UMAP.pdf"), width = 12, height = 8)
print(p2)
dev.off()

pdf(paste0(base_name, args$suffix, "_Combined_UMAP.pdf"), width = 12, height = 8)
print(p3)
dev.off()

cat("=== PDF plots created successfully ===\n")

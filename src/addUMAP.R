#!/usr/bin/env Rscript

library(Matrix)   # Must load first, or ArchR will error later
library(irlba)    # Load irlba first so we can patch it

# FIX: Replace irlba function BEFORE loading ArchR
assignInNamespace("irlba",
  function(A, nv, nu = nv, ...) {
    cat("Using svd() instead of irlba to avoid Matrix conflict\n")
    # svd needs different parameter names
    result <- svd(A, nu = min(nu, min(dim(A))), nv = min(nv, min(dim(A))))
    # Format result to match irlba output
    result$d <- result$d[1:min(nv, length(result$d))]
    result$u <- result$u[, 1:min(nu, ncol(result$u)), drop = FALSE]
    result$v <- result$v[, 1:min(nv, ncol(result$v)), drop = FALSE]
    return(result)
  },
  ns = "irlba"
)

library(ArchR)

# Load other packages
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)
library(pheatmap)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(argparse)

# Parser setup (same as before)
parser <- ArgumentParser(description = 'ArchR single-cell ATAC-seq processing with LSI and clustering')

parser$add_argument('--project_name', type = 'character', required = TRUE)
parser$add_argument('--genome', type = 'character', default = 'mm10')
parser$add_argument('--threads', type = 'integer', default = 4)
parser$add_argument('--suffix', type = 'character', default = '_processed')

parser$add_argument('--lsi_resolution', type = 'double', default = 0.2)
parser$add_argument('--lsi_sample_cells', type = 'integer', default = 10000)
parser$add_argument('--lsi_nstart', type = 'integer', default = 10)
parser$add_argument('--rna_var_features', type = 'integer', default = 2500)

parser$add_argument('--atac_cluster_res', type = 'double', default = 1.0)
parser$add_argument('--rna_cluster_res', type = 'double', default = 1.0)
parser$add_argument('--combined_cluster_res', type = 'double', default = 3.0)
parser$add_argument('--combined_dims', type = 'integer', default = 55)
parser$add_argument('--cluster_dims', type = 'integer', default = 20)

parser$add_argument('--umap_mindist', type = 'double', default = 0.8)

parser$add_argument('--plot_width', type = 'integer', default = 1200)
parser$add_argument('--plot_height', type = 'integer', default = 800)
parser$add_argument('--plot_res', type = 'integer', default = 150)

args <- parser$parse_args()

addArchRThreads(threads = args$threads)
addArchRGenome(args$genome)

input_dir  <- args$project_name
output_dir <- paste0(args$project_name, args$suffix)

proj_ALL <- loadArchRProject(path = input_dir, force = FALSE, showLogo = TRUE)

cat(sprintf("Loaded project with %d cells\n", nCells(proj_ALL)))

proj_ALL <- addIterativeLSI(
  ArchRProj = proj_ALL,
  clusterParams = list(
    resolution = args$lsi_resolution,
    sampleCellsPre = args$lsi_sample_cells,
    n.start = args$lsi_nstart
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

proj_ALL <- addIterativeLSI(
  ArchRProj = proj_ALL,
  clusterParams = list(
    resolution = args$lsi_resolution,
    sampleCells = args$lsi_sample_cells,
    n.start = args$lsi_nstart
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  varFeatures = args$rna_var_features,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

proj_ALL <- addCombinedDims(proj_ALL,
                            reducedDims = c("LSI_ATAC", "LSI_RNA"),
                            name = "LSI_Combined")

proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_ATAC",
                    name = "UMAP_ATAC", minDist = args$umap_mindist, force = TRUE)

proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_RNA",
                    name = "UMAP_RNA", minDist = args$umap_mindist, force = TRUE)

proj_ALL <- addUMAP(proj_ALL, reducedDims = "LSI_Combined",
                    dimsToUse = 1:args$combined_dims, name = "UMAP_Combined",
                    minDist = args$umap_mindist, force = TRUE)

proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_ATAC",
                        name = "Clusters_ATAC",
                        resolution = args$atac_cluster_res, force = TRUE)

proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_RNA",
                        name = "Clusters_RNA",
                        resolution = args$rna_cluster_res, force = TRUE)

proj_ALL <- addClusters(proj_ALL, reducedDims = "LSI_Combined",
                        dimsToUse = 1:args$cluster_dims,
                        name = "Clusters_Combined",
                        resolution = args$combined_cluster_res, force = TRUE)

saveArchRProject(proj_ALL, outputDirectory = output_dir, load = FALSE)

cat("=== ArchR pipeline completed successfully ===\n")

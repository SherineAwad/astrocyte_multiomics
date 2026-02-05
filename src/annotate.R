#!/usr/bin/env Rscript

library(ArchR)
library(ggplot2)
library(argparse)

# Parser setup
parser <- ArgumentParser(description = 'Add cell type annotations to ArchR project')
parser$add_argument('--project_name', type = 'character', required = TRUE, help = 'Path to ArchR project')
parser$add_argument('--annotation_file', type = 'character', required = TRUE, help = 'Path to annotation CSV file WITHOUT header')
parser$add_argument('--annotation_col', type = 'character', default = 'CellType', help = 'Name to assign to cell type column internally')
parser$add_argument('--cluster_col', type = 'character', default = 'Cluster', help = 'Name to assign to cluster ID column internally')
parser$add_argument('--suffix', type = 'character', default = '_annotated', help = 'Suffix for output project directory')
args <- parser$parse_args()

# Load ArchR project
cat("Loading ArchR project...\n")
proj_ALL <- loadArchRProject(path = args$project_name, force = FALSE, showLogo = TRUE)
cat(sprintf("Loaded project with %d cells\n", nCells(proj_ALL)))

# Load annotation file WITHOUT header, expecting exactly 2 columns: cluster, celltype
cat("Reading annotation file WITHOUT header, expecting exactly 2 columns...\n")
annotation_df <- read.csv(args$annotation_file, header = FALSE, stringsAsFactors = FALSE)
colnames(annotation_df) <- c(args$cluster_col, args$annotation_col)
cat(sprintf("Loaded annotation file with %d rows\n", nrow(annotation_df)))

# Convert cluster IDs to character for matching
annotation_df[[args$cluster_col]] <- as.character(annotation_df[[args$cluster_col]])
proj_clusters <- as.character(proj_ALL$Clusters_Combined)

# *** STRIP leading 'C' from project cluster IDs to match annotation ***
proj_clusters <- gsub("^C", "", proj_clusters)

cat("\nUnique clusters in ArchR project after stripping 'C':\n")
print(sort(unique(proj_clusters)))

cat("\nClusters from annotation file:\n")
print(sort(unique(annotation_df[[args$cluster_col]])))

# Find clusters missing in annotation file
missing_clusters <- setdiff(unique(proj_clusters), annotation_df[[args$cluster_col]])
if(length(missing_clusters) > 0){
  cat("\nWarning: The following clusters are missing in the annotation file and will be labeled 'Unknown':\n")
  print(missing_clusters)
}

# Create named vector for cluster -> celltype
cluster_to_celltype <- setNames(annotation_df[[args$annotation_col]], annotation_df[[args$cluster_col]])

# Assign 'Unknown' for missing clusters
for(clust in missing_clusters){
  cluster_to_celltype[clust] <- "Unknown"
}

# Make sure cluster_to_celltype has entries for all clusters in proj_clusters
cluster_to_celltype <- cluster_to_celltype[unique(proj_clusters)]

# Assign CellType annotation to project metadata
proj_ALL$CellType <- cluster_to_celltype[proj_clusters]

cat("\nSummary of assigned CellType annotations:\n")
print(table(proj_ALL$CellType))

# Plot and save UMAP with cell type labels
figure_name <- paste0(basename(args$project_name), "_annotated.png")
p <- plotEmbedding(
  ArchRProj = proj_ALL,
  colorBy = "cellColData",
  name = "CellType",
  embedding = "UMAP_Combined",
  labelAsFactors = FALSE,
  labelMeans = TRUE
)
ggsave(filename = figure_name, plot = p, width = 6, height = 5, dpi = 300)
cat(sprintf("Saved UMAP plot to %s\n", figure_name))

# Save annotated project
output_dir <- paste0(args$project_name, args$suffix)
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = output_dir, load = FALSE)
cat(sprintf("Saved annotated project to %s\n", output_dir))


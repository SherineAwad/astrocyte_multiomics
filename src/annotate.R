#!/usr/bin/env Rscript

library(ArchR)
library(ggplot2)
library(argparse)

# ----------------------------
# Parser setup
# ----------------------------
parser <- ArgumentParser(description = 'Add cell type annotations to ArchR project')
parser$add_argument('--project_name', type = 'character', required = TRUE,
                    help = 'Path to ArchR project')
parser$add_argument('--annotation_file', type = 'character', required = TRUE,
                    help = 'Path to annotation CSV file WITHOUT header')
parser$add_argument('--annotation_col', type = 'character', default = 'CellType',
                    help = 'Name to assign to cell type column internally')
parser$add_argument('--cluster_col', type = 'character', default = 'Cluster',
                    help = 'Name to assign to cluster ID column internally')
parser$add_argument('--suffix', type = 'character', default = '_annotated',
                    help = 'Suffix for output project directory')
parser$add_argument('--remove_clusters', type = 'integer', nargs = '*', default = NULL,
                    help = 'Clusters to remove from project (e.g., 1 2 3)')

args <- parser$parse_args()

# ----------------------------
# Load ArchR project
# ----------------------------
cat("Loading ArchR project...\n")
proj_ALL <- loadArchRProject(path = args$project_name, force = FALSE, showLogo = TRUE)
cat(sprintf("Loaded project with %d cells\n", nCells(proj_ALL)))

# ----------------------------
# Optionally remove clusters
# ----------------------------
if(!is.null(args$remove_clusters)){
  cat("Removing specified clusters:\n")
  print(args$remove_clusters)
  
  proj_clusters <- as.character(proj_ALL$Clusters_Combined)
  proj_clusters_stripped <- gsub("^C", "", proj_clusters)
  
  keep_cells <- !proj_clusters_stripped %in% as.character(args$remove_clusters)
  cat(sprintf("Removing %d cells from clusters %s\n",
              sum(!keep_cells), paste(args$remove_clusters, collapse = ", ")))
  
  proj_ALL <- proj_ALL[keep_cells, ]
  cat(sprintf("Project now has %d cells after removal\n", nCells(proj_ALL)))
}

# ----------------------------
# Load annotation file
# ----------------------------
cat("Reading annotation file WITHOUT header, expecting exactly 2 columns...\n")
annotation_df <- read.csv(args$annotation_file, header = FALSE, stringsAsFactors = FALSE)
colnames(annotation_df) <- c(args$cluster_col, args$annotation_col)
cat(sprintf("Loaded annotation file with %d rows\n", nrow(annotation_df)))

# ----------------------------
# Assign cell type annotations
# ----------------------------
proj_clusters <- as.character(proj_ALL$Clusters_Combined)
proj_clusters <- gsub("^C", "", proj_clusters)

cat("\nUnique clusters in ArchR project after stripping 'C':\n")
print(sort(unique(proj_clusters)))

cat("\nClusters from annotation file:\n")
print(sort(unique(annotation_df[[args$cluster_col]])))

# Find missing clusters
missing_clusters <- setdiff(unique(proj_clusters), annotation_df[[args$cluster_col]])
if(length(missing_clusters) > 0){
  cat("\nWarning: The following clusters are missing in the annotation file and will be labeled 'Unknown':\n")
  print(missing_clusters)
}

# Map clusters -> cell types
cluster_to_celltype <- setNames(annotation_df[[args$annotation_col]], annotation_df[[args$cluster_col]])

# Assign 'Unknown' to missing clusters
for(clust in missing_clusters){
  cluster_to_celltype[clust] <- "Unknown"
}

# Ensure all clusters in proj_clusters have a mapping
cluster_to_celltype <- cluster_to_celltype[unique(proj_clusters)]
proj_ALL$CellType <- cluster_to_celltype[proj_clusters]

cat("\nSummary of assigned CellType annotations:\n")
print(table(proj_ALL$CellType))

# ----------------------------
# Plot UMAP
# ----------------------------
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

# ----------------------------
# Save annotated project
# ----------------------------
output_dir <- paste0(args$project_name, args$suffix)
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = output_dir, load = FALSE)
cat(sprintf("Saved annotated project to %s\n", output_dir))


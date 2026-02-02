#!/usr/bin/env Rscript

library(ArchR)
library(ggplot2)
library(argparse)
library(gridExtra)
library(patchwork)

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

# Helper function to check if embedding exists
embedding_exists <- function(proj, embedding_name) {
  tryCatch({
    emb <- getEmbedding(ArchRProj = proj, embedding = embedding_name)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# Get cell metadata
cell_col_data <- getCellColData(proj_ALL)

# Define embeddings
embeddings <- c("UMAP_RNA", "UMAP_ATAC", "UMAP_Combined")
embedding_names <- c("RNA", "ATAC", "Combined")

# ORIGINAL PLOTS - EXACT SAME AS BEFORE
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

# Save ORIGINAL plots with ORIGINAL names
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

cat("Created original 3 UMAP plots\n")

# 1. NEW: UMAP per sample (with NEW file names)
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating UMAP plots per sample...\n")
  
  for(i in seq_along(embeddings)) {
    if(embedding_exists(proj_ALL, embeddings[i])) {
      p_sample <- plotEmbedding(
        ArchRProj = proj_ALL,
        colorBy = "cellColData",
        name = "Sample",
        embedding = embeddings[i],
        plotAs = "points",
        size = 0.5
      )
      
      pdf(paste0(base_name, args$suffix, "_", embedding_names[i], "_SAMPLE_UMAP.pdf"), 
          width = 12, height = 8)
      print(p_sample)
      dev.off()
      
      cat(sprintf("  Created %s_SAMPLE_UMAP.pdf\n", embedding_names[i]))
    }
  }
}

# 2. NEW: UMAP per sample per modality (with NEW file names)
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating UMAP per sample per modality...\n")
  
  samples <- unique(cell_col_data$Sample)
  
  for(embed_idx in seq_along(embeddings)) {
    if(embedding_exists(proj_ALL, embeddings[embed_idx])) {
      sample_plots <- list()
      
      for(sample in samples) {
        sample_cells <- which(cell_col_data$Sample == sample)
        proj_sample <- proj_ALL[sample_cells, ]
        
        p <- plotEmbedding(
          ArchRProj = proj_sample,
          colorBy = "cellColData",
          name = "Clusters_Combined",
          embedding = embeddings[embed_idx],
          plotAs = "points",
          size = 1.0
        ) + ggtitle(paste(sample, "-", embedding_names[embed_idx]))
        
        sample_plots[[sample]] <- p
      }
      
      n_cols <- min(3, length(samples))
      n_rows <- ceiling(length(samples) / n_cols)
      
      combined_plot <- wrap_plots(sample_plots, ncol = n_cols, nrow = n_rows)
      
      pdf(paste0(base_name, args$suffix, "_", embedding_names[embed_idx], "_SAMPLES_GRID.pdf"), 
          width = 6 * n_cols, height = 5 * n_rows)
      print(combined_plot)
      dev.off()
      
      cat(sprintf("  Created %s_SAMPLES_GRID.pdf\n", embedding_names[embed_idx]))
    }
  }
}

# 3. NEW: QC plots per clusters (with NEW file names)
cat("Creating QC plots per clusters...\n")

cluster_columns <- c("Clusters_RNA", "Clusters_ATAC", "Clusters_Combined")
cluster_columns <- cluster_columns[cluster_columns %in% colnames(cell_col_data)]

qc_metrics <- c("TSSEnrichment", "nFrags", "BlacklistRatio", "NucleosomeRatio")
qc_metrics <- qc_metrics[qc_metrics %in% colnames(cell_col_data)]

if(length(cluster_columns) > 0 && length(qc_metrics) > 0) {
  for(cluster_col in cluster_columns) {
    qc_plot_list <- list()
    
    for(qc_metric in qc_metrics) {
      df <- data.frame(
        Cluster = cell_col_data[[cluster_col]],
        Value = cell_col_data[[qc_metric]]
      )
      
      df <- df[complete.cases(df), ]
      
      if(nrow(df) > 0) {
        p <- ggplot(df, aes(x = Cluster, y = Value, fill = Cluster)) +
          geom_violin(scale = "width", trim = TRUE) +
          geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
          theme_classic() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
          ) +
          labs(title = paste(qc_metric, "by", cluster_col))
        
        qc_plot_list[[qc_metric]] <- p
      }
    }
    
    if(length(qc_plot_list) > 0) {
      n_cols <- min(2, length(qc_plot_list))
      n_rows <- ceiling(length(qc_plot_list) / n_cols)
      
      combined_qc_plot <- wrap_plots(qc_plot_list, ncol = n_cols, nrow = n_rows)
      
      pdf(paste0(base_name, args$suffix, "_", cluster_col, "_QC.pdf"), 
          width = 8 * n_cols, height = 6 * n_rows)
      print(combined_qc_plot)
      dev.off()
      
      cat(sprintf("  Created %s_QC.pdf\n", cluster_col))
    }
  }
}

cat("=== ALL PLOTS CREATED SUCCESSFULLY ===\n")
cat("Original files (unchanged):\n")
cat("  ", base_name, args$suffix, "_RNA_UMAP.pdf\n", sep="")
cat("  ", base_name, args$suffix, "_ATAC_UMAP.pdf\n", sep="")
cat("  ", base_name, args$suffix, "_Combined_UMAP.pdf\n", sep="")
cat("\nNew files added:\n")
cat("  *_SAMPLE_UMAP.pdf (colored by sample)\n")
cat("  *_SAMPLES_GRID.pdf (each sample separately)\n")
cat("  *_QC.pdf (QC metrics per cluster)\n")

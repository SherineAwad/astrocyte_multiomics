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
    # Try to access the embedding
    emb <- getEmbedding(ArchRProj = proj, embedding = embedding_name)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# Get cell metadata for QC plots
cell_col_data <- getCellColData(proj_ALL)

# 1. UMAP per sample
# Check if 'Sample' column exists in cellColData
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating UMAP plots per sample...\n")
  
  # Get unique samples
  samples <- unique(cell_col_data$Sample)
  
  # Create sample-colored UMAPs for each embedding type
  embeddings <- c("UMAP_RNA", "UMAP_ATAC", "UMAP_Combined")
  embedding_names <- c("RNA", "ATAC", "Combined")
  
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
      
      base_name <- basename(args$project_name)
      pdf(paste0(base_name, args$suffix, "_", embedding_names[i], "_UMAP_per_sample.pdf"), 
          width = 12, height = 8)
      print(p_sample)
      dev.off()
      
      cat(sprintf("  Created %s UMAP colored by sample\n", embedding_names[i]))
    } else {
      cat(sprintf("  Skipping %s embedding - not found\n", embedding_names[i]))
    }
  }
} else {
  cat("Warning: 'Sample' column not found in cellColData. Skipping UMAP per sample plots.\n")
}

# 2. UMAP per sample per modality (ATAC/RNA/Combined)
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating UMAP plots per sample per modality...\n")
  
  samples <- unique(cell_col_data$Sample)
  
  # Loop through each embedding type
  for(embed_idx in seq_along(embeddings)) {
    if(embedding_exists(proj_ALL, embeddings[embed_idx])) {
      # Create a list to store plots for each sample
      sample_plots <- list()
      
      for(sample in samples) {
        # Get cell indices for this sample
        sample_cells <- which(cell_col_data$Sample == sample)
        proj_sample <- proj_ALL[sample_cells, ]
        
        # Create UMAP for this sample
        p <- plotEmbedding(
          ArchRProj = proj_sample,
          colorBy = "cellColData",
          name = "Clusters_Combined",  # Using combined clusters as default
          embedding = embeddings[embed_idx],
          plotAs = "points",
          size = 1.0,
          baseSize = 10
        ) + 
          ggtitle(paste(sample, "-", embedding_names[embed_idx])) +
          theme(plot.title = element_text(hjust = 0.5, size = 12))
        
        sample_plots[[sample]] <- p
      }
      
      # Arrange all sample plots in a grid
      n_cols <- min(3, length(samples))
      n_rows <- ceiling(length(samples) / n_cols)
      
      combined_plot <- wrap_plots(sample_plots, ncol = n_cols, nrow = n_rows) + 
        plot_annotation(title = paste(embedding_names[embed_idx], "UMAP - Per Sample"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
      
      # Save combined plot
      base_name <- basename(args$project_name)
      pdf(paste0(base_name, args$suffix, "_", embedding_names[embed_idx], 
                 "_UMAP_per_sample_grid.pdf"), 
          width = 6 * n_cols, height = 5 * n_rows)
      print(combined_plot)
      dev.off()
      
      cat(sprintf("  Created %s UMAP grid with %d samples\n", 
                  embedding_names[embed_idx], length(samples)))
    } else {
      cat(sprintf("  Skipping %s embedding grid - embedding not found\n", embedding_names[embed_idx]))
    }
  }
}

# Original UMAP plots (RNA, ATAC, Combined)
cat("Creating standard UMAP plots...\n")

# Check which embeddings exist and create plots
for(modality in c("RNA", "ATAC", "Combined")) {
  embedding_name <- paste0("UMAP_", modality)
  cluster_name <- paste0("Clusters_", modality)
  
  if(embedding_exists(proj_ALL, embedding_name) && 
     cluster_name %in% colnames(cell_col_data)) {
    
    p <- plotEmbedding(
      ArchRProj = proj_ALL,
      colorBy = "cellColData",
      name = cluster_name,
      embedding = embedding_name,
      plotAs = "points",
      size = 0.5
    )
    
    base_name <- basename(args$project_name)
    pdf(paste0(base_name, args$suffix, "_", modality, "_UMAP.pdf"), width = 12, height = 8)
    print(p)
    dev.off()
    
    cat(sprintf("  Created %s UMAP\n", modality))
  } else {
    if(!embedding_exists(proj_ALL, embedding_name)) {
      cat(sprintf("  Skipping %s UMAP - embedding not found\n", modality))
    } else if(!cluster_name %in% colnames(cell_col_data)) {
      cat(sprintf("  Skipping %s UMAP - %s not found in cellColData\n", modality, cluster_name))
    }
  }
}

# 3. QC plots per clusters
cat("Creating QC plots per clusters...\n")

# Define cluster columns to check
cluster_columns <- c("Clusters_RNA", "Clusters_ATAC", "Clusters_Combined")
cluster_columns <- cluster_columns[cluster_columns %in% colnames(cell_col_data)]

# Define QC metrics to plot
qc_metrics <- c("TSSEnrichment", "nFrags", "BlacklistRatio", "NucleosomeRatio")
# Check which QC metrics exist
qc_metrics <- qc_metrics[qc_metrics %in% colnames(cell_col_data)]

if(length(cluster_columns) > 0 && length(qc_metrics) > 0) {
  # Create QC plots for each cluster column
  for(cluster_col in cluster_columns) {
    cat(sprintf("  Creating QC plots for %s...\n", cluster_col))
    
    # Create a list to store QC plots
    qc_plot_list <- list()
    
    for(qc_metric in qc_metrics) {
      # Create data frame for plotting
      df <- data.frame(
        Cluster = cell_col_data[[cluster_col]],
        Value = cell_col_data[[qc_metric]]
      )
      
      # Remove NAs if any
      df <- df[complete.cases(df), ]
      
      # Only plot if we have data
      if(nrow(df) > 0) {
        p <- ggplot(df, aes(x = Cluster, y = Value, fill = Cluster)) +
          geom_violin(scale = "width", trim = TRUE) +
          geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
          theme_classic() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5)
          ) +
          labs(title = paste(qc_metric, "by", cluster_col),
               x = "Cluster",
               y = qc_metric)
        
        # Handle outliers for better visualization
        if(qc_metric %in% c("nFrags", "TSSEnrichment")) {
          p <- p + scale_y_log10() + annotation_logticks(sides = "l")
        }
        
        qc_plot_list[[qc_metric]] <- p
      } else {
        cat(sprintf("    Warning: No data for %s vs %s\n", qc_metric, cluster_col))
      }
    }
    
    # Only proceed if we have plots
    if(length(qc_plot_list) > 0) {
      # Arrange all QC plots in a grid
      n_cols <- min(2, length(qc_plot_list))
      n_rows <- ceiling(length(qc_plot_list) / n_cols)
      
      combined_qc_plot <- wrap_plots(qc_plot_list, ncol = n_cols, nrow = n_rows) + 
        plot_annotation(title = paste("QC Metrics by", cluster_col),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
      
      # Save combined QC plot
      base_name <- basename(args$project_name)
      pdf(paste0(base_name, args$suffix, "_", cluster_col, "_QC_plots.pdf"), 
          width = 8 * n_cols, height = 6 * n_rows)
      print(combined_qc_plot)
      dev.off()
      
      # Also create individual PDFs for each QC metric
      for(qc_metric in names(qc_plot_list)) {
        pdf(paste0(base_name, args$suffix, "_", cluster_col, "_", qc_metric, "_QC.pdf"), 
            width = 12, height = 8)
        print(qc_plot_list[[qc_metric]])
        dev.off()
      }
      
      cat(sprintf("    Created QC plots for %s\n", cluster_col))
    } else {
      cat(sprintf("    No QC plots created for %s - no valid data\n", cluster_col))
    }
  }
} else {
  cat("Warning: No cluster columns or QC metrics found for QC plots\n")
  if(length(cluster_columns) == 0) {
    cat("  Available cluster columns:", colnames(cell_col_data)[grep("Cluster", colnames(cell_col_data))], "\n")
  }
  if(length(qc_metrics) == 0) {
    cat("  Available QC metrics:", colnames(cell_col_data)[colnames(cell_col_data) %in% 
           c("TSSEnrichment", "nFrags", "BlacklistRatio", "NucleosomeRatio",
             "log10(nFrags)", "log10(nFrags+1)", "DoubletScore", "DoubletEnrichment")], "\n")
  }
}

cat("=== All PDF plots created successfully ===\n")
cat("Output files saved with base name:", basename(args$project_name), "\n")

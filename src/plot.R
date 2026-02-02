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

# Define embeddings to check
embeddings <- c("UMAP_RNA", "UMAP_ATAC", "UMAP_Combined")
embedding_names <- c("RNA", "ATAC", "Combined")

# 1. UMAP PER SAMPLE - All samples together, colored by sample
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating UMAP plots colored by sample...\n")
  
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
      pdf(paste0(base_name, args$suffix, "_", embedding_names[i], "_by_sample.pdf"), 
          width = 12, height = 8)
      print(p_sample)
      dev.off()
      
      cat(sprintf("  Created %s UMAP colored by sample\n", embedding_names[i]))
    } else {
      cat(sprintf("  Skipping %s embedding - not found\n", embedding_names[i]))
    }
  }
} else {
  cat("Warning: 'Sample' column not found in cellColData. Skipping UMAP by sample plots.\n")
}

# 2. UMAP PER SAMPLE PER MODALITY - Separate plots for each sample
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating separate UMAP plots for each sample...\n")
  
  samples <- unique(cell_col_data$Sample)
  
  for(embed_idx in seq_along(embeddings)) {
    if(embedding_exists(proj_ALL, embeddings[embed_idx])) {
      # Create a list to store plots for each sample
      sample_plots <- list()
      
      for(sample in samples) {
        # Get cell indices for this sample
        sample_cells <- which(cell_col_data$Sample == sample)
        proj_sample <- proj_ALL[sample_cells, ]
        
        # Try different cluster columns
        cluster_cols_to_try <- c("Clusters_Combined", "Clusters_RNA", "Clusters_ATAC")
        cluster_col <- NULL
        
        for(cc in cluster_cols_to_try) {
          if(cc %in% colnames(cell_col_data)) {
            cluster_col <- cc
            break
          }
        }
        
        if(is.null(cluster_col)) {
          cluster_col <- "Sample"  # Fallback to sample if no clusters
        }
        
        # Create UMAP for this sample
        p <- plotEmbedding(
          ArchRProj = proj_sample,
          colorBy = "cellColData",
          name = cluster_col,
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
        plot_annotation(title = paste(embedding_names[embed_idx], "UMAP - Separate by Sample"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
      
      # Save combined plot
      base_name <- basename(args$project_name)
      pdf(paste0(base_name, args$suffix, "_", embedding_names[embed_idx], 
                 "_per_sample_separate.pdf"), 
          width = 6 * n_cols, height = 5 * n_rows)
      print(combined_plot)
      dev.off()
      
      cat(sprintf("  Created %s UMAP separate plots for %d samples\n", 
                  embedding_names[embed_idx], length(samples)))
    } else {
      cat(sprintf("  Skipping %s embedding - not found\n", embedding_names[embed_idx]))
    }
  }
}

# 3. STANDARD UMAPS - Colored by clusters (the original plots)
cat("Creating standard cluster-colored UMAP plots...\n")

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
    pdf(paste0(base_name, args$suffix, "_", modality, "_clusters_UMAP.pdf"), width = 12, height = 8)
    print(p)
    dev.off()
    
    cat(sprintf("  Created %s UMAP colored by %s clusters\n", modality, modality))
  } else {
    if(!embedding_exists(proj_ALL, embedding_name)) {
      cat(sprintf("  Skipping %s UMAP - embedding not found\n", modality))
    } else if(!cluster_name %in% colnames(cell_col_data)) {
      cat(sprintf("  Skipping %s UMAP - %s not found in cellColData\n", modality, cluster_name))
    }
  }
}

# 4. QC PLOTS PER CLUSTERS
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
cat("Output files will include:\n")
cat("  1. *_RNA_by_sample.pdf - All cells colored by sample\n")
cat("  2. *_RNA_per_sample_separate.pdf - Each sample separately\n")
cat("  3. *_RNA_clusters_UMAP.pdf - Colored by RNA clusters\n")
cat("  4. QC plots for each cluster type\n")
cat("Output files saved with base name:", basename(args$project_name), "\n")

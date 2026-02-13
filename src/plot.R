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

# 1. UMAP per sample (with NEW file names)
if("Sample" %in% colnames(cell_col_data)) {
  cat("Creating UMAP plots per sample...\n")
  
  # RNA UMAP colored by sample
  if(embedding_exists(proj_ALL, "UMAP_RNA")) {
    p_sample_rna <- plotEmbedding(
      ArchRProj = proj_ALL,
      colorBy = "cellColData",
      name = "Sample",
      embedding = "UMAP_RNA",
      plotAs = "points",
      size = 0.5
    )
    
    pdf(paste0(base_name, args$suffix, "_RNA_SAMPLE_UMAP.pdf"), width = 12, height = 8)
    print(p_sample_rna)
    dev.off()
    cat("  Created RNA_SAMPLE_UMAP.pdf\n")
  }
  
  # ATAC UMAP colored by sample
  if(embedding_exists(proj_ALL, "UMAP_ATAC")) {
    p_sample_atac <- plotEmbedding(
      ArchRProj = proj_ALL,
      colorBy = "cellColData",
      name = "Sample",
      embedding = "UMAP_ATAC",
      plotAs = "points",
      size = 0.5
    )
    
    pdf(paste0(base_name, args$suffix, "_ATAC_SAMPLE_UMAP.pdf"), width = 12, height = 8)
    print(p_sample_atac)
    dev.off()
    cat("  Created ATAC_SAMPLE_UMAP.pdf\n")
  }
  
  # Combined UMAP colored by sample
  if(embedding_exists(proj_ALL, "UMAP_Combined")) {
    p_sample_combined <- plotEmbedding(
      ArchRProj = proj_ALL,
      colorBy = "cellColData",
      name = "Sample",
      embedding = "UMAP_Combined",
      plotAs = "points",
      size = 0.5
    )
    
    pdf(paste0(base_name, args$suffix, "_Combined_SAMPLE_UMAP.pdf"), width = 12, height = 8)
    print(p_sample_combined)
    dev.off()
    cat("  Created Combined_SAMPLE_UMAP.pdf\n")
  }
}

# 2. UMAP per sample for COMBINED embedding only
if("Sample" %in% colnames(cell_col_data) && embedding_exists(proj_ALL, "UMAP_Combined")) {
  cat("Creating UMAP per sample for Combined embedding...\n")
  
  samples <- unique(cell_col_data$Sample)
  sample_plots <- list()
  
  for(sample in samples) {
    sample_cells <- which(cell_col_data$Sample == sample)
    proj_sample <- proj_ALL[sample_cells, ]
    
    p <- plotEmbedding(
      ArchRProj = proj_sample,
      colorBy = "cellColData",
      name = "Clusters_Combined",
      embedding = "UMAP_Combined",
      plotAs = "points",
      size = 1.0
    ) + ggtitle(sample)
    
    sample_plots[[sample]] <- p
  }
  
  n_cols <- min(3, length(samples))
  n_rows <- ceiling(length(samples) / n_cols)
  
  combined_plot <- wrap_plots(sample_plots, ncol = n_cols, nrow = n_rows)
  
  pdf(paste0(base_name, args$suffix, "_Combined_SAMPLES_GRID.pdf"),
      width = 6 * n_cols, height = 5 * n_rows)
  print(combined_plot)
  dev.off()
  
  cat("  Created Combined_SAMPLES_GRID.pdf\n")
}

# 3. BOTH SAMPLES TOGETHER with clusters for annotation
if("Sample" %in% colnames(cell_col_data) && embedding_exists(proj_ALL, "UMAP_Combined")) {
  cat("Creating both samples together with clusters...\n")
  
  # Create plot with both samples, colored by clusters
  p_both_samples <- plotEmbedding(
    ArchRProj = proj_ALL,
    colorBy = "cellColData",
    name = "Clusters_Combined",
    embedding = "UMAP_Combined",
    plotAs = "points",
    size = 0.5
  )
  
  pdf(paste0(base_name, args$suffix, "_BothSamples_CombinedClusters_UMAP.pdf"), width = 12, height = 8)
  print(p_both_samples)
  dev.off()
  
  cat("  Created BothSamples_CombinedClusters_UMAP.pdf\n")
}

# 4. FIXED: QC plots per clusters - COMBINED NOW SHOWS BOTH ATAC AND RNA METRICS
cat("Creating QC plots per clusters with proper assay-specific metrics...\n")

# ATAC-specific QC metrics
qc_metrics_atac <- c("TSSEnrichment", "nFrags", "BlacklistRatio", "NucleosomeRatio")
qc_metrics_atac <- qc_metrics_atac[qc_metrics_atac %in% colnames(cell_col_data)]

# RNA QC metrics
qc_metrics_rna <- c("Gex_nUMI", "Gex_nGenes", "Gex_Log10_nUMI", "Gex_Log10_nGenes")
qc_metrics_rna <- qc_metrics_rna[qc_metrics_rna %in% colnames(cell_col_data)]

# For Combined clusters - show BOTH ATAC and RNA QC metrics
if("Clusters_Combined" %in% colnames(cell_col_data)) {
  qc_plot_list <- list()
  
  # Add ATAC-specific QC plots
  for(qc_metric in qc_metrics_atac) {
    df <- data.frame(
      Cluster = cell_col_data[["Clusters_Combined"]],
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
        labs(title = paste(qc_metric, "(ATAC)"), 
             y = qc_metric, x = "Cluster")
      qc_plot_list[[qc_metric]] <- p
    }
  }
  
  # Add RNA QC plots
  for(qc_metric in qc_metrics_rna) {
    df <- data.frame(
      Cluster = cell_col_data[["Clusters_Combined"]],
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
        labs(title = paste(qc_metric, "(RNA)"), 
             y = qc_metric, x = "Cluster")
      qc_plot_list[[qc_metric]] <- p
    }
  }
  
  if(length(qc_plot_list) > 0) {
    n_cols <- min(2, length(qc_plot_list))
    n_rows <- ceiling(length(qc_plot_list) / n_cols)
    
    combined_qc_plot <- wrap_plots(qc_plot_list, ncol = n_cols, nrow = n_rows)
    
    pdf(paste0(base_name, args$suffix, "_Clusters_Combined_QC.pdf"),
        width = 8 * n_cols, height = 6 * n_rows)
    print(combined_qc_plot)
    dev.off()
    cat("  Created Clusters_Combined_QC.pdf with ATAC and RNA metrics\n")
  }
}

# For RNA clusters - RNA QC metrics
if("Clusters_RNA" %in% colnames(cell_col_data) && length(qc_metrics_rna) > 0) {
  qc_plot_list <- list()
  
  for(qc_metric in qc_metrics_rna) {
    df <- data.frame(
      Cluster = cell_col_data[["Clusters_RNA"]],
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
        labs(title = qc_metric, y = qc_metric, x = "Cluster")
      qc_plot_list[[qc_metric]] <- p
    }
  }
  
  if(length(qc_plot_list) > 0) {
    n_cols <- min(2, length(qc_plot_list))
    n_rows <- ceiling(length(qc_plot_list) / n_cols)
    
    combined_qc_plot <- wrap_plots(qc_plot_list, ncol = n_cols, nrow = n_rows)
    
    pdf(paste0(base_name, args$suffix, "_Clusters_RNA_QC.pdf"),
        width = 8 * n_cols, height = 6 * n_rows)
    print(combined_qc_plot)
    dev.off()
    cat("  Created Clusters_RNA_QC.pdf\n")
  }
}

# For ATAC clusters - only ATAC-specific QC metrics
if("Clusters_ATAC" %in% colnames(cell_col_data) && length(qc_metrics_atac) > 0) {
  qc_plot_list <- list()
  
  for(qc_metric in qc_metrics_atac) {
    df <- data.frame(
      Cluster = cell_col_data[["Clusters_ATAC"]],
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
        labs(title = qc_metric, y = qc_metric, x = "Cluster")
      qc_plot_list[[qc_metric]] <- p
    }
  }
  
  if(length(qc_plot_list) > 0) {
    n_cols <- min(2, length(qc_plot_list))
    n_rows <- ceiling(length(qc_plot_list) / n_cols)
    
    combined_qc_plot <- wrap_plots(qc_plot_list, ncol = n_cols, nrow = n_rows)
    
    pdf(paste0(base_name, args$suffix, "_Clusters_ATAC_QC.pdf"),
        width = 8 * n_cols, height = 6 * n_rows)
    print(combined_qc_plot)
    dev.off()
    cat("  Created Clusters_ATAC_QC.pdf\n")
  }
}

cat("=== ALL PLOTS CREATED SUCCESSFULLY ===\n")
cat("Original files (unchanged):\n")
cat("  ", base_name, args$suffix, "_RNA_UMAP.pdf\n", sep="")
cat("  ", base_name, args$suffix, "_ATAC_UMAP.pdf\n", sep="")
cat("  ", base_name, args$suffix, "_Combined_UMAP.pdf\n", sep="")
cat("\nNew files added:\n")
cat("  *_SAMPLE_UMAP.pdf (colored by sample)\n")
cat("  Combined_SAMPLES_GRID.pdf (each sample separately - Combined only)\n")
cat("  BothSamples_CombinedClusters_UMAP.pdf (both samples with clusters)\n")
cat("  Clusters_Combined_QC.pdf (ATAC and RNA metrics)\n")
cat("  Clusters_ATAC_QC.pdf (ATAC metrics only)\n")
cat("  Clusters_RNA_QC.pdf (RNA metrics only)\n")

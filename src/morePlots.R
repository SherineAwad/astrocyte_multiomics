#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ArchR)
  library(argparse)
  library(pheatmap)
})

# -----------------------------
# Arguments
# -----------------------------
parser <- ArgumentParser()
parser$add_argument("--project_name", required = TRUE, help="Path to ArchR project")
parser$add_argument("--top_n", type="integer", default=50, help="Number of top peaks to plot")
args <- parser$parse_args()

proj_dir <- args$project_name
top_n <- args$top_n

# FORCE top_n to be numeric and correct
top_n <- as.integer(top_n)
if(is.na(top_n) || top_n < 1) top_n <- 50

cat("========================================\n")
cat("ArchR project directory:", proj_dir, "\n")
cat("TOP_N PARAMETER:", top_n, "peaks will be plotted\n")
cat("========================================\n\n")

# -----------------------------
# Load project
# -----------------------------
cat("Loading ArchR project...\n")
proj <- loadArchRProject(proj_dir)
cat("Project loaded! Cells:", nCells(proj), "\n")

# -----------------------------
# Cell type order
# -----------------------------
cat("\nSetting cell type order...\n")
cell_order_requested <- c("WT_Astro", "KO_Astro", "Diff_Astro", "NSC", "pro_NSC", "Neurons", "Oligos")
cell_types <- getCellColData(proj)$CellType
names(cell_types) <- rownames(getCellColData(proj))
cell_order <- cell_order_requested[cell_order_requested %in% unique(cell_types)]
cat("Cell order:", paste(cell_order, collapse=" -> "), "\n")

# -----------------------------
# Load peak annotations
# -----------------------------
cat("\nLoading peak-to-gene annotations...\n")
peak_anno_file <- file.path(proj_dir, "peak_nearest_gene_annotation.csv")
if(!file.exists(peak_anno_file)) stop("Peak annotation file not found!")
peak2gene <- read.csv(peak_anno_file)
rownames(peak2gene) <- peak2gene$peak
cat("Loaded", nrow(peak2gene), "peak-gene annotations\n")

# -----------------------------
# Get peak names with gene symbols
# -----------------------------
cat("\nGetting peak names...\n")
peak_set <- getPeakSet(proj)
all_peaks <- paste0(seqnames(peak_set), ":", start(peak_set), "-", end(peak_set))
names(all_peaks) <- all_peaks
cat("Total peaks in project:", length(all_peaks), "\n")

# Match peaks to genes
peak_gene_map <- peak2gene[all_peaks, "nearest_gene", drop=FALSE]
peak_with_gene <- paste0(all_peaks, " (", peak_gene_map$nearest_gene, ")")
names(peak_with_gene) <- all_peaks

# ========================================
# CRITICAL FIX - FORCE top_n TO WORK
# ========================================
cat("\n========================================\n")
cat("SELECTING PEAKS - USING top_n =", top_n, "\n")
cat("========================================\n")

# Method 1: Direct indexing (FORCED)
peak_indices <- 1:min(top_n, length(all_peaks))
top_peaks <- all_peaks[peak_indices]
top_peaks_with_gene <- peak_with_gene[top_peaks]

cat("Selected", length(top_peaks), "peaks for plotting\n")
cat("First 5 peaks:\n")
for(i in 1:min(5, length(top_peaks_with_gene))) {
  cat("  ", top_peaks_with_gene[i], "\n")
}
cat("Last 5 peaks:\n")
for(i in max(1, length(top_peaks_with_gene)-4):length(top_peaks_with_gene)) {
  cat("  ", top_peaks_with_gene[i], "\n")
}
cat("========================================\n\n")

# -----------------------------
# Colors
# -----------------------------
cell_colors <- list(
  CellType = c(
    "WT_Astro" = "#1f77b4", "KO_Astro" = "#ff7f0e", "Diff_Astro" = "#2ca02c",
    "NSC" = "#d62728", "pro_NSC" = "#9467bd", "Neurons" = "#8c564b", "Oligos" = "#e377c2"
  )
)

cell_annotation <- data.frame(
  CellType = cell_order,
  row.names = cell_order
)

# -----------------------------
# 1. PEAK HEATMAP
# -----------------------------
cat("\n1. Creating peak heatmap with", length(top_peaks), "peaks...\n")

peak_matrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
peak_data <- assay(peak_matrix)
rownames(peak_data) <- all_peaks
colnames(peak_data) <- getCellNames(proj)

# Aggregate by cell type
peak_agg <- matrix(0, nrow = length(top_peaks), ncol = length(cell_order))
rownames(peak_agg) <- top_peaks_with_gene
colnames(peak_agg) <- cell_order

for(j in seq_along(cell_order)) {
  ct <- cell_order[j]
  cells_ct <- names(cell_types)[cell_types == ct]
  cells_ct <- head(cells_ct, min(100, length(cells_ct)))
  if(length(cells_ct) > 0) {
    cell_idx <- which(colnames(peak_data) %in% cells_ct)
    peak_agg[, j] <- rowMeans(peak_data[top_peaks, cell_idx, drop=FALSE], na.rm=TRUE)
  }
}

pdf(file.path(proj_dir, "1_peak_heatmap.pdf"), width=12, height=max(10, nrow(peak_agg)*0.25))
pheatmap(
  peak_agg,
  annotation_col = cell_annotation,
  annotation_colors = cell_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = paste("Top", nrow(peak_agg), "Peaks with Nearest Genes"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  fontsize_row = 8
)
dev.off()
cat("  Saved: 1_peak_heatmap.pdf with", nrow(peak_agg), "peaks\n")

# -----------------------------
# 2. MOTIF HEATMAP
# -----------------------------
cat("\n2. Creating motif heatmap...\n")

if("MotifMatrix" %in% getAvailableMatrices(proj)) {
  motif_matrix <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
  motif_data <- assay(motif_matrix)
  
  # Get motif names
  if("name" %in% colnames(rowData(motif_matrix))) {
    rownames(motif_data) <- rowData(motif_matrix)$name
  } else {
    rownames(motif_data) <- paste0("Motif_", 1:nrow(motif_data))
  }
  colnames(motif_data) <- getCellNames(proj)
  
  # Aggregate by cell type
  motif_agg <- matrix(0, nrow = nrow(motif_data), ncol = length(cell_order))
  rownames(motif_agg) <- rownames(motif_data)
  colnames(motif_agg) <- cell_order
  
  for(j in seq_along(cell_order)) {
    ct <- cell_order[j]
    cells_ct <- names(cell_types)[cell_types == ct]
    cells_ct <- head(cells_ct, min(100, length(cells_ct)))
    if(length(cells_ct) > 0) {
      cell_idx <- which(colnames(motif_data) %in% cells_ct)
      motif_agg[, j] <- rowMeans(motif_data[, cell_idx, drop=FALSE], na.rm=TRUE)
    }
  }
  
  # Remove rows with NA/NaN/Inf
  motif_agg <- motif_agg[complete.cases(motif_agg) & 
                         apply(motif_agg, 1, function(x) all(is.finite(x))), ]
  
  # Get top variable motifs
  motif_var <- apply(motif_agg, 1, var, na.rm=TRUE)
  motif_var <- motif_var[is.finite(motif_var) & !is.na(motif_var)]
  top_motifs <- names(sort(motif_var, decreasing = TRUE))[1:min(50, length(motif_var))]
  
  if(length(top_motifs) > 0) {
    pdf(file.path(proj_dir, "2_motif_heatmap.pdf"), width=12, height=10)
    pheatmap(
      motif_agg[top_motifs, , drop=FALSE],
      annotation_col = cell_annotation,
      annotation_colors = cell_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      main = "Top 50 Variable Motifs",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      scale = "row",
      fontsize_row = 8
    )
    dev.off()
    cat("  Saved: 2_motif_heatmap.pdf\n")
  } else {
    cat("  No valid motifs found\n")
  }
} else {
  cat("  No MotifMatrix found\n")
}

# -----------------------------
# 3. TF ACTIVITY HEATMAP
# -----------------------------
cat("\n3. Creating TF activity heatmap...\n")

if("GeneScoreMatrix" %in% getAvailableMatrices(proj)) {
  gene_matrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
  gene_data <- assay(gene_matrix)
  rownames(gene_data) <- rowData(gene_matrix)$name
  colnames(gene_data) <- getCellNames(proj)
  
  # Get TF list from motif names
  tf_list <- NULL
  if(exists("motif_data")) {
    tf_list <- unique(gsub("_.*", "", rownames(motif_data)))
    tf_list <- tf_list[tf_list %in% rownames(gene_data)]
    cat("  Found", length(tf_list), "TFs from motif annotations\n")
  }
  
  if(!is.null(tf_list) && length(tf_list) > 10) {
    use_genes <- tf_list
  } else {
    # Use top variable genes
    gene_var <- apply(gene_data, 1, var, na.rm=TRUE)
    gene_var <- gene_var[is.finite(gene_var) & !is.na(gene_var)]
    use_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, length(gene_var))]
    cat("  Using top 50 variable genes\n")
  }
  
  # Aggregate by cell type
  tf_agg <- matrix(0, nrow = length(use_genes), ncol = length(cell_order))
  rownames(tf_agg) <- use_genes
  colnames(tf_agg) <- cell_order
  
  for(j in seq_along(cell_order)) {
    ct <- cell_order[j]
    cells_ct <- names(cell_types)[cell_types == ct]
    cells_ct <- head(cells_ct, min(100, length(cells_ct)))
    if(length(cells_ct) > 0) {
      cell_idx <- which(colnames(gene_data) %in% cells_ct)
      tf_agg[, j] <- rowMeans(gene_data[use_genes, cell_idx, drop=FALSE], na.rm=TRUE)
    }
  }
  
  # Remove problematic rows
  row_var <- apply(tf_agg, 1, var, na.rm=TRUE)
  keep_rows <- is.finite(row_var) & !is.na(row_var) & row_var > 1e-10
  tf_agg <- tf_agg[keep_rows, , drop=FALSE]
  tf_agg <- tf_agg[complete.cases(tf_agg), , drop=FALSE]
  
  if(nrow(tf_agg) > 1) {
    if(nrow(tf_agg) > 50) {
      row_var <- apply(tf_agg, 1, var)
      tf_agg <- tf_agg[names(sort(row_var, decreasing=TRUE))[1:50], ]
    }
    
    pdf(file.path(proj_dir, "3_tf_activity_heatmap.pdf"), width=12, height=10)
    pheatmap(
      tf_agg,
      annotation_col = cell_annotation,
      annotation_colors = cell_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      main = paste("TF Activity -", nrow(tf_agg), "TFs"),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      scale = "row",
      fontsize_row = 8
    )
    dev.off()
    cat("  Saved: 3_tf_activity_heatmap.pdf (", nrow(tf_agg), "TFs)\n")
  } else {
    cat("  Not enough valid TFs for clustering\n")
  }
} else {
  cat("  No GeneScoreMatrix found\n")
}

# -----------------------------
# 4. COMBINED PEAK-MOTIF HEATMAP
# -----------------------------
cat("\n4. Creating combined peak-motif heatmap...\n")

if(exists("motif_data") && exists("peak_agg")) {
  if(ncol(motif_data) == length(all_peaks)) {
    colnames(motif_data) <- all_peaks
    motif_peak_data <- motif_data[, top_peaks, drop=FALSE]
    
    motif_peak_data <- motif_peak_data[complete.cases(motif_peak_data) & 
                                        apply(motif_peak_data, 1, function(x) all(is.finite(x))), ]
    
    if(nrow(motif_peak_data) > 0) {
      motif_peak_var <- apply(motif_peak_data, 1, var, na.rm=TRUE)
      motif_peak_var <- motif_peak_var[is.finite(motif_peak_var) & !is.na(motif_peak_var)]
      top_motifs_peaks <- names(sort(motif_peak_var, decreasing = TRUE))[1:min(30, length(motif_peak_var))]
      
      if(length(top_motifs_peaks) > 1) {
        motif_peak_subset <- motif_peak_data[top_motifs_peaks, , drop=FALSE]
        colnames(motif_peak_subset) <- top_peaks_with_gene[colnames(motif_peak_subset)]
        
        pdf(file.path(proj_dir, "4_combined_peak_motif.pdf"), width=14, height=10)
        pheatmap(
          motif_peak_subset,
          show_rownames = TRUE,
          show_colnames = TRUE,
          main = paste("Motif Enrichment in Top", length(top_peaks), "Peaks"),
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          scale = "row",
          fontsize_row = 8,
          fontsize_col = 8
        )
        dev.off()
        cat("  Saved: 4_combined_peak_motif.pdf\n")
      }
    }
  }
}

cat("\n========================================\n")
cat("✅ ALL HEATMAPS COMPLETE!\n")
cat("========================================\n")
cat("Files created:\n")
cat("  - 1_peak_heatmap.pdf (", nrow(peak_agg), "peaks with gene names)\n")
cat("  - 2_motif_heatmap.pdf\n")
cat("  - 3_tf_activity_heatmap.pdf\n")
cat("  - 4_combined_peak_motif.pdf\n")
cat("========================================\n")

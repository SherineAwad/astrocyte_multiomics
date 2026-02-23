suppressPackageStartupMessages({
  library(ArchR)
  library(argparse)
  library(ComplexHeatmap)
  library(circlize)
})

# -----------------------------
# Arguments
# -----------------------------
parser <- ArgumentParser()
parser$add_argument("--project_name", required = TRUE)
parser$add_argument("--top_peaks", default = 300, type = "integer")
parser$add_argument("--top_TFs", default = 50, type = "integer")
args <- parser$parse_args()

proj_dir <- args$project_name
top_peaks <- args$top_peaks
top_TFs <- args$top_TFs

# -----------------------------
# Load project
# -----------------------------
cat("Loading ArchR project...\n")
proj <- loadArchRProject(path = proj_dir)
cat("ArchR project loaded successfully!\n")

# -----------------------------
# Cell type order (DESIRED)
# -----------------------------
cell_order_requested <- c(
  "WT_Astro",
  "KO_Astro",
  "Diff_Astro",
  "NSC",
  "pro_NSC",
  "Neurons",
  "Oligos"
)

# ============================================================
# Helper function to get available cell types
# ============================================================
get_available_celltypes <- function(proj, requested) {
  available <- getCellColData(proj)$CellType %>% unique() %>% as.character()
  intersect(requested, available)
}

# ============================================================
# 1) MARKER PEAKS HEATMAP (TOP N)
# ============================================================
cat("Plotting marker peaks heatmap...\n")

marker_file <- file.path(proj_dir, "markerPeaks_object.rds")
if(file.exists(marker_file)) {
  markers <- readRDS(marker_file)
  
  # Plot and save marker peaks heatmap
  tryCatch({
    hm <- plotMarkerHeatmap(
      seMarker = markers,
      cutOff = "FDR <= 0.01 & Log2FC >= 1",
      labelMarkers = NULL,
      transpose = FALSE
    )
    
    png(
      filename = file.path(proj_dir, "MarkerPeaks_Heatmap_Top.png"),
      width = 2200,
      height = 2000,  # Increased height
      res = 200
    )
    # Draw with better margins
    draw(hm, 
         heatmap_legend_side = "right",
         padding = unit(c(2, 2, 2, 2), "cm"))  # Add padding
    dev.off()
    cat("✓ Marker peaks heatmap saved\n")
  }, error = function(e) {
    cat("Error in marker peaks heatmap: ", e$message, "\n")
  })
} else {
  cat("Warning: Marker peaks file not found\n")
}

# ============================================================
# 2) MOTIF ENRICHMENT HEATMAP
# ============================================================
cat("Plotting motif enrichment heatmap...\n")

motif_file <- file.path(proj_dir, "motifEnrichment.rds")
if(file.exists(motif_file)) {
  tryCatch({
    # Read motif enrichment data
    motif_data <- readRDS(motif_file)
    motif_mat <- as.matrix(assay(motif_data))
    
    # Get available cell types
    available_cells <- intersect(cell_order_requested, colnames(motif_mat))
    
    if(length(available_cells) > 0) {
      # Subset matrix
      motif_mat <- motif_mat[, available_cells, drop = FALSE]
      
      # Remove rows with all zeros or NA
      motif_mat <- motif_mat[rowSums(is.na(motif_mat)) == 0, , drop = FALSE]
      motif_mat <- motif_mat[rowSums(motif_mat) != 0, , drop = FALSE]
      
      if(nrow(motif_mat) > 0) {
        # Select top motifs by variance (if too many)
        if(nrow(motif_mat) > 100) {
          motif_var <- apply(motif_mat, 1, var)
          top_motifs <- names(sort(motif_var, decreasing = TRUE))[1:100]
          motif_mat <- motif_mat[top_motifs, , drop = FALSE]
        }
        
        png(
          filename = file.path(proj_dir, "MotifEnrichment_Heatmap.png"),
          width = 2400,  # Wider
          height = 2800,  # Taller
          res = 200
        )
        
        # Create heatmap with better formatting
        ht <- Heatmap(
          motif_mat,
          name = "Motif Enrichment\n(z-score)",
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          show_row_names = FALSE,  # Hide row names if too many
          show_row_dend = TRUE,
          column_names_rot = 45,  # Rotate column names
          column_names_gp = gpar(fontsize = 10),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10),
            labels_gp = gpar(fontsize = 8)
          ),
          border = TRUE
        )
        
        draw(ht, 
             padding = unit(c(2, 2, 2, 4), "cm"),  # Add right padding for legend
             heatmap_legend_side = "right")
        
        dev.off()
        cat("✓ Motif enrichment heatmap saved\n")
      } else {
        cat("Warning: No valid motifs after filtering\n")
      }
    } else {
      cat("Warning: No requested cell types found in motif matrix\n")
    }
  }, error = function(e) {
    cat("Error in motif enrichment: ", e$message, "\n")
    if(length(dev.list()) > 0) dev.off()
  })
} else {
  cat("Warning: Motif enrichment file not found\n")
}

# ============================================================
# 3) TF ACTIVITY HEATMAP (CORRECTED)
# ============================================================
cat("Computing TF activity heatmap...\n")

tryCatch({
  # Correct way to get motif matrix in ArchR
  motif_matrix <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "MotifMatrix"
  )
  
  # Get the deviation (z-score) matrix - this is the TF activity
  # In ArchR, the MotifMatrix contains "deviations" (z-scores) and "matrix" (raw)
  if("deviations" %in% names(assays(motif_matrix))) {
    tf_mat <- assay(motif_matrix, "deviations")
  } else {
    tf_mat <- assay(motif_matrix)  # Fallback to first assay
  }
  
  # Get cell types from project
  cell_types <- getCellColData(proj)$CellType
  available_cells <- intersect(cell_order_requested, unique(cell_types))
  
  if(length(available_cells) > 0) {
    # Calculate mean TF activity per cell type
    tf_by_celltype <- matrix(0, 
                             nrow = nrow(tf_mat), 
                             ncol = length(available_cells),
                             dimnames = list(rownames(tf_mat), available_cells))
    
    for(ct in available_cells) {
      cells_ct <- which(cell_types == ct)
      if(length(cells_ct) > 0) {
        tf_by_celltype[, ct] <- rowMeans(tf_mat[, cells_ct, drop = FALSE], na.rm = TRUE)
      }
    }
    
    # Remove rows with all NA/zero
    tf_by_celltype <- tf_by_celltype[ rowSums(is.na(tf_by_celltype)) == 0, , drop = FALSE]
    tf_by_celltype <- tf_by_celltype[ rowSums(tf_by_celltype) != 0, , drop = FALSE]
    
    if(nrow(tf_by_celltype) > 0) {
      # Select top TFs by variance
      tf_var <- apply(tf_by_celltype, 1, var, na.rm = TRUE)
      tf_var <- tf_var[!is.na(tf_var) & tf_var > 0]
      
      if(length(tf_var) > 0) {
        n_top <- min(top_TFs, length(tf_var))
        top_tf_names <- names(sort(tf_var, decreasing = TRUE))[1:n_top]
        tf_mat_final <- tf_by_celltype[top_tf_names, , drop = FALSE]
        
        # Plot TF activity heatmap
        png(
          filename = file.path(proj_dir, "TFactivity_Heatmap.png"),
          width = 2400,
          height = 2800,
          res = 200
        )
        
        ht <- Heatmap(
          tf_mat_final,
          name = "TF Activity\n(z-score)",
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          row_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 10),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10),
            labels_gp = gpar(fontsize = 8)
          ),
          border = TRUE
        )
        
        draw(ht, 
             padding = unit(c(2, 2, 2, 4), "cm"),
             heatmap_legend_side = "right")
        
        dev.off()
        cat("✓ TF activity heatmap saved\n")
      } else {
        cat("Warning: No variance in TF activity\n")
      }
    } else {
      cat("Warning: TF activity matrix empty after filtering\n")
    }
  } else {
    cat("Warning: No requested cell types found\n")
  }
}, error = function(e) {
  cat("Error in TF activity: ", e$message, "\n")
  if(length(dev.list()) > 0) dev.off()
})

cat("\n=== ALL PLOTS COMPLETED ===\n")
cat("Check the project directory for:\n")
cat("  - MarkerPeaks_Heatmap_Top.png\n")
cat("  - MotifEnrichment_Heatmap.png\n")
cat("  - TFactivity_Heatmap.png\n")

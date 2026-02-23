suppressPackageStartupMessages({
  library(ArchR)
  library(argparse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
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

cat("Project directory:", proj_dir, "\n")
cat("Top peaks requested:", top_peaks, "\n")
cat("Top TFs requested:", top_TFs, "\n")

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
# Helper function to get available cell types in correct order
# ============================================================
get_available_celltypes <- function(proj, requested) {
  available <- unique(getCellColData(proj)$CellType)
  requested[requested %in% available]
}

available_celltypes <- get_available_celltypes(proj, cell_order_requested)
cat("Available cell types in order:", paste(available_celltypes, collapse = " -> "), "\n")

# ============================================================
# 1) MARKER PEAKS HEATMAP - WITH TOP_PEAKS LIMIT
# ============================================================
cat("\nPlotting marker peaks heatmap...\n")

marker_file <- file.path(proj_dir, "markerPeaks_object.rds")
if(file.exists(marker_file)) {
  markers <- readRDS(marker_file)
  
  tryCatch({
    # Get the marker peaks matrix
    marker_mat <- assay(markers)
    
    # Filter for significant markers first
    # Get marker indices that pass cutoff
    marker_list <- lapply(markers, function(x) {
      if(!is.null(x$FDR) && !is.null(x$Log2FC)) {
        which(x$FDR <= 0.01 & x$Log2FC >= 1)
      } else {
        integer(0)
      }
    })
    
    # Get unique significant peaks across all cell types
    sig_peaks <- unique(unlist(lapply(names(marker_list), function(ct) {
      if(length(marker_list[[ct]]) > 0) {
        rownames(markers[[ct]])[marker_list[[ct]]]
      }
    })))
    
    cat("Significant peaks found:", length(sig_peaks), "\n")
    
    # Apply top_peaks limit
    if(length(sig_peaks) > top_peaks) {
      # Calculate variance for significant peaks only
      sig_mat <- marker_mat[rownames(marker_mat) %in% sig_peaks, , drop = FALSE]
      peak_var <- apply(sig_mat, 1, var, na.rm = TRUE)
      top_peaks_use <- names(sort(peak_var, decreasing = TRUE))[1:top_peaks]
      cat("Selecting top", top_peaks, "peaks by variance\n")
    } else {
      top_peaks_use <- sig_peaks
      cat("Using all", length(sig_peaks), "significant peaks\n")
    }
    
    # Create a custom marker subset
    markers_subset <- markers
    for(ct in names(markers_subset)) {
      if(!is.null(markers_subset[[ct]])) {
        # Keep only top peaks in each cell type's marker set
        keep <- rownames(markers_subset[[ct]]) %in% top_peaks_use
        markers_subset[[ct]] <- markers_subset[[ct]][keep, , drop = FALSE]
      }
    }
    
    # Plot with subset
    hm <- plotMarkerHeatmap(
      seMarker = markers_subset,
      cutOff = "FDR <= 0.01 & Log2FC >= 1",
      labelMarkers = NULL,
      transpose = FALSE
    )
    
    png(
      filename = file.path(proj_dir, "MarkerPeaks_Heatmap_Top.png"),
      width = 2200,
      height = min(3000, 800 + length(top_peaks_use) * 15),  # Dynamic height
      res = 200
    )
    draw(hm,
         heatmap_legend_side = "right",
         padding = unit(c(2, 2, 2, 2), "cm"))
    dev.off()
    cat("✓ Marker peaks heatmap saved with", length(top_peaks_use), "peaks\n")
    
  }, error = function(e) {
    cat("Error in marker peaks heatmap: ", e$message, "\n")
    if(length(dev.list()) > 0) dev.off()
  })
} else {
  cat("Warning: Marker peaks file not found\n")
}

# ============================================================
# 2) MOTIF ENRICHMENT HEATMAP - WITH ROW NAMES FIXED
# ============================================================
cat("\nPlotting motif enrichment heatmap...\n")

motif_file <- file.path(proj_dir, "motifEnrichment.rds")
if(file.exists(motif_file)) {
  tryCatch({
    motif_data <- readRDS(motif_file)
    motif_mat <- as.matrix(assay(motif_data))
    
    # CRITICAL: Extract and preserve motif names
    motif_names <- rownames(motif_mat)
    cat("Total motifs:", length(motif_names), "\n")
    cat("First few motif names:", paste(head(motif_names, 3), collapse=", "), "\n")
    
    # Ensure columns are in correct order
    cols_exist <- available_celltypes[available_celltypes %in% colnames(motif_mat)]
    
    if(length(cols_exist) > 0) {
      motif_mat <- motif_mat[, cols_exist, drop = FALSE]
      
      # Clean matrix
      keep_rows <- complete.cases(motif_mat) & rowSums(motif_mat != 0) > 0
      motif_mat <- motif_mat[keep_rows, , drop = FALSE]
      motif_names <- rownames(motif_mat)  # Update names after filtering
      
      cat("Motifs after filtering:", length(motif_names), "\n")
      
      if(nrow(motif_mat) > 0) {
        # Limit number of motifs shown
        if(nrow(motif_mat) > 100) {
          motif_var <- apply(motif_mat, 1, var, na.rm = TRUE)
          top_motifs <- names(sort(motif_var, decreasing = TRUE))[1:100]
          motif_mat <- motif_mat[top_motifs, , drop = FALSE]
          motif_names <- rownames(motif_mat)
          cat("Showing top 100 motifs by variance\n")
        }
        
        # Scale but PRESERVE ROW NAMES
        motif_mat_scaled <- t(scale(t(motif_mat)))
        motif_mat_scaled[is.nan(motif_mat_scaled)] <- 0
        rownames(motif_mat_scaled) <- motif_names  # CRITICAL: Restore names
        
        # Calculate dimensions
        n_rows <- nrow(motif_mat_scaled)
        n_cols <- ncol(motif_mat_scaled)
        width <- max(2400, n_cols * 150 + 600)
        height <- max(2800, n_rows * 20 + 600)  # More height per row for names
        
        png(
          filename = file.path(proj_dir, "MotifEnrichment_Heatmap.png"),
          width = width,
          height = height,
          res = 200
        )
        
        # Create heatmap with row names forced to show
        ht <- Heatmap(
          motif_mat_scaled,
          name = "Motif Enrichment\n(z-score)",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          show_row_names = TRUE,  # FORCE show row names
          row_names_gp = gpar(fontsize = max(6, min(10, 200/n_rows))),  # Dynamic font size
          row_names_max_width = unit(15, "cm"),  # Allow long names
          row_names_side = "right",  # Put names on right for better visibility
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 12),
          border = TRUE,
          rect_gp = gpar(col = "white", lwd = 0.5),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 9)
          )
        )
        
        draw(ht,
             padding = unit(c(2, 2, 2, 4), "cm"),
             heatmap_legend_side = "right")
        
        dev.off()
        cat("✓ Motif enrichment heatmap saved with", n_rows, "motifs\n")
        cat("  Row names are displayed (e.g.,", paste(head(motif_names, 3), collapse=", "), "...)\n")
      }
    }
  }, error = function(e) {
    cat("Error in motif enrichment: ", e$message, "\n")
    if(length(dev.list()) > 0) dev.off()
  })
} else {
  cat("Warning: Motif enrichment file not found\n")
}

# ============================================================
# 3) TF ACTIVITY HEATMAP
# ============================================================
cat("\nComputing TF activity heatmap...\n")

tryCatch({
  motif_matrix <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "MotifMatrix"
  )
  
  # Get deviation scores
  if("deviations" %in% names(assays(motif_matrix))) {
    tf_mat <- assay(motif_matrix, "deviations")
  } else {
    tf_mat <- assay(motif_matrix)
  }
  
  # Get cell types
  cell_types <- getCellColData(proj)$CellType
  
  if(length(available_celltypes) > 0) {
    # Calculate mean per cell type
    tf_by_celltype <- matrix(0,
                            nrow = nrow(tf_mat),
                            ncol = length(available_celltypes),
                            dimnames = list(rownames(tf_mat), available_celltypes))
    
    for(ct in available_celltypes) {
      cells_ct <- which(cell_types == ct)
      if(length(cells_ct) > 0) {
        tf_by_celltype[, ct] <- rowMeans(tf_mat[, cells_ct, drop = FALSE], na.rm = TRUE)
      }
    }
    
    # Clean matrix
    tf_by_celltype <- tf_by_celltype[complete.cases(tf_by_celltype), , drop = FALSE]
    tf_by_celltype <- tf_by_celltype[rowSums(tf_by_celltype != 0) > 0, , drop = FALSE]
    
    if(nrow(tf_by_celltype) > 0) {
      # Select top TFs by variance
      tf_var <- apply(tf_by_celltype, 1, var, na.rm = TRUE)
      tf_var <- tf_var[!is.na(tf_var) & tf_var > 0]
      
      if(length(tf_var) > 0) {
        n_top <- min(top_TFs, length(tf_var))
        top_tf_names <- names(sort(tf_var, decreasing = TRUE))[1:n_top]
        tf_mat_final <- tf_by_celltype[top_tf_names, , drop = FALSE]
        
        # Scale
        tf_mat_scaled <- t(scale(t(tf_mat_final)))
        tf_mat_scaled[is.nan(tf_mat_scaled)] <- 0
        
        # Calculate dimensions
        n_rows <- nrow(tf_mat_scaled)
        n_cols <- ncol(tf_mat_scaled)
        width <- max(2400, n_cols * 150 + 600)
        height <- max(2800, n_rows * 25 + 600)
        
        png(
          filename = file.path(proj_dir, "TFactivity_Heatmap.png"),
          width = width,
          height = height,
          res = 200
        )
        
        ht <- Heatmap(
          tf_mat_scaled,
          name = "TF Activity\n(z-score)",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          row_names_gp = gpar(fontsize = 10),
          row_names_side = "right",
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 12),
          border = TRUE,
          rect_gp = gpar(col = "white", lwd = 0.5)
        )
        
        draw(ht,
             padding = unit(c(2, 2, 2, 4), "cm"),
             heatmap_legend_side = "right")
        
        dev.off()
        cat("✓ TF activity heatmap saved with", n_rows, "TFs\n")
      }
    }
  }
}, error = function(e) {
  cat("Error in TF activity: ", e$message, "\n")
  if(length(dev.list()) > 0) dev.off()
})

cat("\n========================================\n")
cat("ALL PLOTS COMPLETED\n")
cat("========================================\n")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(argparse)
  library(ComplexHeatmap)
  library(matrixStats)
  library(Matrix)
  library(GenomicRanges)
  library(grid)
})

# ----------------------------
# Argument parsing
# ----------------------------
parser <- ArgumentParser(description = "ArchR marker peaks + motif enrichment heatmaps")
parser$add_argument('--project_name', required = TRUE)
parser$add_argument('--annotation_col', default = 'CellType')
parser$add_argument('--motif_set', default = 'cisbp')
parser$add_argument('--cutoff_fdr', type='numeric', default = 0.1)
parser$add_argument('--cutoff_log2fc', type='numeric', default = 0.5)
parser$add_argument('--suffix', default = '_analysis')
parser$add_argument('--top_n', type='integer', default=50)
args <- parser$parse_args()

# ----------------------------
# Output folder
# ----------------------------
out_dir <- file.path(args$project_name, paste0("analysis", args$suffix))
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ----------------------------
# Load project
# ----------------------------
message("Loading ArchR project...")
proj <- loadArchRProject(path = args$project_name, force = FALSE)

if (!args$annotation_col %in% names(proj@cellColData)) {
  stop("Annotation column not found: ", args$annotation_col)
}

# ----------------------------
# Peak calling
# ----------------------------
proj <- addGroupCoverages(proj, groupBy=args$annotation_col, force=TRUE)
proj <- addReproduciblePeakSet(proj, groupBy=args$annotation_col, force=TRUE)
proj <- addPeakMatrix(proj, force=TRUE)

# ----------------------------
# Marker peaks
# ----------------------------
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = args$annotation_col,
  bias = c("TSSEnrichment","log10(nFrags)"),
  testMethod = "wilcoxon"
)

# ----------------------------
# Motif annotations
# ----------------------------
proj <- addMotifAnnotations(
  proj,
  motifSet=args$motif_set,
  name="Motif",
  force=TRUE
)

cutoff_expr <- paste0(
  "FDR <= ", args$cutoff_fdr,
  " & Log2FC >= ", args$cutoff_log2fc
)

motifEnrich <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = proj,
  peakAnnotation="Motif",
  cutOff=cutoff_expr
)

# ----------------------------
# SAVE PROJECT BEFORE PLOTTING
# ----------------------------
message("Saving project with all computations BEFORE plotting...")
proj@projectMetadata$markerPeaks <- markerPeaks
proj@projectMetadata$motifEnrichment <- motifEnrich

tryCatch({
  saveArchRProject(
    proj,
    outputDirectory=args$project_name,
    overwrite=TRUE,
    load=FALSE
  )
}, error=function(e) {
  message("Warning: Could not save project: ", e$message)
})

# ----------------------------
# Save RDS objects
# ----------------------------
saveRDS(markerPeaks, file.path(out_dir, "markerPeaks_object.rds"))
saveRDS(motifEnrich, file.path(out_dir, "motifEnrichment.rds"))

# ----------------------------
# Extract peak IDs properly with NULL check
# ----------------------------
peak_gr <- rowRanges(markerPeaks)

# FIX: Check if peak_gr is NULL and handle appropriately
if(is.null(peak_gr)) {
  message("WARNING: rowRanges(markerPeaks) returned NULL. Using fallback method to get peak IDs.")

  # Get PeakMatrix directly and extract peak IDs from row names
  peakSE <- getMatrixFromProject(proj, "PeakMatrix")
  peak_mat <- assay(peakSE)
  peak_rows <- rownames(peak_mat)

  # Create peak_ids from row names (assuming they're in format "chr:start-end")
  peak_ids <- peak_rows

  # For FDR, create a dummy vector since we don't have real FDR values
  peak_fdr <- matrix(1, nrow=length(peak_ids), ncol=1)
  best_fdr <- rep(1, length(peak_ids))
  names(best_fdr) <- peak_ids

} else {
  # Original code works
  peak_ids <- paste0(
    seqnames(peak_gr), ":",
    start(peak_gr), "-",
    end(peak_gr)
  )

  # ----------------------------
  # Get best FDR
  # ----------------------------
  peak_fdr <- assays(markerPeaks)$FDR
  best_fdr <- matrixStats::rowMins(as.matrix(peak_fdr), na.rm=TRUE)
  names(best_fdr) <- peak_ids
}

best_fdr <- sort(best_fdr)

top_n <- min(args$top_n, length(best_fdr))
top_peaks <- names(best_fdr)[1:top_n]

message("Selected ", length(top_peaks), " top peaks")

# ----------------------------
# Get PeakMatrix
# ----------------------------
peakSE <- getMatrixFromProject(proj, "PeakMatrix")
peak_mat <- assay(peakSE)

# Ensure rownames match peak IDs
peak_rows <- rownames(peak_mat)

# Intersect safely
top_peaks <- intersect(top_peaks, peak_rows)

# fallback if intersect small
if(length(top_peaks) < 5){

  message("Few intersect peaks found — fallback to top matrix peaks")

  peak_var <- matrixStats::rowVars(as.matrix(peak_mat))

  peak_var <- sort(peak_var, decreasing=TRUE)

  top_n <- min(args$top_n, length(peak_var))

  top_peaks <- names(peak_var)[1:top_n]
}

# ----------------------------
# Aggregate peaks by group
# ----------------------------
cell_types <- as.character(proj@cellColData[[args$annotation_col]])
names(cell_types) <- rownames(proj@cellColData)

groups <- unique(cell_types)

agg_mat <- matrix(
  0,
  nrow=length(top_peaks),
  ncol=length(groups)
)

rownames(agg_mat) <- top_peaks
colnames(agg_mat) <- groups

for(i in seq_along(groups)){

  cells <- names(cell_types)[cell_types == groups[i]]

  idx <- which(colnames(peak_mat) %in% cells)

  if(length(idx)>0){
    agg_mat[,i] <- Matrix::rowMeans(
      peak_mat[top_peaks, idx, drop=FALSE]
    )
  }
}

# ----------------------------
# Scale
# ----------------------------
agg_scaled <- t(scale(t(agg_mat)))
agg_scaled[is.na(agg_scaled)] <- 0

# ----------------------------
# Peak heatmap - FIXED with tryCatch
# ----------------------------
message("Generating peak heatmap...")
tryCatch({
  png(
    file.path(out_dir,"MarkerPeak_Heatmap.png"),
    width=2000,
    height=1600,
    res=220,
    type="cairo"  # Add cairo type for better compatibility
  )
  
  ht <- Heatmap(
    agg_scaled,
    name="Accessibility",
    cluster_rows=TRUE,
    cluster_columns=TRUE,
    show_row_names=TRUE,
    show_column_names=TRUE,
    row_names_gp=gpar(fontsize=6),
    column_names_gp=gpar(fontsize=10),
    use_raster=FALSE  # Disable rasterization to avoid memory issues
  )
  
  draw(ht)
  dev.off()
  message("Peak heatmap saved successfully")
}, error=function(e) {
  message("Warning: Could not generate peak heatmap: ", e$message)
  if(length(dev.list()) > 0) dev.off()
})

# ----------------------------
# Motif ranking with NULL check
# ----------------------------
motif_fdr <- tryCatch({
  assays(motifEnrich)$FDR
}, error=function(e) {
  message("WARNING: Could not extract FDR from motifEnrich: ", e$message)
  return(NULL)
})

if(is.null(motif_fdr)) {
  message("WARNING: motif_fdr is NULL. Using fallback method with motif matrix variance.")
  
  # Get motif matrix and rank by variance instead
  motifSE <- getMatrixFromProject(proj, "MotifMatrix")
  motif_mat <- assay(motifSE)
  
  motif_var <- matrixStats::rowVars(as.matrix(motif_mat))
  motif_var <- sort(motif_var, decreasing=TRUE)
  
  top_n <- min(args$top_n, length(motif_var))
  top_motifs <- names(motif_var)[1:top_n]
  
} else {
  # Original FDR-based ranking
  best_motif_fdr <- matrixStats::rowMins(
    as.matrix(motif_fdr),
    na.rm=TRUE
  )
  
  best_motif_fdr <- sort(best_motif_fdr)
  
  top_n <- min(args$top_n, length(best_motif_fdr))
  
  top_motifs <- names(best_motif_fdr)[1:top_n]
}

# ----------------------------
# Motif matrix
# ----------------------------
motifSE <- getMatrixFromProject(proj, "MotifMatrix")
motif_mat <- assay(motifSE)

top_motifs <- intersect(top_motifs, rownames(motif_mat))

if(length(top_motifs) < 5){

  message("Few motif matches — using most variable motifs")

  motif_var <- matrixStats::rowVars(as.matrix(motif_mat))

  motif_var <- sort(motif_var, decreasing=TRUE)

  top_n <- min(args$top_n, length(motif_var))

  top_motifs <- names(motif_var)[1:top_n]
}

motif_scaled <- t(scale(t(motif_mat[top_motifs,,drop=FALSE])))
motif_scaled[is.na(motif_scaled)] <- 0

# ----------------------------
# Motif heatmap - FIXED with tryCatch
# ----------------------------
message("Generating motif heatmap...")
tryCatch({
  png(
    file.path(out_dir,"MotifEnrichment_Heatmap.png"),
    width=1800,
    height=1400,
    res=220,
    type="cairo"  # Add cairo type for better compatibility
  )
  
  ht <- Heatmap(
    motif_scaled,
    name="Motif Enrichment",
    cluster_rows=TRUE,
    cluster_columns=TRUE,
    show_row_names=TRUE,
    row_names_gp=gpar(fontsize=6),
    column_names_gp=gpar(fontsize=10),
    use_raster=FALSE  # Disable rasterization to avoid memory issues
  )
  
  draw(ht)
  dev.off()
  message("Motif heatmap saved successfully")
}, error=function(e) {
  message("Warning: Could not generate motif heatmap: ", e$message)
  if(length(dev.list()) > 0) dev.off()
})

message("✓ ANALYSIS COMPLETE")
message("Output folder: ", out_dir)

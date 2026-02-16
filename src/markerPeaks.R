#!/usr/bin/env Rscript

# ----------------------------
# Full ArchR marker peaks pipeline from annotated project
# ----------------------------

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(argparse)
  library(ComplexHeatmap)
})

# ----------------------------
# Argument Parsing
# ----------------------------
parser <- ArgumentParser(description = 'Compute marker peaks for an annotated ArchR project')
parser$add_argument('--project_name', type = 'character', required = TRUE,
                    help = 'Path to annotated ArchR project')
parser$add_argument('--annotation_col', type = 'character', default = 'CellType',
                    help = 'Column in cellColData to use for grouping')
parser$add_argument('--suffix', type = 'character', default = '_peaks',
                    help = 'Suffix for output ArchR project directory')

args <- parser$parse_args()

# ----------------------------
# Load existing ArchR project
# ----------------------------
proj <- loadArchRProject(path = args$project_name)

# ----------------------------
# Step 1: Add group coverages
# ----------------------------
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = args$annotation_col
)

# ----------------------------
# Step 2: Add reproducible peak set
# ----------------------------
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = args$annotation_col
)

# ----------------------------
# Step 3: Add PeakMatrix
# ----------------------------
proj <- addPeakMatrix(proj)

# ----------------------------
# Step 4: Compute marker peaks
# ----------------------------
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = args$annotation_col,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# ----------------------------
# Step 5: Build summary table for downstream usage
# ----------------------------
groups <- colnames(markerPeaks)
peak_summary_list <- list()

for (i in seq_along(groups)) {
  group <- groups[i]
  
  group_data <- assay(markerPeaks)[, i, drop = FALSE]
  
  if (nrow(group_data) > 0) {
    df <- as.data.frame(group_data)
    colnames(df) <- "Value"
    df$Group <- group
    df$Peak <- rownames(df)
    
    if ("Log2FC" %in% names(assays(markerPeaks))) {
      df$Log2FC <- assay(markerPeaks, "Log2FC")[, i]
    }
    if ("FDR" %in% names(assays(markerPeaks))) {
      df$FDR <- assay(markerPeaks, "FDR")[, i]
    }
    if ("Mean" %in% names(assays(markerPeaks))) {
      df$Mean <- assay(markerPeaks, "Mean")[, i]
    }
    
    peak_summary_list[[group]] <- df
  } else {
    message(paste("No marker peaks found for group:", group))
  }
}

if (length(peak_summary_list) > 0) {
  peak_summary <- do.call(rbind, peak_summary_list)
  proj@projectMetadata$markerPeaks <- peak_summary
  message(paste("Found marker peaks for", length(peak_summary_list), "groups"))
} else {
  message("No marker peaks found for any group")
  peak_summary <- data.frame()
}

# ----------------------------
# Step 6: SAVE PLOTS AND DATA
# ----------------------------
output_dir <- paste0(args$project_name, args$suffix)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# SAVE PNG ONLY
if (length(peak_summary_list) > 0) {
  tryCatch({
    png(file.path(output_dir, "MarkerPeaks_Heatmap.png"), width = 1200, height = 1000, res = 150)
    plotMarkerHeatmap(markerPeaks, n = 20, transpose = TRUE)
    dev.off()
    message("Heatmap saved to: ", file.path(output_dir, "MarkerPeaks_Heatmap.png"))
    
    # SAVE THE MARKER PEAKS SUMMARY CSV
    write.csv(peak_summary, file.path(output_dir, "marker_peaks_summary.csv"), row.names = FALSE)
    message("Marker peaks summary saved to: ", file.path(output_dir, "marker_peaks_summary.csv"))
    
    # SAVE MARKER PEAKS OBJECT FOR LATER USE
    saveRDS(markerPeaks, file.path(output_dir, "markerPeaks_object.rds"))
    message("Marker peaks R object saved to: ", file.path(output_dir, "markerPeaks_object.rds"))
    
  }, error = function(e) {
    message("ERROR saving plots/data: ", e$message)
  })
} else {
  message("No marker peaks - skipping plots and CSV")
}

# ----------------------------
# Step 7: Save ArchR project with suffix
# ----------------------------
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = output_dir,
  load = FALSE
)

message("======================================")
message("ArchR project with marker peaks saved to: ", output_dir)
message("CHECK THIS DIRECTORY FOR YOUR FILES:")
message("  - MarkerPeaks_Heatmap.png")
message("  - marker_peaks_summary.csv")
message("  - markerPeaks_object.rds")
message("======================================")

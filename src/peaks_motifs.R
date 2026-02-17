#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(argparse)
  library(ComplexHeatmap)
})

# ----------------------------
# Argument parsing
# ----------------------------
parser <- ArgumentParser(description = "ArchR marker peaks + motif enrichment")
parser$add_argument('--project_name', required = TRUE,
                    help = 'Path to ArchR project directory')
parser$add_argument('--annotation_col', default = 'CellType',
                    help = 'Column in cellColData to group cells')
parser$add_argument('--motif_set', default = 'cisbp',
                    help = 'Motif set to use')
parser$add_argument('--cutoff_fdr', type = 'numeric', default = 0.1,
                    help = 'FDR cutoff for motifs')
parser$add_argument('--cutoff_log2fc', type = 'numeric', default = 0.5,
                    help = 'Log2FC cutoff for motifs')
parser$add_argument('--suffix', default = '_analysis',
                    help = 'Suffix for output folder')

args <- parser$parse_args()

# ----------------------------
# Create output directory
# ----------------------------
out_dir <- paste0(args$project_name, args$suffix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------------------------
# Load project
# ----------------------------
message("Loading ArchR project...")
proj <- loadArchRProject(path = args$project_name, force = FALSE)

# ----------------------------
# Validate annotation column
# ----------------------------
if (!args$annotation_col %in% names(proj@cellColData)) {
  stop("ERROR: annotation column not found: ", args$annotation_col)
}

# ----------------------------
# FIX FOR PEAKSET - Load from disk if NULL
# ----------------------------
if (is.null(proj@peakSet)) {
  message("peakSet is NULL - loading from disk...")
  peak_files <- list.files(
    file.path(args$project_name, "PeakCalls", args$annotation_col),
    pattern = "reproduciblePeaks.gr.rds",
    full.names = TRUE
  )
  if (length(peak_files) > 0) {
    proj@peakSet <- readRDS(peak_files[1])
    message("✓ Loaded peakSet with ", length(proj@peakSet), " peaks")
  } else {
    stop("ERROR: No reproducible peaks found in PeakCalls directory")
  }
}

# ----------------------------
# Add group coverages
# ----------------------------
message("Adding group coverages...")
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = args$annotation_col,
  force = TRUE
)

# ----------------------------
# Add reproducible peak set
# ----------------------------
message("Adding reproducible peak set...")
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = args$annotation_col,
  force = TRUE
)

# ----------------------------
# Add PeakMatrix
# ----------------------------
message("Adding PeakMatrix...")
proj <- addPeakMatrix(proj)

# ----------------------------
# Get marker peaks
# ----------------------------
message("Computing marker peaks...")
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = args$annotation_col,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# ----------------------------
# Save marker peaks heatmap
# ----------------------------
message("Generating marker peaks heatmap...")
png(file.path(out_dir, "MarkerPeaks_Heatmap.png"),
    width = 1800, height = 1400, res = 200)
hm_peaks <- plotMarkerHeatmap(
  seMarker = markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  nLabel = 20,
  labelRows = TRUE
)
draw(hm_peaks)
dev.off()

# ----------------------------
# Add motif annotations
# ----------------------------
message("Adding motif annotations...")
proj <- addMotifAnnotations(
  ArchRProj = proj,
  motifSet = args$motif_set,
  name = "Motif",
  force = TRUE
)

# ----------------------------
# Motif enrichment
# ----------------------------
message("Performing motif enrichment...")
cutoff_expr <- paste0("FDR <= ", args$cutoff_fdr,
                      " & Log2FC >= ", args$cutoff_log2fc)

motifEnrich <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = cutoff_expr
)

# ----------------------------
# Save motif enrichment heatmap
# ----------------------------
message("Generating motif enrichment heatmap...")
png(file.path(out_dir, "MotifEnrichment_Heatmap.png"),
    width = 1600, height = 1200, res = 200)
plotEnrichHeatmap(motifEnrich, transpose = TRUE)
dev.off()

# ----------------------------
# Save RDS objects
# ----------------------------
message("Saving RDS objects...")
saveRDS(markerPeaks, file.path(out_dir, "markerPeaks_object.rds"))
saveRDS(motifEnrich, file.path(out_dir, "motifEnrichment.rds"))

# ----------------------------
# Store results in project metadata
# ----------------------------
proj@projectMetadata$markerPeaks <- markerPeaks
proj@projectMetadata$motifEnrichment <- motifEnrich

# ----------------------------
# Save ArchR project
# ----------------------------
message("Saving ArchR project...")
tryCatch({
  saveArchRProject(
    ArchRProj = proj,
    outputDirectory = out_dir,
    load = FALSE,
    overwrite = TRUE
  )
  message("✓ Project saved")
}, error = function(e) {
  message("⚠ Project save failed: ", e$message)
  message("✓ RDS objects and heatmaps saved")
})

# ----------------------------
# Final summary
# ----------------------------
message("\n======================================")
message("ANALYSIS COMPLETE")
message("Output: ", out_dir)
message("\nFiles:")
message("  ✓ MarkerPeaks_Heatmap.png")
message("  ✓ MotifEnrichment_Heatmap.png")
message("  ✓ markerPeaks_object.rds")
message("  ✓ motifEnrichment.rds")
message("======================================")

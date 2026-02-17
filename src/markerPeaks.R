#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(argparse)
  library(ComplexHeatmap)
})

# ----------------------------
# Arguments
# ----------------------------
parser <- ArgumentParser(description = "ArchR marker peaks + heatmap (SAFE)")
parser$add_argument("--project_name", required = TRUE,
                    help = "Path to annotated ArchR project")
parser$add_argument("--annotation_col", default = "CellType",
                    help = "Grouping column")
parser$add_argument("--suffix", default = "_peaks",
                    help = "Suffix for output project")

args <- parser$parse_args()

# ----------------------------
# Load project
# ----------------------------
proj <- loadArchRProject(
  path = args$project_name,
  force = FALSE
)

# ----------------------------
# Add group coverages (idempotent)
# ----------------------------
if (!args$annotation_col %in% names(proj@cellColData)) {
  stop("ERROR: annotation column not found: ", args$annotation_col)
}

proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = args$annotation_col,
  force = TRUE
)

# ----------------------------
# Add reproducible peak set (idempotent)
# ----------------------------
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = args$annotation_col,
  force = TRUE
)

# ----------------------------
# Add PeakMatrix (idempotent)
# ----------------------------
proj <- addPeakMatrix(proj)

# ----------------------------
# Compute marker peaks (THIS IS THE REAL OBJECT)
# ----------------------------
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = args$annotation_col,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# ----------------------------
# Store marker peaks INSIDE ArchR project
# ----------------------------
proj@projectMetadata$markerPeaks <- markerPeaks

# ----------------------------
# Output directory
# ----------------------------
output_dir <- paste0(args$project_name, args$suffix)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ----------------------------
# SAVE HEATMAP (ABSOLUTELY SAFE)
# ----------------------------
png(
  filename = file.path(output_dir, "MarkerPeaks_Heatmap.png"),
  width = 1800,
  height = 1400,
  res = 200
)

hm <- plotMarkerHeatmap(
  seMarker = markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(hm)
dev.off()

# ----------------------------
# Save marker peaks object (optional but useful)
# ----------------------------
saveRDS(
  markerPeaks,
  file.path(output_dir, "markerPeaks_object.rds")
)

# ----------------------------
# Save ArchR project COPY (THIS IS CRITICAL)
# ----------------------------
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = output_dir,
  load = FALSE
)

# ----------------------------
# Final messages
# ----------------------------
message("======================================")
message("MARKER PEAKS PIPELINE COMPLETE")
message("Saved ArchR project to:")
message("  ", output_dir)
message("Files created:")
message("  - MarkerPeaks_Heatmap.png")
message("  - markerPeaks_object.rds")
message("Marker peaks are stored INSIDE ArchR:")
message("  proj@projectMetadata$markerPeaks")
message("======================================")


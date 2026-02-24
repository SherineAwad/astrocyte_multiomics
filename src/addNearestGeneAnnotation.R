#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ArchR)
  library(argparse)
  library(GenomicRanges)
})

# -----------------------------
# Arguments
# -----------------------------
parser <- ArgumentParser()
parser$add_argument("--project_name", required = TRUE, help="Path to existing ArchR project")
args <- parser$parse_args()
proj_dir <- args$project_name

cat("ArchR project directory:", proj_dir, "\n\n")

# -----------------------------
# Load ArchR project
# -----------------------------
cat("Loading ArchR project...\n")
proj <- loadArchRProject(proj_dir)
cat("Project loaded! Cells:", nCells(proj), "\n\n")

# -----------------------------
# Get peaks and genes
# -----------------------------
cat("Getting peaks and gene annotations...\n")
peaks <- getPeakSet(proj)               # GRanges of peaks
genes_list <- getGeneAnnotation(proj)   # This may be GRangesList
cat("Number of peaks:", length(peaks), "\n")
cat("Number of gene sets:", length(genes_list), "\n\n")

# -----------------------------
# Flatten gene annotation
# -----------------------------
cat("Flattening gene annotations...\n")
genes <- unlist(genes_list, use.names = FALSE)

# Ensure gene symbol column exists
if(!"symbol" %in% colnames(mcols(genes))){
  stop("Gene GRanges does not have a 'symbol' column. Check metadata names with colnames(mcols(genes))")
}

cat("Number of genes after flattening:", length(genes), "\n\n")

# -----------------------------
# Compute nearest gene for each peak
# -----------------------------
cat("Computing nearest gene for each peak...\n")
nearest_idx <- nearest(peaks, genes)

peak2gene <- data.frame(
  peak = paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks)),
  nearest_gene = genes$symbol[nearest_idx]
)
cat("Nearest-gene annotation complete! Total peaks annotated:", nrow(peak2gene), "\n\n")

# -----------------------------
# Save as CSV
# -----------------------------
out_csv <- file.path(proj_dir, "peak_nearest_gene_annotation.csv")
write.csv(peak2gene, out_csv, row.names = FALSE)
cat("Peak → nearest gene table saved to:", out_csv, "\n\n")

# -----------------------------
# Save ArchR project (unchanged)
# -----------------------------
cat("Saving ArchR project back to the same folder...\n")
saveArchRProject(proj, outputDirectory = proj_dir, load = FALSE)
cat("Project saved ✅\n")

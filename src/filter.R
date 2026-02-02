library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = 'ArchR single-cell ATAC-seq processing with filtering')

# Required argument
parser$add_argument('--project_name', type = 'character', required = TRUE,
                    help = 'Name of the ArchR project in current directory')

# Optional arguments with defaults
parser$add_argument('--genome', type = 'character', default = 'mm10',
                    help = 'Genome assembly (default: mm10)')
parser$add_argument('--threads', type = 'integer', default = 4,
                    help = 'Number of threads to use (default: 4)')
parser$add_argument('--suffix', type = 'character', default = '_filtered',
                    help = 'Suffix to add to project name for filtered output (default: _filtered)')

# Filtering thresholds with defaults
parser$add_argument('--min_tss_enrichment', type = 'double', default = 10,
                    help = 'Minimum TSS enrichment score (default: 10)')
parser$add_argument('--min_nfrags', type = 'integer', default = 5000,
                    help = 'Minimum number of fragments (default: 5000)')
parser$add_argument('--min_gex_ngenes', type = 'integer', default = 1000,
                    help = 'Minimum number of genes in RNA-seq (default: 1000)')
parser$add_argument('--max_gex_ngenes', type = 'integer', default = 7000,
                    help = 'Maximum number of genes in RNA-seq (default: 7000)')
parser$add_argument('--min_gex_numi', type = 'integer', default = 1500,
                    help = 'Minimum number of UMIs in RNA-seq (default: 1500)')
parser$add_argument('--max_gex_numi', type = 'integer', default = 30000,
                    help = 'Maximum number of UMIs in RNA-seq (default: 30000)')

# Parse arguments
args <- parser$parse_args()

# Initialize ArchR
addArchRThreads(threads = args$threads)
addArchRGenome(args$genome)

# Define input and output paths
input_dir <- args$project_name
output_dir <- paste0(args$project_name, args$suffix)

# Load ArchR project from original directory
cat(sprintf("Loading ArchR project from: %s\n", input_dir))
proj_ALL <- loadArchRProject(path = input_dir, force = FALSE, showLogo = TRUE)

# Display initial cell counts
cat("\n=== Initial Cell Counts ===\n")
cat(sprintf("Total cells before filtering: %d\n", nCells(proj_ALL)))
cat("Cells per sample:\n")
print(table(proj_ALL$Sample))

# Create PRE-FILTERING QC plots
cat("\n=== Generating Pre-filtering QC Plots ===\n")
prefilter_figure_name <- paste0(input_dir, "_preFilterQC.pdf")
pdf(file = prefilter_figure_name, width = 12, height = 8)

# Create each pre-filter plot separately and print them in the PDF
# Remove ALL NA values before plotting to avoid quantile error
proj_ALL_clean <- proj_ALL[
    !is.na(proj_ALL$TSSEnrichment) & 
    !is.na(proj_ALL$nFrags) & 
    !is.na(proj_ALL$Gex_nUMI) & 
    !is.na(proj_ALL$Gex_nGenes)
]

pre_p1 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(pre_p1)

pre_p2 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(pre_p2)

pre_p3 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nUMI",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(pre_p3)

pre_p4 <- plotGroups(
    ArchRProj = proj_ALL_clean,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nGenes",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(pre_p4)

dev.off()

cat(sprintf("Pre-filtering QC plots saved to: %s\n", prefilter_figure_name))

# Apply filtering (keep original code exactly)
cat(sprintf("\n=== Applying Filters ===\n"))
cat(sprintf("TSSEnrichment > %.2f\n", args$min_tss_enrichment))
cat(sprintf("nFrags > %d\n", args$min_nfrags))
cat(sprintf("Gex_nGenes: %d - %d\n", args$min_gex_ngenes, args$max_gex_ngenes))
cat(sprintf("Gex_nUMI: %d - %d\n", args$min_gex_numi, args$max_gex_numi))

proj_ALL <- proj_ALL[
    proj_ALL$TSSEnrichment > args$min_tss_enrichment &
    proj_ALL$nFrags > args$min_nfrags &
    !is.na(proj_ALL$Gex_nUMI)
]

proj_ALL <- proj_ALL[
    proj_ALL$Gex_nGenes > args$min_gex_ngenes &
    proj_ALL$Gex_nGenes < args$max_gex_ngenes &
    proj_ALL$Gex_nUMI > args$min_gex_numi &
    proj_ALL$Gex_nUMI < args$max_gex_numi
]

# Display post-filtering cell counts
cat("\n=== Post-filtering Cell Counts ===\n")
cat(sprintf("Total cells after filtering: %d\n", nCells(proj_ALL)))
cat("Cells per sample:\n")
print(table(proj_ALL$Sample))

# Create POST-FILTERING QC plots
cat("\n=== Generating Post-filtering QC Plots ===\n")
figure_name <- paste0(input_dir, "_postFilterQC.pdf")

pdf(file = figure_name, width = 12, height = 8)

# Create each plot separately and print them in the PDF
p1 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(p1)

p2 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(p2)

p3 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nUMI",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(p3)

p4 <- plotGroups(
    ArchRProj = proj_ALL,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "Gex_nGenes",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
print(p4)

dev.off()

cat(sprintf("Post-filtering QC plots saved to: %s\n", figure_name))

# Save filtered project to new directory
cat(sprintf("\n=== Saving filtered project to: %s ===\n", output_dir))
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = output_dir, load = FALSE)

cat("\n=== Script completed successfully ===\n")
cat(sprintf("Original project remains at: %s\n", input_dir))
cat(sprintf("Filtered project saved to: %s\n", output_dir))
cat(sprintf("Pre-filtering plots saved to: %s\n", prefilter_figure_name))
cat(sprintf("Post-filtering plots saved to: %s\n", figure_name))

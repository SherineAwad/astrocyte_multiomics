# Astrocyte multiomics data: preprocessing for scenic+ 


### Filtering Thresholds

We used the following filtering criteria: 

* **Minimum TSS enrichment**: Ensures sufficient transcription start site (TSS) signal for ATAC-seq data. Default is **10**.
* **Minimum number of fragments**: Filters out cells with too few ATAC-seq fragments. Default is **5000**.
* **Minimum number of genes detected in RNA-seq**: Excludes low-quality cells with sparse transcript coverage. Default is **1000**.
* **Maximum number of genes detected in RNA-seq**: Removes potential multiplets or highly overexpressed cells. Default is **7000**.
* **Minimum number of UMIs in RNA-seq**: Ensures sufficient sequencing depth per cell. Default is **1500**.
* **Maximum number of UMIs in RNA-seq**: Filters out cells with abnormally high counts that may indicate doublets. Default is **30000**.


[pre filter](Astrocyte_preFilterQC.pdf)

[post filter](Astrocyte_filtered_postFilterQC.pdf)


## UMAPs

![](Astrocyte_filtered_UMAP_Combined_UMAP.png)
![](Astrocyte_filtered_UMAP_RNA_UMAP.png)
![](Astrocyte_filtered_UMAP_ATAC_UMAP.png)



<!-- QC plots -->
<img src="figures/Astrocyte_filtered_UMAP_Clusters_RNA_QC_plots.png" alt="RNA QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_ATAC_QC_plots.png" alt="ATAC QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_Combined_QC_plots.png" alt="Combined QC" width="33%">

<!-- nFrags QC -->
<img src="figures/Astrocyte_filtered_UMAP_Clusters_RNA_nFrags_QC.png" alt="RNA nFrags QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_ATAC_nFrags_QC.png" alt="ATAC nFrags QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_Combined_nFrags_QC.png" alt="Combined nFrags QC" width="33%">

<!-- TSSEnrichment QC -->
<img src="figures/Astrocyte_filtered_UMAP_Clusters_RNA_TSSEnrichment_QC.png" alt="RNA TSS QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_ATAC_TSSEnrichment_QC.png" alt="ATAC TSS QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_Combined_TSSEnrichment_QC.png" alt="Combined TSS QC" width="33%">

<!-- NucleosomeRatio QC -->
<img src="figures/Astrocyte_filtered_UMAP_Clusters_RNA_NucleosomeRatio_QC.png" alt="RNA Nuc QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_ATAC_NucleosomeRatio_QC.png" alt="ATAC Nuc QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_Combined_NucleosomeRatio_QC.png" alt="Combined Nuc QC" width="33%">

<!-- BlacklistRatio QC -->
<img src="figures/Astrocyte_filtered_UMAP_Clusters_RNA_BlacklistRatio_QC.png" alt="RNA Blacklist QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_ATAC_BlacklistRatio_QC.png" alt="ATAC Blacklist QC" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Clusters_Combined_BlacklistRatio_QC.png" alt="Combined Blacklist QC" width="33%">

<!-- UMAP per sample -->
<img src="figures/Astrocyte_filtered_UMAP_RNA_UMAP_per_sample.png" alt="RNA UMAP per sample" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_ATAC_UMAP_per_sample.png" alt="ATAC UMAP per sample" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Combined_UMAP_per_sample.png" alt="Combined UMAP per sample" width="33%">

<!-- UMAP per sample grid -->
<img src="figures/Astrocyte_filtered_UMAP_RNA_UMAP_per_sample_grid.png" alt="RNA UMAP grid" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_ATAC_UMAP_per_sample_grid.png" alt="ATAC UMAP grid" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Combined_UMAP_per_sample_grid.png" alt="Combined UMAP grid" width="33%">

<!-- UMAP -->
<img src="figures/Astrocyte_filtered_UMAP_RNA_UMAP.png" alt="RNA UMAP" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_ATAC_UMAP.png" alt="ATAC UMAP" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Combined_UMAP.png" alt="Combined UMAP" width="33%">




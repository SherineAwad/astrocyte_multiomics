# Astrocyte multiomics data: preprocessing for scenic+ 


### Samples 

- Control,13005-TH-1/
- KO: 13005-TH-2/


### Filtering Thresholds

We used the following filtering criteria: 

* **Minimum TSS enrichment**: Ensures sufficient transcription start site (TSS) signal for ATAC-seq data. Default is **10**.
* **Minimum number of fragments**: Filters out cells with too few ATAC-seq fragments. Default is **5000**.
* **Minimum number of genes detected in RNA-seq**: Excludes low-quality cells with sparse transcript coverage. Default is **1000**.
* **Maximum number of genes detected in RNA-seq**: Removes potential multiplets or highly overexpressed cells. Default is **7000**.
* **Minimum number of UMIs in RNA-seq**: Ensures sufficient sequencing depth per cell. Default is **1500**.
* **Maximum number of UMIs in RNA-seq**: Filters out cells with abnormally high counts that may indicate doublets. Default is **30000**.


[pre filter](figures/Astrocyte_preFilterQC.pdf)

[post filter](figures/Astrocyte_postFilterQC.pdf)


## UMAPs

#### RNA 
![](figures/Astrocyte_filtered_UMAP_RNA_UMAP.png?v=3)
#### ATAC
![](figures/Astrocyte_filtered_UMAP_ATAC_UMAP.png?v=3)
#### Combined 
![](figures/Astrocyte_filtered_UMAP_Combined_UMAP.png?v=3)

## QCs 

#### RNA 
![](figures/Astrocyte_filtered_UMAP_Clusters_RNA_QC.png?v=3)
#### ATAC 
![](figures/Astrocyte_filtered_UMAP_Clusters_ATAC_QC.png?v=3)
#### Combined 
![](figures/Astrocyte_filtered_UMAP_Clusters_Combined_QC.png?v=3)


## Per sample UMAP 

#### RNA
![](figures/Astrocyte_filtered_UMAP_RNA_SAMPLE_UMAP.png?v=3)
#### ATAC
![](figures/Astrocyte_filtered_UMAP_ATAC_SAMPLE_UMAP.png?v=3)
#### Combined
![](figures/Astrocyte_filtered_UMAP_Combined_SAMPLE_UMAP.png?v=3)


## Clusters 

![](figures/Astrocyte_filtered_UMAP_Combined_SAMPLES_GRID.png?v=3)

#### Combined 

![](figures/Astrocyte_filtered_UMAP_BothSamples_CombinedClusters_UMAP.png?v=3)






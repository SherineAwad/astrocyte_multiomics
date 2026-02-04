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
![](figures/Astrocyte_filtered_UMAP_RNA_UMAP.png?v=4)
#### ATAC
![](figures/Astrocyte_filtered_UMAP_ATAC_UMAP.png?v=4)
#### Combined 
![](figures/Astrocyte_filtered_UMAP_Combined_UMAP.png?v=4)

## QCs 

#### RNA 
![](figures/Astrocyte_filtered_UMAP_Clusters_RNA_QC.png?v=4)
#### ATAC 
![](figures/Astrocyte_filtered_UMAP_Clusters_ATAC_QC.png?v=4)
#### Combined 
![](figures/Astrocyte_filtered_UMAP_Clusters_Combined_QC.png?v=4)


## Per sample UMAP 

#### RNA
![](figures/Astrocyte_filtered_UMAP_RNA_SAMPLE_UMAP.png?v=4)
#### ATAC
![](figures/Astrocyte_filtered_UMAP_ATAC_SAMPLE_UMAP.png?v=4)
#### Combined
![](figures/Astrocyte_filtered_UMAP_Combined_SAMPLE_UMAP.png?v=4)


## Clusters 

![](figures/Astrocyte_filtered_UMAP_Combined_SAMPLES_GRID.png?v=4)

#### Combined 

![](figures/Astrocyte_filtered_UMAP_BothSamples_CombinedClusters_UMAP.png?v=4)

## Marker genes 

<img src="figures/Astrocyte_filtered_UMAP_Aldh1l1_UMAP.png?v=2" alt="aldh1l1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Aldoc_UMAP.png?v=2" alt="aldoc" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Apoe_UMAP.png?v=2" alt="apoe" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Aqp4_UMAP.png?v=2" alt="aqp4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Ascl1_UMAP.png?v=2" alt="ascl1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Bcl11a_UMAP.png?v=2" alt="bcl11a" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Cdk1_UMAP.png?v=2" alt="cdk1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Clu_UMAP.png?v=2" alt="clu" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Dcx_UMAP.png?v=2" alt="dcx" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Elavl3_UMAP.png?v=2" alt="elavl3" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Elavl4_UMAP.png?v=2" alt="elavl4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Glul_UMAP.png?v=2" alt="glul" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Hes5_UMAP.png?v=2" alt="hes5" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Malat1_UMAP.png?v=2" alt="malat1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Mbp_UMAP.png?v=2" alt="mbp" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Mki67_UMAP.png?v=2" alt="mki67" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Mpeg1_UMAP.png?v=2" alt="mpeg1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Neurog2_UMAP.png?v=2" alt="neurog2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Npy_UMAP.png?v=2" alt="npy" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Ntsr2_UMAP.png?v=2" alt="ntsr2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Olig1_UMAP.png?v=2" alt="olig1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Olig2_UMAP.png?v=2" alt="olig2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pdgfra_UMAP.png?v=2" alt="pdgfra" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Plp1_UMAP.png?v=2" alt="plp1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Prox1_UMAP.png?v=2" alt="prox1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_S100b_UMAP.png?v=2" alt="s100b" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Scrt1_UMAP.png?v=2" alt="scrt1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Slc1a3_UMAP.png?v=2" alt="slc1a3" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox10_UMAP.png?v=2" alt="sox10" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox2_UMAP.png?v=2" alt="sox2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Sox9_UMAP.png?v=2" alt="sox9" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Tcf7l2_UMAP.png?v=2" alt="tcf7l2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Tubb3_UMAP.png?v=2" alt="tubb3" width="33%">






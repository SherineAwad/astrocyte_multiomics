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

<img src="figures/Astrocyte_filtered_UMAP_Malat1_UMAP.png?v=3" alt="malat1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Slc1a3_UMAP.png?v=3" alt="slc1a3" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox2_UMAP.png?v=3" alt="sox2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Notch1_UMAP.png?v=3" alt="notch1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Glul_UMAP.png?v=3" alt="glul" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Elavl4_UMAP.png?v=3" alt="elavl4" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Prdx6_UMAP.png?v=3" alt="prdx6" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox11_UMAP.png?v=3" alt="sox11" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Lhx2_UMAP.png?v=3" alt="lhx2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Apoe_UMAP.png?v=3" alt="apoe" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Ascl1_UMAP.png?v=3" alt="ascl1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Olig2_UMAP.png?v=3" alt="olig2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Slc6a9_UMAP.png?v=3" alt="slc6a9" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Aldh1l1_UMAP.png?v=3" alt="aldh1l1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pax6_UMAP.png?v=3" alt="pax6" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Ntsr2_UMAP.png?v=3" alt="ntsr2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Aqp4_UMAP.png?v=3" alt="aqp4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_S100b_UMAP.png?v=3" alt="s100b" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Aldoc_UMAP.png?v=3" alt="aldoc" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox10_UMAP.png?v=3" alt="sox10" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Dcx_UMAP.png?v=3" alt="dcx" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Bcl11a_UMAP.png?v=3" alt="bcl11a" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Bsn_UMAP.png?v=3" alt="bsn" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Tubb3_UMAP.png?v=3" alt="tubb3" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Vim_UMAP.png?v=3" alt="vim" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Hes1_UMAP.png?v=3" alt="hes1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Clu_UMAP.png?v=3" alt="clu" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Gad1_UMAP.png?v=3" alt="gad1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Mbp_UMAP.png?v=3" alt="mbp" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Hes5_UMAP.png?v=3" alt="hes5" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Prox1_UMAP.png?v=3" alt="prox1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Plp1_UMAP.png?v=3" alt="plp1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Rpe65_UMAP.png?v=3" alt="rpe65" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Mki67_UMAP.png?v=3" alt="mki67" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Insm1_UMAP.png?v=3" alt="insm1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Rlbp1_UMAP.png?v=3" alt="rlbp1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Cdk1_UMAP.png?v=3" alt="cdk1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Neurog2_UMAP.png?v=3" alt="neurog2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Scrt1_UMAP.png?v=3" alt="scrt1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Gfap_UMAP.png?v=3" alt="gfap" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Rbfox3_UMAP.png?v=3" alt="rbfox3" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Abca8a_UMAP.png?v=3" alt="abca8a" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Mpeg1_UMAP.png?v=3" alt="mpeg1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Csf1r_UMAP.png?v=3" alt="csf1r" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pdgfra_UMAP.png?v=3" alt="pdgfra" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Atoh7_UMAP.png?v=3" alt="atoh7" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Calb1_UMAP.png?v=3" alt="calb1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Calb2_UMAP.png?v=3" alt="calb2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Npy_UMAP.png?v=3" alt="npy" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Slc17a7_UMAP.png?v=3" alt="slc17a7" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Prdm1_UMAP.png?v=3" alt="prdm1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Tfap2a_UMAP.png?v=3" alt="tfap2a" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pax2_UMAP.png?v=3" alt="pax2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Cabp5_UMAP.png?v=3" alt="cabp5" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Kcnj8_UMAP.png?v=3" alt="kcnj8" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Chat_UMAP.png?v=3" alt="chat" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Ccr2_UMAP.png?v=3" alt="ccr2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Foxn4_UMAP.png?v=3" alt="foxn4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Nrl_UMAP.png?v=3" alt="nrl" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Emx1_UMAP.png?v=3" alt="emx1" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Rho_UMAP.png?v=3" alt="rho" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Lhx4_UMAP.png?v=3" alt="lhx4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Otx2_UMAP.png?v=3" alt="otx2" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Lhx1_UMAP.png?v=3" alt="lhx1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Isl1_UMAP.png?v=3" alt="isl1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Arr3_UMAP.png?v=3" alt="arr3" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Acta2_UMAP.png?v=3" alt="acta2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Tie1_UMAP.png?v=3" alt="tie1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sebox_UMAP.png?v=3" alt="sebox" width="33%">

<img src="figures/Astrocyte_filtered_UMAP_Pou4f2_UMAP.png?v=3" alt="pou4f2" width="33%">







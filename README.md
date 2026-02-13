# Astrocyte multiomics data: preprocessing for scenic+ 


### Samples 

| sample_name | atac_file                             | rna_file                                 |
| ----------- | ------------------------------------- | ---------------------------------------- |
| KO1         | 13005-TH-2/outs/atac_fragments.tsv.gz | 13005-TH-2/outs/raw_feature_bc_matrix.h5 |
| Control     | 13784-TH-1/outs/atac_fragments.tsv.gz | 13784-TH-1/outs/raw_feature_bc_matrix.h5 |
| KO2         | 13784-TH-2/outs/atac_fragments.tsv.gz | 13784-TH-2/outs/raw_feature_bc_matrix.h5 |


### Filtering Thresholds

We used the following filtering criteria: 

### Cell Filtering Criteria

Cells were retained if they met the following thresholds: TSS enrichment ≥ 10, number of ATAC fragments ≥ 1,000, number of detected genes between 1,000 and 7,000, and number of RNA UMIs between 1,500 and 30,000. Only cells passing all criteria were included in downstream analyses.

# Filtering Thresholds

* min_tss_enrichment: 10
* min_nfrags: 1000
* min_gex_ngenes: 1000
* max_gex_ngenes: 7000
* min_gex_numi: 1500
* max_gex_numi: 30000



### Per sample filtering

| Sample  | Before | After | % Kept |
|---------|--------|-------|--------|
| KO2     | 59,364 | 26,915 | 45.3% |
| Control | 60,328 | 27,778 | 46.0% |
| KO1     | 17,295 | 7,529  | 43.5% |

---


[pre filter](figures/Astrocyte_preFilterQC.pdf)

[post filter](figures/Astrocyte_postFilterQC.pdf)


## UMAPs

#### RNA 
![](figures/Astrocyte_filtered_UMAP_RNA_UMAP.png?v=7)
#### ATAC
![](figures/Astrocyte_filtered_UMAP_ATAC_UMAP.png?v=7)
#### Combined 
![](figures/Astrocyte_filtered_UMAP_Combined_UMAP.png?v=7)

## QCs 

#### RNA
![](figures/Astrocyte_filtered_UMAP_Clusters_RNA_QC.png?v=7)
#### ATAC 
![](figures/Astrocyte_filtered_UMAP_Clusters_ATAC_QC.png?v=7)
#### Combined 
![](figures/Astrocyte_filtered_UMAP_Clusters_Combined_QC.png?v=7)


## Per sample UMAP 

#### RNA
![](figures/Astrocyte_filtered_UMAP_RNA_SAMPLE_UMAP.png?v=7)
#### ATAC
![](figures/Astrocyte_filtered_UMAP_ATAC_SAMPLE_UMAP.png?v=7)
#### Combined
![](figures/Astrocyte_filtered_UMAP_Combined_SAMPLE_UMAP.png?v=7)


## Clusters 

![](figures/Astrocyte_filtered_UMAP_Combined_SAMPLES_GRID.png?v=7)

#### Combined 

![](figures/Astrocyte_filtered_UMAP_BothSamples_CombinedClusters_UMAP.png?v=7)

## Marker genes 

<img src="figures/Astrocyte_filtered_UMAP_Aldh1l1_UMAP.png?v=5" alt="aldh1l1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Hes5_UMAP.png?v=5" alt="hes5" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Prox1_UMAP.png?v=5" alt="prox1" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Aldoc_UMAP.png?v=5" alt="aldoc" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Clu_UMAP.png?v=5" alt="clu" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Lhx8_UMAP.png?v=5" alt="lhx8" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Reln_UMAP.png?v=5" alt="reln" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Apoe_UMAP.png?v=5" alt="apoe" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Malat1_UMAP.png?v=5" alt="malat1" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Aqp4_UMAP.png?v=5" alt="aqp4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Mbp_UMAP.png?v=5" alt="mbp" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pou3f2_UMAP.png?v=5" alt="pou3f2" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Arx_UMAP.png?v=5" alt="arx" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Mpeg1_UMAP.png?v=5" alt="mpeg1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_S100b_UMAP.png?v=5" alt="s100b" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Ascl1_UMAP.png?v=5" alt="ascl1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Dcx_UMAP.png?v=5" alt="dcx" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Necab1_UMAP.png?v=5" alt="necab1" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Scrt1_UMAP.png?v=5" alt="scrt1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Dlx1_UMAP.png?v=5" alt="dlx1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Necab2_UMAP.png?v=5" alt="necab2" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Scrt2_UMAP.png?v=5" alt="scrt2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Tubb3_UMAP.png?v=5" alt="tubb3" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Dlx2_UMAP.png?v=5" alt="dlx2" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Neurog2_UMAP.png?v=5" alt="neurog2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Slc17a6_UMAP.png?v=5" alt="slc17a6" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Bcl11a_UMAP.png?v=5" alt="bcl11a" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Dlx5_UMAP.png?v=5" alt="dlx5" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Npy_UMAP.png?v=5" alt="npy" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Slc1a3_UMAP.png?v=5" alt="slc1a3" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Bcl11b_UMAP.png?v=5" alt="bcl11b" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Dlx6_UMAP.png?v=5" alt="dlx6" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Ntsr2_UMAP.png?v=5" alt="nts2r2" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Slc32a1_UMAP.png?v=5" alt="slc32a1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Elavl3_UMAP.png?v=5" alt="elavl3" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Olig1_UMAP.png?v=5" alt="olig1" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Sox10_UMAP.png?v=5" alt="sox10" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Calb1_UMAP.png?v=5" alt="calb1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Elavl4_UMAP.png?v=5" alt="elavl4" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Olig2_UMAP.png?v=5" alt="olig2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox2_UMAP.png?v=5" alt="sox2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Calb2_UMAP.png?v=5" alt="calb2" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Foxp2_UMAP.png?v=5" alt="foxp2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pcp4_UMAP.png?v=5" alt="pcp4" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Sox9_UMAP.png?v=5" alt="sox9" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Cdk1_UMAP.png?v=5" alt="cdk1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Gad1_UMAP.png?v=5" alt="gad1" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Pdgfra_UMAP.png?v=5" alt="pdgfra" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Sst_UMAP.png?v=5" alt="sst" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Gad2_UMAP.png?v=5" alt="gad2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Plp1_UMAP.png?v=5" alt="plp1" width="33%">
<img src="figures/Astrocyte_filtered_UMAP_Tcf7l2_UMAP.png?v=5" alt="tcf7l2" width="33%"><img src="figures/Astrocyte_filtered_UMAP_Glul_UMAP.png?v=5" alt="glul" width="33%">




## Annotations: NEED MORE WORK HERE 







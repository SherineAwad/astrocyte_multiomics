# Snakefile
configfile: "config.yaml"

rule all:
    input:
        directory(config["project_name"]),
        directory(config["project_name"] + config["filter_suffix"]),
        directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"]),
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_Combined_UMAP.pdf",
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_ATAC_UMAP.pdf",
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_RNA_UMAP.pdf",
        "markers_done.txt",
        directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"]) + config["annotation_suffix"],
        directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + config["annotation_suffix"]+config['peaks_suffix']),
        #directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + config["annotation_suffix"]+config['peaks_suffix']+config['motif_suffix'])
rule preprocess:
    input:
        config_csv=config["samples_csv"]
    output:
        directory(config["project_name"])
    params:
        project=config["project_name"],
        genome=config["genome"],
        threads=config["threads"]
    shell:
        "Rscript src/preprocess.R {input.config_csv} {params.project} {params.genome} {params.threads}"


rule filter:
    input:
        directory(config["project_name"])
    output:
        directory(config["project_name"] + config["filter_suffix"])
    params:
        project=config["project_name"],
        genome=config["genome"],
        threads=config["threads"],
        suffix=config["filter_suffix"],
        min_tss_enrichment =config['min_tss_enrichment'], 
        min_nfrags = config['min_nfrags'], 
        min_gex_ngenes = config['min_gex_ngenes'], 
        max_gex_ngenes =config['max_gex_ngenes'],
        min_gex_numi =config['min_gex_numi'],
        max_gex_numi =config['max_gex_numi']
    shell:
        "Rscript src/filter.R --project_name {params.project} --genome {params.genome} --threads {params.threads} --suffix {params.suffix} --min_tss_enrichment {params.min_tss_enrichment} --min_nfrags {params.min_nfrags} --min_gex_ngenes {params.min_gex_ngenes} --max_gex_ngenes {params.max_gex_ngenes} --min_gex_numi {params.min_gex_numi} --max_gex_numi {params.max_gex_numi}"



rule addUMAP:
    input:
        directory(config["project_name"] + config["filter_suffix"])
    output:
        directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"])
    params:
        project=config["project_name"] + config["filter_suffix"],
        genome=config["genome"],
        threads=config["threads"],
        suffix=config["umap_suffix"]
    shell:
        "Rscript src/addUMAP.R --project_name {params.project} --genome {params.genome} --threads {params.threads} --suffix {params.suffix}"


rule plot:
    input:
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"]
    params:
        project=config["project_name"] + config["filter_suffix"] + config["umap_suffix"]
    output:
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_Combined_UMAP.pdf",
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_ATAC_UMAP.pdf",
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_RNA_UMAP.pdf"
    shell:
        "Rscript src/plot.R --project_name {params.project}"


rule plotMarkers:
    input:
        dir=directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"]),
    output:
        out= "markers_done.txt" 
    params:
        project=config["project_name"] + config["filter_suffix"] + config["umap_suffix"],
        markers=config['MarkerGenes']
    shell:
        """
        Rscript src/plotMarkers.R --project_name {params.project} --markers {params.markers} 
        touch {output}
        """

rule annotate:
    input:
       dir=directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"])
    output:
       dir=directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + config["annotation_suffix"]),
       annotated_png=config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + "_annotated.png"
    params:
          annotations=config['ANNOTATIONS'],
          suffix=config['annotation_suffix']
    shell:
      "Rscript src/annotate.R --project_name {input.dir} --annotation_file {params.annotations} --suffix {params.suffix} --remove_clusters 1 2 3 4 14"


rule markerPeaks: 
    input: 
     dir = config["project_name"] + config["filter_suffix"] + \
              config["umap_suffix"] + config["annotation_suffix"] 
    params: 
       config['peaks_suffix'] 
    output: 
      dir=directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + config["annotation_suffix"]+config['peaks_suffix']) 
    shell: 
        "Rscript src/markerPeaks.R  --project_name {input} --suffix {params}" 



rule motifEnrich: 
    input: 
        dir = config["project_name"] + config["filter_suffix"] + \
              config["umap_suffix"] + config["annotation_suffix"] + \
              config["peaks_suffix"]
    params: 
      suffix= config['motif_suffix']
    output: 
      dir = directory(config["project_name"] + config["filter_suffix"] + config["umap_suffix"] + config["annotation_suffix"]+config['peaks_suffix']+config['motif_suffix']) 
    shell: 
       "Rscript src/motifEnrich.R --project_name {input.dir} --motif_set cisbp --suffix {params.suffix}"        

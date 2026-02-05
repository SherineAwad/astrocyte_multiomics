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
        config["project_name"] + config["filter_suffix"] + config["umap_suffix"]+ "_annotated.png"

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
        suffix=config["filter_suffix"]
    shell:
        "Rscript src/filter.R --project_name {params.project} --genome {params.genome} --threads {params.threads} --suffix {params.suffix}"


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
      "Rscript src/annotate.R --project_name {input.dir} --annotation_file {params.annotations} --suffix {params.suffix}"





 

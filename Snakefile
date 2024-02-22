configfile: "scripts/setup/config.json"
configfile: "scripts/setup/fract_config.json"
import os
import re

READS = ["_1", "_2"]
F_SAMPLES = config["samples_fract"]

rule all:
    input:
        "data/raw/fastq/multiqc/multiqc_report.html",
        "data/processed/bulk/rau_fractions_gse.RData",
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png",
        "results/7_plot_comps/sample_comps_relative.png",
        "results/10_plot_de/volcano_adjusted.png",
        "results/5_findMarkers/cell_clusters.png",
        "results/5_findMarkers/marker_specificity.png",
        "data/raw/anno/gencode.vM34.annotation.gtf.gz"

rule load_index:
    output: 
        "data/raw/anno/gencode.vM34.transcripts.fa.gz",
    shell: 
        "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.transcripts.fa.gz -P data/raw/anno"
rule load_gtf:
    output: 
        "data/raw/anno/gencode.vM34.annotation.gtf.gz"
    shell: 
        "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz -P data/raw/anno"
rule salmon_index:
    input:
        "data/raw/anno/gencode.vM34.transcripts.fa.gz"
    output:
        "data/raw/anno/gencode.vM34.salmon"
    shell:
        "salmon index --gencode -p 12 -t {input} -i {output}"
rule salmon_quant:
    input:
        r1 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_1.fastq.gz",
        r2 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_2.fastq.gz",
        index = "data/raw/anno/gencode.vM34.salmon"
    output:
        "data/raw/fastq/{F_SAMPLES}/quant.sf"
    params:
        dir = "data/raw/fastq/{F_SAMPLES}"
    shell:
        "salmon quant -i {input.index} -l A -p 12 --gcBias "
        "--numGibbsSamples 20 --thinningFactor 100 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"
rule fastqc:
    input:
        "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}.fastq.gz"
    output:
        html = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}_fastqc.html",
        zip = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}_fastqc.zip"
    shell:
        "fastqc {input} --outdir=$(dirname {input})"
rule multiqc:
    input:
        expand(["data/raw/fastq/{sample}/quant.sf",
                "data/raw/fastq/{sample}/{sample}{read}_fastqc.html"],
                sample=F_SAMPLES, read=READS)
    output:
        "data/raw/fastq/multiqc/multiqc_report.html"
    shell:
        "multiqc . -o data/raw/fastq/multiqc"
rule tximport:
    input:
        "data/raw/anno/gencode.vM34.annotation.gtf.gz",
        expand(["data/raw/fastq/{sample}/quant.sf",
                "data/raw/fastq/{sample}/{sample}{read}_fastqc.html"],
                sample=F_SAMPLES, read=READS)
    output:
        "data/processed/bulk/rau_fractions_ensembl.csv"
    shell:
        "Rscript scripts/0_transcripts_to_genes.R"
rule ensb2gene:
    input:
        "data/processed/bulk/rau_fractions_ensembl.csv"
    output:
        "data/processed/bulk/rau_fractions_gene.csv"
    shell:
        "Rscript scripts/4_1_ens_to_gene.R"
rule load_sn:
    output: 
        "data/processed/single_cell/unprocessed/{samples}.h5seurat"
    shell:
        "Rscript scripts/1_load_sn.R {wildcards.samples}"
rule ambient_doublets:
    input:
        "data/processed/single_cell/unprocessed/{samples}.h5seurat"
    output: 
        "data/processed/single_cell/no_doublets/{samples}_no_doublets.h5seurat"
    shell:
        "Rscript scripts/2_ambient_doublets.R {wildcards.samples}"
rule merge_sn:
    input:
        expand("data/processed/single_cell/no_doublets/{samples}_no_doublets.h5seurat", 
               samples = config["samples"])
    output: 
        "results/3_merge_sn/cluster_features_3.png",
        "data/processed/single_cell/merged_no_doublets.h5seurat"
    shell:
        "Rscript scripts/3_merge_sn.R"
rule clean_bulk:
    input:
        "data/processed/bulk/rau_fractions_gene.csv",
        "data/raw/jensen/jensen_counts_correct.xlsx"
    output: 
        "data/processed/bulk/all_counts.csv",
        "data/processed/bulk/jensen_pheno.csv",
        "data/processed/bulk/jensen_bulk_clean.csv"
    shell:
        "Rscript scripts/4_clean_bulk.R"
rule findMarkers:
    input:
        "data/processed/single_cell/merged_no_doublets.h5seurat",
        "data/processed/bulk/all_counts.csv"
    output: 
        "results/5_findMarkers/cell_clusters.png"
    shell:
        "Rscript scripts/5_findMarkers.R"
rule plotMarkers:
    input:
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "data/processed/single_cell/cluster_markers.csv",
    output: 
        "results/5_findMarkers/marker_specificity.png"
    shell:
        "Rscript scripts/5_1_plot_markers.R"
rule deconvolute:
    input:
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "data/raw/rau_fractions/celltype_pheno.csv",
        "data/processed/bulk/jensen_pheno.csv",
        "data/processed/single_cell/cluster_markers.csv",
        "data/processed/bulk/all_counts.csv"
    output: 
        "data/processed/compositions/whole_samples.csv",
        "data/processed/compositions/fraction_samples.csv"
    shell:
        "Rscript scripts/6_deconvolute.R"
rule plot_comps:
    input:
        "data/processed/compositions/whole_samples.csv",
        "data/processed/compositions/fraction_samples.csv"
    output: 
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png",
        "results/7_plot_comps/sample_comps_relative.png"
    shell:
        "Rscript scripts/7_plot_comps.R"
rule dirichlet:
    input:
        "data/processed/compositions/whole_samples.csv"
    output: 
        "data/processed/models/dirichelet_coefficients.csv",
        "results/8_dirichlet/dirichlet_coeff.png"
    shell:
        "Rscript scripts/8_dirichlet.R"
rule diffential_expression:
    input:
        "data/processed/compositions/whole_samples.csv",
        "data/processed/bulk/jensen_pheno.csv",
        "data/processed/bulk/jensen_bulk_clean.csv"
    output: 
        "data/processed/models/adjusted_de.csv"
    shell:
        "Rscript scripts/9_differential_expression.R"
rule plot_de:
    input:
        "data/processed/models/adjusted_de.csv"
    output: 
        "results/10_plot_de/volcano_adjusted.png"
    shell:
        "Rscript scripts/10_plot_de.R"

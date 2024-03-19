configfile: "scripts/setup/config.json"
configfile: "scripts/setup/fastq_config.json"
import os
import re

READS = ["_1", "_2"]
F_SAMPLES = config["samples_fastq"]
SN_SAMPLES = ["b6_1", "b6_2"]

rule all:
    input:
        "data/raw/multiqc/multiqc_report.html",
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png",
        "results/8_dirichlet/dirichlet_coeff.png",
        "results/10_plot_de/volcano_adjusted.png",
        "results/5_findMarkers/cell_clusters.png",
        "results/5_findMarkers/marker_specificity.png",
        "data/raw/anno/gencode.vM34.annotation.gtf.gz",
        "data/processed/bulk/all_counts.csv",
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "results/11_clusterProfiler/adj_interaction_clusters.png",
        expand("data/processed/single_cell/unprocessed/{sn_sample}.h5seurat", sn_sample=SN_SAMPLES)

#rule load_index:
#    output: 
#        "data/raw/anno/gencode.vM34.transcripts.fa.gz"
#    shell: 
#        "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.transcripts.fa.gz -P data/raw/anno"
rule load_decoy:
    output: 
        "data/raw/anno/salmon_sa_index"
    shell: 
        "wget http://refgenomes.databio.org/v3/assets/archive/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/salmon_sa_index?tag=default -P data/raw/anno -O salmon_sa_index"
rule load_gtf:
    output: 
        "data/raw/anno/gencode.vM34.annotation.gtf.gz"
    shell: 
        "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz -P data/raw/anno"
rule salmon_index:
#    input:
#        trans = "data/raw/anno/decoy/gentrome.fa",
#        decoy = "data/raw/anno/decoy/decoys.txt"
    output:
        directory("data/raw/anno/decoyaware.vM34.salmon")
    shell:
        "salmon index --gencode -t data/raw/anno/decoy/gentrome.fa -d data/raw/anno/decoy/decoys.txt -i {output} -k 25"
rule salmon_quant:
    input:
        r1 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_1.fastq.gz",
        r2 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_2.fastq.gz",
        index = "data/raw/anno/decoyaware.vM34.salmon"
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
        "data/raw/multiqc/multiqc_report.html"
    shell:
        "multiqc . -o data/raw/multiqc -f"
rule tximport:
    input:
        "data/raw/anno/gencode.vM34.annotation.gtf.gz",
        expand(["data/raw/fastq/{sample}/quant.sf",
                "data/raw/fastq/{sample}/{sample}{read}_fastqc.html"],
                sample=F_SAMPLES, read=READS)
    output:
        "data/processed/bulk/all_bulk_ensembl.csv"
    shell:
        "Rscript scripts/0_transcripts_to_genes.R"
rule ensb2gene:
    input:
        "data/processed/bulk/all_bulk_ensembl.csv"
    output:
        "data/processed/bulk/all_bulk_gene.csv"
    shell:
        "Rscript scripts/4_1_ens_to_gene.R"
rule load_sn:
    output: 
        "data/processed/single_cell/unprocessed/{sn_sample}.h5seurat"
    resources:
        mem_mb=6000
    shell:
        "Rscript scripts/1_load_sn.R {wildcards.sn_sample}"
rule ambient_doublets:
    input:
        "data/processed/single_cell/unprocessed/{sn_sample}.h5seurat"
    output: 
        "data/processed/single_cell/no_doublets/{sn_sample}_no_doublets.h5seurat"
    resources:
        mem_mb=8000
    shell:
        "Rscript scripts/2_ambient_doublets.R {wildcards.sn_sample}"
rule merge_sn:
    input:
        expand(["data/processed/single_cell/no_doublets/{sn_sample}_no_doublets.h5seurat"],
                sn_sample=SN_SAMPLES)
    output:
        "data/processed/single_cell/merged_no_doublets.h5seurat"
    resources:
        mem_mb=16000
    shell:
        "Rscript scripts/3_merge_sn.R"
rule clean_bulk:
    input:
        "data/processed/bulk/all_bulk_gene.csv"
    output: 
        "data/processed/bulk/all_counts.csv",
        "data/processed/bulk/pheno_table.csv"
    shell:
        "Rscript scripts/4_clean_bulk.R"
rule findMarkers:
    input:
        "data/processed/single_cell/merged_no_doublets.h5seurat",
        "data/processed/bulk/all_counts.csv"
    output: 
        "results/5_findMarkers/cell_clusters.png",
        "data/processed/single_cell/cluster_markers.csv",
        "data/processed/single_cell/celltype_labeled.h5seurat"
    shell:
        "Rscript scripts/5_findMarkers.R"
rule plotMarkers:
    input:
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "data/processed/single_cell/cluster_markers.csv"
    output: 
        "results/5_findMarkers/marker_specificity.png"
    shell:
        "Rscript scripts/5_1_plot_markers.R"
rule deconvolute:
    input:
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "data/processed/bulk/pheno_table.csv",
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
        "results/7_plot_comps/sample_comps.png"
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
        "data/processed/bulk/pheno_table.csv",
        "data/processed/bulk/all_counts.csv"
    output: 
        "data/processed/models/adjusted_de_interaction.RDS",
        "data/processed/models/unadjusted_de_interaction.RDS"
    shell:
        "Rscript scripts/9_differential_expression.R"
rule plot_de:
    input:
        "data/processed/models/adjusted_de_interaction.RDS"
    output: 
        "results/10_plot_de/volcano_adjusted.png"
    shell:
        "Rscript scripts/10_plot_de.R"
rule plot_upset:
    input:
        "data/processed/models/unadjusted_de_interaction.RDS"
    output:
        "results/10_plot_de/upset_unadj.png"
    shell:
        "Rscript scripts/10_1_upset_plot.R"
rule gene_ont:
    input:
        "data/processed/models/adjusted_de_interaction.RDS",
        "data/processed/models/adjusted_de_interaction.RDS"
    output: 
        "results/11_clusterProfiler/unadj_interaction_clusters.png",
        "results/11_clusterProfiler/adj_interaction_clusters.png"
    shell:
        "Rscript scripts/11_clusterProfiler.R"
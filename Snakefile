configfile: "scripts/setup/config.json"
configfile: "scripts/setup/fastq_config.json"
import os
import re

READS = ["_1", "_2"]
F_SAMPLES = config["samples_fastq"]
SN_SAMPLES = ["b6_1", "b6_2"]
MODEL_TYPE = ["adjusted", "unadjusted"]

rule all:
    input:
        "data/raw/multiqc/multiqc_report.html",
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png",
        "results/8_dirichlet/dirichlet_coeff.png",
        "results/10_plot_de/volcano_adjusted.png",
        #"results/10_plot_de/upset_unadj.png",
        "results/5_findMarkers/cell_clusters.png",
        "results/5_findMarkers/marker_specificity.png",
        "data/processed/single_cell/celltype_labeled.h5seurat",
        expand("data/processed/pathway_genesets/go_{model_type}_any_p.RDS", model_type=MODEL_TYPE),
        expand("data/processed/single_cell/unprocessed/{sn_sample}.h5seurat", sn_sample=SN_SAMPLES)

rule load_transcript:
    output:
        "data/raw/anno/Mus_musculus.GRCm39.cdna.all.fa.gz"
    shell:
        "wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz -P data/raw/anno"
rule load_toplevel:
    output:
        "data/raw/anno/Mus_musculus.GRCm39.dna.toplevel.fa.gz"
    shell:
        "wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz -P data/raw/anno/"
rule make_decoy_1:
    input:
        fa="data/raw/anno/Mus_musculus.GRCm39.dna.toplevel.fa.gz"
    output:
        decoy="data/raw/anno/decoy.txt"
    resources:
        mem_mb=10000
    shell:
        """
        grep '^>' <(gunzip -c {input.fa}) | cut -d ' ' -f 1 > {output.decoy}
        """
rule make_decoy_2:
    input:
        "data/raw/anno/decoy.txt"
    output:
        "data/raw/anno/decoy.cleaned.txt"
    shell:
        """
        sed 's/>//g' {input} > {output}
        """
rule make_gentrome:
    input:
        top = "data/raw/anno/Mus_musculus.GRCm39.dna.toplevel.fa.gz",
        transc ="data/raw/anno/Mus_musculus.GRCm39.cdna.all.fa.gz"
    output:
        "data/raw/anno/gentrome.fa.gz"
    shell:
        "cat {input.transc} {input.top} > data/raw/anno/gentrome.fa.gz"
rule salmon_index:
    input:
        fa = "data/raw/anno/gentrome.fa.gz",
        decoy = "data/raw/anno/decoy.cleaned.txt"
    resources:
        mem_mb=32000
    output:
        directory("data/raw/anno/decoy_Mus_musculus.GRCm39.salmon")
    shell:
        "salmon index -t {input.fa} -k 31 --keepFixedFasta -p 16 -i {output} -d {input.decoy}"
rule salmon_quant:
    input:
        r1 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_1.fastq.gz",
        r2 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_2.fastq.gz",
        index = "data/raw/anno/decoy_Mus_musculus.GRCm39.salmon"
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
        "Rscript scripts/04_1_ens_to_gene.R"
rule load_sn:
    output: 
        "data/processed/single_cell/unprocessed/{sn_sample}.h5seurat"
    resources:
        mem_mb=6000
    shell:
        "Rscript scripts/01_load_sn.R {wildcards.sn_sample}"
rule ambient_doublets:
    input:
        "data/processed/single_cell/unprocessed/{sn_sample}.h5seurat"
    output: 
        "data/processed/single_cell/no_doublets/{sn_sample}_no_doublets.h5seurat"
    resources:
        mem_mb=8000
    shell:
        "Rscript scripts/02_ambient_doublets.R {wildcards.sn_sample}"
rule merge_sn:
    input:
        expand(["data/processed/single_cell/no_doublets/{sn_sample}_no_doublets.h5seurat"],
                sn_sample=SN_SAMPLES)
    output:
        "data/processed/single_cell/merged_no_doublets.h5seurat"
    resources:
        mem_mb=16000
    shell:
        "Rscript scripts/03_merge_sn.R"
rule clean_bulk:
    input:
        "data/processed/bulk/all_bulk_gene.csv"
    output: 
        "data/processed/bulk/all_counts.csv",
        "data/processed/bulk/pheno_table.csv"
    shell:
        "Rscript scripts/04_clean_bulk.R"
rule findMarkers:
    input:
        "data/processed/single_cell/merged_no_doublets.h5seurat",
        "data/processed/bulk/all_counts.csv"
    output: 
        "results/5_findMarkers/cell_clusters.png",
        "data/processed/single_cell/cluster_markers.csv",
        "data/processed/single_cell/celltype_labeled.h5seurat"
    shell:
        "Rscript scripts/05_findMarkers.R"
rule plotMarkers:
    input:
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "data/processed/single_cell/cluster_markers.csv"
    output: 
        "results/5_findMarkers/marker_specificity.png"
    shell:
        "Rscript scripts/05_1_plot_markers.R"
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
        "Rscript scripts/06_deconvolute.R"
rule plot_comps:
    input:
        "data/processed/compositions/whole_samples.csv",
        "data/processed/compositions/fraction_samples.csv"
    output: 
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png"
    shell:
        "Rscript scripts/07_plot_comps.R"
rule dirichlet:
    input:
        "data/processed/compositions/whole_samples.csv"
    output: 
        "data/processed/models/dirichelet_coefficients.csv",
        "results/8_dirichlet/dirichlet_coeff.png"
    shell:
        "Rscript scripts/08_dirichlet.R"
rule diffential_expression:
    input:
        "data/processed/compositions/whole_samples.csv",
        "data/processed/bulk/pheno_table.csv",
        "data/processed/bulk/all_counts.csv"
    output: 
        "data/processed/models/adjusted_de_interaction.RDS",
        "data/processed/models/unadjusted_de_interaction.RDS"
    shell:
        "Rscript scripts/09_differential_expression.R"
rule plot_de:
    input:
        "data/processed/models/adjusted_de_interaction.RDS"
    output: 
        "results/10_plot_de/volcano_adjusted.png"
    shell:
        "Rscript scripts/10_plot_de.R"
#rule plot_upset:
#    input:
#        "data/processed/models/unadjusted_de_interaction.RDS"
#    output:
#        "results/10_plot_de/upset_unadj.png"
#    shell:
#        "Rscript scripts/10_1_upset_plot.R"
rule gene_ont:
    input:
        "data/processed/models/adjusted_de_interaction.RDS",
        "data/processed/models/unadjusted_de_interaction.RDS"
    output: 
        "data/processed/pathway_genesets/go_{model_type}_any_p.RDS"
    shell:
        "Rscript scripts/11_clusterProfiler.R {wildcards.model_type}"
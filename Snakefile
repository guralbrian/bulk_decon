configfile: "config.json"
import os
import re

READS = ["_R1", "_R2"]

# This is a function to flexibly list paths of samples following Christoph Rau's fractions nomenclature
def get_sample_paths():
    sample_paths = []
    for root, dirs, files in os.walk("data/raw/fastq/"):
        for file in files:
            if re.match(r".*_S\d+_L002_R[12]_001\.fastq\.gz$", file):
                # Get the full path without the .fastq.gz extension
                full_path = os.path.join(root, file)
                sample_path = re.sub(r"\.fastq\.gz$", "", full_path)
                sample_paths.append(sample_path)
    return sample_paths

def get_sample_paths():
    sample_set = set()
    pattern = re.compile(r"(.*/.*_S\d+_L002_R*")
    for root, dirs, files in os.walk("data/raw/fastq/"):
        for file in files:
            match = pattern.match(file)
            if match:
                # Extracts the directory and the part of the filename up to _R1 or _R2
                sample_set.add(match.group(1))
    return list(sample_set)

# Generate the expected HTML file names
sample_paths = get_sample_paths()
expected_html_files = [sample_path + "_fastqc.html" for sample_path in sample_paths]

rule all:
    input:
        "data/raw/fastq/multiqc/multiqc_report.html",
        expand("data/processed/single_cell/no_doublets/{samples}_no_doublets.h5seurat", 
               samples = config["samples"]),
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png",
        "results/7_plot_comps/sample_comps_relative.png",
        "data/processed/models/dirichelet_coefficients.csv",
        "results/10_plot_de/volcano_adjusted.png",
        "data/raw/anno/gencode.vM34.salmon/"

rule fastqc:
    input:
        "data/raw/fastq/{sample_path}.fastq.gz"
    output:
        html = "{sample_path}_fastqc.html",
        zip = "{sample_path}_fastqc.zip"
    shell:
        "fastqc {input} --outdir=$(dirname {input})"
rule load_index:
    output: 
        "data/raw/anno/gencode.vM34.transcripts.fa.gz"
    shell: 
        "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.transcripts.fa.gz -P data/raw/anno"
rule salmon_index:
    input:
        "data/raw/anno/gencode.vM34.transcripts.fa.gz"
    output:
        "data/raw/anno/gencode.vM34.salmon"
    shell:
        "salmon index --gencode -p 12 -t {input} -i {output}"
rule salmon_quant:
    input:
        r1 = "fastq/{sample}_1.fastq.gz",
        r2 = "fastq/{sample}_2.fastq.gz",
        index = "/proj/milovelab/anno/gencode.v38-salmon_1.4.0"
    output:
        "quants/{sample}/quant.sf"
    params:
        dir = "quants/{sample}"
    shell:
        "{SALMON} quant -i {input.index} -l A -p 12 --gcBias "
        "--numGibbsSamples 20 --thinningFactor 100 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"
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
        "data/raw/rau_fractions/celltype_counts.csv",
        "data/raw/rau_fractions/celltype_pheno.csv",
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
        "data/processed/single_cell/celltype_labeled.h5seurat",
        "data/processed/single_cell/cluster_markers.csv"
    shell:
        "Rscript scripts/5_findMarkers.R"
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

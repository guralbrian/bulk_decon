configfile: "config.json"
rule all:
    input:
        expand("data/processed/single_cell/no_doublets/{samples}_no_doublets.h5seurat", 
               samples = config["samples"]),
        "results/7_plot_comps/pure_cell_types.png",
        "results/7_plot_comps/sample_comps.png",
        "results/7_plot_comps/sample_comps_relative.png",
        "data/processed/models/dirichelet_coefficients.csv",
        "results/10_plot_de/volcano_adjusted.png"
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

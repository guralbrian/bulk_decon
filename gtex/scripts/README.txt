1_gtex_load.R

Loads, cleans, and saves the snRNAseq and bulk RNAseq datasets
Accesses data by URL
Converts AnnData file to .H5Seurat file


2_gtex_de_genes.R

Aggregates subject-matched single nucleus data into pseudo-bulk counts matrices, then compares them to actual bulk.
Identifies differentially expressed genes between bulk and single nucleus RNAseq
Can adjust several parameters including prefiltering of nuclei by RNA count/features, mitocondrial/ribosomal gene content, and doublet tags. Also can change log fold change threshold for DE and directionality (i.e. greaterAbs, lessAbs, etc.)
Saves a seurat obj w/o DE genes


3_gtex_decon.R

Preforms bulk deconvolution on several hundred GTEx samples after processing snRNAseq
Clusters snRNAseq, excludes clusters on adjustable criteria (i.e. mitocondrial content, doublet scores, etc.), and repeats clustering after exclusion of nuclei
Outputs stacked barcharts of composition estimates

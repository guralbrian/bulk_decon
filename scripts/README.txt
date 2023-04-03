OVERVIEW OF /SCRIPTS
/general contains scripts which were not intended to be used as functions. 
  They are used to do the actual heavy lifting to achieve the goals of this project: to estimate cellular composition from bulk RNAseq
  
/functions contains scripts which define functions. The scripts in /general depend upon these functions. 


WORKFLOW

gtex_load.R
 -load bulk and sn datasets
 -converts H5AD to Seurat
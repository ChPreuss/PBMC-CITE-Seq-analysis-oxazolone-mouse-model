# PBMC-CITE-SEQ-data-oxazolone-mouse-model

## Single cell CITE-Seq data from human PBMC cells in oxazolone treated mouse ears for Ucenprubart
The provided code was used to perform CITE-Seq analysis of FACS-sorted live human CD45+ cells isolated from the ears of huNOG-EXL mice one day after the third oxazolone challenge in a model of atopic dermatitis. Mice were previously treated with CD200R agonist mAb (Ucenprubart), dupilumab, or an isotype control mAb (two biological replicates per treatment group). The data used for the analysis has been deposited under the NCBI GEO submission number GSE220685.

## Data analysis 
This repository contains two scripts and a folder with cell classification files for human and mouse cells: <br>
1.) The first script was used to process the raw CITE-Seq data and generate an annotated Seurat object for downstream analysis. <br>
2.) The second script was used to perform differentially gene expression and pathway analysis of the dataset. It contains several functions to plot the results and provides summary statistics.<br>
3.) A folder with cell classifications for mouse and human cells for each sample from the ears of the huNOG-EXL mice is provided. The cell ranger pipeline was used to perform a multi-genome analysis in which the algorithm classified human cells as the 10th percentile of all barcodes where human UMI counts are greater than mouse UMI counts.  The generated cell ranger secondary analysis output ‘gem_classification.csv’ files were utilized to remove mouse cells from each sample. 

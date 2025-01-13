# Define R libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library("org.Hs.eg.db")
library(psych)
library(corrplot)
library(eulerr)
theme_set(theme_cowplot())


############################################################################
# Importing the data from cellranger into R
############################################################################

# Define sample names to be imported 
sample_names  <- c("CD200REar-1",
                   "CD200REar-2",
                   "DupilumabEar-1",
                   "DupilumabEar-2",
                   "IsotypeEar-1",
                   "IsotypeEar-2")


# Import samples into R by providing the directory containing the count matrices
# The directory should contain count matrix files from cell ranger and the cell classifications
# The raw CITE-Seq expression data can be found in the NBI-GEO dataset accompanying the manuscript

cite_seq_data_dirs <-
  Map(function(x) glue::glue("Read in count matrices by samples", sample_name=x), sample_names)

# Import cell classfications from cell ranger for the "barnyard" experiment where cells from mouse and human can be mixed and analyzed together.
# This allows a subset of multiplets to be detected on the basis that some reads with a given cell barcode align to the human reference genome, 
# and a different set of reads with the same cell barcode will align to the mouse reference genome.
# The actual classifications from the experiment can be found in the gem_classification.csv files provided in the github repo folder.
# The cell range algorithm classifies human cells as the 10th percentile of all barcodes where human UMI counts > mouse UMI counts.

cell_classifications_by_mouse_or_human <-
  Map(function(x) glue::glue("Read in count matrices cell classsifications by sample_gem_classification.csv", sample_name=x), sample_names)


## Map samples to human and mouse genomes
input_data <- Map(Seurat::Read10X, cite_seq_data_dirs)

cell_classifications_by_mouse_or_human_list <- Map(readr::read_csv, cell_classifications_by_mouse_or_human)

make_seurat_objects <- function(input, classification, Cell_annot){
  so <- Seurat::CreateSeuratObject(counts=input$`Gene Expression`,min.cells = 3)
  adt_assay<-CreateAssay5Object(counts = input$`Antibody Capture`)
  so[["ADT"]] <- adt_assay
  classification <- as.data.frame(classification)
  rownames(classification) <- classification$barcode
  classification$barcode  <- NULL
  so <- Seurat::AddMetaData(so, metadata=classification)
  return(so)
}

seurat_objects <- Map(make_seurat_objects, input_data, cell_classifications_by_mouse_or_human_list)

filter_so <- function(so, type){
  if(type=="hg19"){
    human_rna_features <- rownames(so@assays$RNA@features)[rownames(so@assays$RNA@features) %>% stringr::str_detect("^hg19-")]
    human_protein_features <-rownames(so@assays$ADT@features)[rownames(so@assays$ADT@features) %>%  stringr::str_detect("^Human-")]
    features_to_filter <-  append(human_rna_features, human_protein_features)
    seurat_object_filtered <- subset(so, subset=call=="hg19", features=features_to_filter)
  }
  if(type=="mm10"){
    mouse_rna_features <- rownames(so@assays$RNA@features)[rownames(so@assays$RNA@features) %>% stringr::str_detect("^mm10-")]
    mouse_protein_features <- rownames(so@assays$ADT@features)
    features_to_filter <- append(mouse_rna_features, mouse_protein_features)
    seurat_object_filtered <- subset(so, subset=call=="mm10", features=features_to_filter)
  }
  return(seurat_object_filtered)
}

# Filter human and mouse samples

seurat_objects_human <- Map(function(x) filter_so(x, "hg19"), seurat_objects)
seurat_objects_mouse <- Map(function(x) filter_so(x, "mm10"), seurat_objects)

# Plot Human and Mouse Cell Counts
human.mouse.CD200R<-bind_rows(cell_classifications_by_mouse_or_human_list, .id = "column_label")

human.mouse.plot<-human.mouse.CD200R %>% 
  group_by(column_label) %>%
  summarize(Human_Mouse = table(call),Anno=names(table((call))))

ggplot(human.mouse.plot, aes(fill=Anno, y=Human_Mouse, x=column_label)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette="Dark2") +
  theme_minimal()


## Recreate gene names without prepended hg19 and retain only human cells
rename_so <- function(so){
  count.data <- GetAssayData(object = so[["RNA"]], slot = "counts")
  count.data <- as.matrix(count.data)
  adt.data <- GetAssayData(object = so[["ADT"]], slot = "counts") 
  #  adt.data <- as.matrix(adt.data)
  ## Rename rows
  rownames(count.data) <- rownames(count.data) %>% stringr::str_remove("hg19-")
  ## Generate new Seurat object.
  meta_data  <- so@meta.data
  new_obj <- CreateSeuratObject(
    count.data,
    project = "SeuratProject",
    assay = "RNA",
    min.cells = 3,
    min.features = 0,
    names.field = 1,
    names.delim = "_",
    meta.data = meta_data)
  new_obj[["ADT"]] <- CreateAssay5Object(counts=adt.data)
  return(new_obj)
}

seurat_objects_human_renamed <- Map(rename_so, seurat_objects_human)
seurat_objects_human <- seurat_objects_human_renamed


## Filter cells with MT DNA content > 75%
filter_cells_with_mt_content <- function(so, max_mt_percent=25){
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(
    so,
    subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 25 )
  return(so)
}


seurat_objects_human <- Map(filter_cells_with_mt_content, seurat_objects_human)


##  Merging Seurat Objects and doing simple QC steps
# Merge datasets into one single seurat object and assign labels for treatment
seurat_objects_human[[1]]$stim<-"Ucenprubart"
seurat_objects_human[[2]]$stim<-"Ucenprubart"
seurat_objects_human[[3]]$stim<-"Dupilumab"
seurat_objects_human[[4]]$stim<-"Dupilumab"
seurat_objects_human[[5]]$stim<-"Isotype"
seurat_objects_human[[6]]$stim<-"Isotype"


alldata <- merge(seurat_objects_human[[1]], c(seurat_objects_human[[2]], seurat_objects_human[[3]],seurat_objects_human[[4]],seurat_objects_human[[5]],seurat_objects_human[[6]]), add.cell.ids = c("CD200REar-1", "CD200REar-2", "DupilumabEar-1","DupilumabEar-2","IsotypeEar-1","IsotypeEar-2"))

alldata$orig.ident<-gsub("\\_.*","",names(alldata$orig.ident))
# Read from ribosomal proteins and mitochondrial proteins
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
feats <- c("nFeature_RNA", "nCount_RNA","nFeature_Protein","nCount_Protein", "percent_mito", "percent_ribo", "percent_hb")

# Plot features
#VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + NoLegend()

# Plot number of genes vs. number of molecules detected in a cell
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)


## Integrate using highest dimensional data. In this case: RNA
## Taken from:https://satijalab.org/seurat/articles/integration_introduction.html

# split the RNA measurements into three layers one for each stimulation
alldata<-JoinLayers(alldata)
alldata[["RNA"]] <- split(alldata[["RNA"]], f = alldata$stim)

# run standard Seurat analysis workflow
seurat_objects_human <- NormalizeData(alldata)
seurat_objects_human <- FindVariableFeatures(seurat_objects_human)
seurat_objects_human <- ScaleData(seurat_objects_human)

# Run PCA for the first 30 principal components and identify cell clusters
seurat_objects_human <- RunPCA(seurat_objects_human) 
seurat_objects_human <- FindNeighbors(seurat_objects_human, dims = 1:30, reduction = "pca")
seurat_objects_human <- FindClusters(seurat_objects_human, resolution = 0.5, cluster.name = "unintegrated_clusters")

# Run UMAP
seurat_objects_human <- RunUMAP(seurat_objects_human, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Plot unintegrated clusters
DimPlot(seurat_objects_human, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

# Run the Seurat v5 integration procedure to return a single dimensional reduction that captures the shared sources of variance across the three stimulation layers, 
# This returns a dimensional reduction object (i.e. integrated.cca) which can be used for visualization and unsupervised clustering analysis. 
seurat_objects_human <- IntegrateLayers(object = seurat_objects_human, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                        verbose = FALSE)

# Re-join RNA layers after integration
seurat_objects_human[["RNA"]] <- JoinLayers(seurat_objects_human[["RNA"]])
seurat_objects_human[["ADT"]] <- JoinLayers(seurat_objects_human[["ADT"]])


# Re-run clustering and UMAP for integrated dataset
seurat_objects_human <- FindNeighbors(seurat_objects_human, reduction = "integrated.cca", dims = 1:30)
seurat_objects_human <- FindClusters(seurat_objects_human, resolution = 0.5)
seurat_objects_human <- RunUMAP(seurat_objects_human, dims = 1:30, reduction = "integrated.cca")

#Clean up objects to save space for next steps
rm(alldata)
rm(seurat_objects)
gc()


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. 
# Based on this markers, cell cycling genes are regressed out from analysis
# Signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase among proliferating cells 
#(which are often uninteresting), will be regressed out of the data
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_objects_human <- CellCycleScoring(seurat_objects_human, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_objects_human$CC.Difference <- seurat_objects_human$S.Score - seurat_objects_human$G2M.Score
#seurat_objects_human <- ScaleData(seurat_objects_human, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_objects_human))
seurat_objects_human <- ScaleData(seurat_objects_human, vars.to.regress = "CC.Difference", features = rownames(seurat_objects_human))

# Visualization of Seurat object highlighting cell clusters grouped by stimulation
DimPlot(seurat_objects_human, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# Normalize and scale protein data from ADT layer and rename Seurat Object after integration
human_data_integrated<-seurat_objects_human
# Assign data slot based on ADT counts
human_data_integrated@assays$ADT$data<-human_data_integrated@assays$ADT$counts
human_data_integrated <- NormalizeData(human_data_integrated, assay = "ADT", normalization.method = "CLR")
human_data_integrated <- ScaleData(human_data_integrated, assay = "ADT")

############################################################################
# Annotating Cells using SingleR Human Primary Cell Atlas Data as Reference
############################################################################
# See https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

library(SingleR) # Required for cell type annotations from Human Primary Cell Atls
hpca.se <- HumanPrimaryCellAtlasData()

DefaultAssay(human_data_integrated) <- "RNA"

pred.hesc_label_fine <-
  SingleR(test =
            as.SingleCellExperiment(human_data_integrated),
          ref = hpca.se,
          assay.type.test=1,
          labels = hpca.se$label.fine)

pred.hesc_label <-
  SingleR(test =
            as.SingleCellExperiment(human_data_integrated),
          ref = hpca.se,
          assay.type.test=1,
          labels = hpca.se$label.main)

pred_cells <- pred.hesc_label
pred_cells_fine <- pred.hesc_label_fine

############################################################################
# Cell type annotation
############################################################################
# Curated annotation of seruat clusters based on SingleR.labels / CITE-Seq / Publicly Available data
Idents(human_data_integrated)<-"seurat_clusters"
human_data_integrated <- RenameIdents(human_data_integrated, 
                                      `0` = "CD4+_T-cells", 
                                      `1` = "CD4+_T-cells",
                                      `2` = "CD8+_T-cells",
                                      `3` = "CD4+_T-cells", 
                                      `4` = "Macrophages", 
                                      `5` = "Monocytes", 
                                      `6` = "Tregs", 
                                      `7` = "Gamma-delta_T-cells", 
                                      `8` = "Basophils",
                                      `9` = "NK_cells", 
                                      `10` = "GMP",
                                      `11` = "B_cells",
                                      `12` = "DC")

human_data_integrated$celltypes.author.annotation <- Idents(human_data_integrated)


### Visualization of Seurat object
DimPlot(human_data_integrated, label = TRUE)

#### Plotting cell type specific markers that refer to annotated Seurat clusters

markers.to.plot <- c("CCL17","LAMP3","IRF4","MS4A1","BANK1","S100B","NAPSB","CTSW","TRDC","GATA2","TCN1","IL1RL1","TOP2A","MKI67","FOXP3","IKZF2","CD93","S100A9","S100A8","MMP9","CLEC10A","C1QA","CD8A","CD8B","TSHZ2")
DotPlot(human_data_integrated, features = markers.to.plot, cols = c("red", "green"), dot.scale = 8) & coord_flip() & RotatedAxis()

#### Save seurat object
saveRDS(human_data_integrated,file="My_directory")


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


## Import processed single cell dataset from Oxazolone treated mice containing cells from human PBMC donors as Seurat object

human_data_integrated <- readRDS("Read_annotated_Seurat_File_for_PBMC_data")


## Visualization of UMAP 
DimPlot(human_data_integrated, label = TRUE)


# CD200R surface expression across celltypes annotated by authors
VlnPlot(human_data_integrated, features = c("Human-CD200R"), group.by = "celltypes.author.annotation",pt.size = 0, combine = FALSE,assay = "ADT")

# CD200R surface expression across different treatments for cell types annotated by authors
VlnPlot(human_data_integrated, features = c("Human-CD200R"), group.by = "celltypes.author.annotation",split.by = "stim",pt.size = 0, combine = FALSE,assay = "ADT")


#########################################
# Differential Gene Expression analysis #
#########################################

#Because we have identified common cell types across condition, we can ask what genes change in different conditions for cells of the same type. First, we create a column in #the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column. Then we use FindMarkers() to find the genes that are different between #LY stimulated and control  cells. N

human_data_integrated$celltype.stim <- paste(human_data_integrated$celltypes.author.annotation, human_data_integrated$stim, sep = "_")
Idents(human_data_integrated) <- "celltype.stim"


diff.expression.wilcoxon <- function(treatment1,treatment2) {
  
  # Define data frame to store celltype and gene expression information  
  celltype.de.wilcox<-data.frame()
  
  for (i in unique((human_data_integrated$celltypes.author.annotation))){
    cell.response <- FindMarkers(human_data_integrated, ident.1 = paste(i,treatment1,sep = "_"), ident.2 = paste(i,treatment2,sep = "_"), verbose = FALSE,logfc.threshold = 0)
    cell.response$abslogFC<-abs(cell.response$avg_log2FC)
    cell.response$Gene<-rownames(cell.response)
    cell.response$Celltype<-i
    cell.response$Comparison<-paste(treatment1,treatment2,sep="_")
    celltype.de.wilcox<-rbind(celltype.de.wilcox,cell.response)
  }
  return(celltype.de.wilcox)
}

# Perform differential expression analysis using FindMarker function which uses a Wilcoxon rank sum test to compare cells labeled with CD200R vs. Isotype
CD200R.cells.summary<-diff.expression.wilcoxon("Ucenprubart","Isotype")
colnames(CD200R.cells.summary)<-c("p_val_Dupi","avg_log2FC_CD200R","pct.1_CD200R","pct.2_CD200R","p_val_adj_CD200R","abslogFC_CD200R","Gene","Celltype_CD200R","Comparison")

# Perform differential expression analysis using FindMarker function which uses a Wilcoxon rank sum test to compare cells labeled with Dupilumab vs. Isotype
Dupi.cells.summary<-diff.expression.wilcoxon("Dupilumab","Isotype")
colnames(Dupi.cells.summary)<-c("p_val_Dupi","avg_log2FC_Dupi","pct.1_Dupi","pct.2_Dupi","p_val_adj_Dupi","abslogFC_Dupi","Gene","Celltype_Dupi","Comparison")


#########################################################################
### Overview of differentially expressed genes across cells types
### after Dupilumab and Ucenprubart stimulated cells
### Figure Main Text
#########################################################################

# Barplots of differentially expressed genes using ggplot2

# Define differential expressed genes after CD200R/Ucenprubart stimulation
de.Genes.CD200R<-CD200R.cells.summary %>%filter (p_val_adj_CD200R < 0.05) %>% 
  dplyr::count(Celltype_CD200R) %>% 
  dplyr::rename("Celltype" = "Celltype_CD200R")

de.Genes.CD200R$TRT<-"Ucenprubart vs. Isotype"

# Define differential expressed genes after Dupilumab stimulation
de.Genes.Dupi<-Dupi.cells.summary %>%filter (p_val_adj_Dupi < 0.05)  %>% 
  dplyr::count(Celltype_Dupi) %>% 
  dplyr::rename("Celltype" = "Celltype_Dupi")

de.Genes.Dupi$TRT<-"Dupilumab vs. Isotype"

# Combine differentially expressed genes for Dupilumab and Ucenprubart stimulated cells
de.genes<-rbind(de.Genes.CD200R,de.Genes.Dupi)

level_order <- c( 'CD4+_T-cells','CD8+_T-cells','Gamma-delta_T-cells','Tregs','NK_cells','Macrophages','Monocytes','DC','GMP','B_cells','Basophils') #this vector might be useful for other plots/analyses


## Use color scheme

cvi_cols = c("#D5D2CA","#00A1DE")

# Plot barplot using color scheme with number of differential expressed gene per treatment and celltype
ggplot(data=de.genes, aes(fill=TRT, y=n, x=factor(Celltype,level=level_order))) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values= cvi_cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("") + ylab("Differentially expressed genes, No")


##############################################################
#### Proportional Venn diagram of differentially expressed  
#### genes in macrophages for Figure in main text
##############################################################


# Select differentially expressed genes in macrophages treated with Ucenprubart
d.e.macrophage.genes.CD200R<-CD200R.cells.summary %>% 
  filter (Celltype_CD200R == "Macrophages",p_val_adj_CD200R < 0.05)  %>% 
  dplyr::select (Gene)

# Select differentially expressed genes in macrophages treated with Dupilumab
d.e.macrophage.genes.Dupi<-all.CD4.Tcell.genes.Dupi<-Dupi.cells.summary %>% 
  filter (Celltype_Dupi == "Macrophages", p_val_adj_Dupi < 0.05) %>% 
  dplyr::select (Gene)

# Create matrix from euler package to draw the Venn diagram
fit<-euler(c( "Dupilumab"= length(setdiff(d.e.macrophage.genes.Dupi$Gene,d.e.macrophage.genes.CD200R$Gene)),
              "Ucenprubart" = length(setdiff(d.e.macrophage.genes.CD200R$Gene,d.e.macrophage.genes.Dupi$Gene)),
              "Dupilumab&Ucenprubart"= length(intersect(d.e.macrophage.genes.Dupi$Gene,d.e.macrophage.genes.CD200R$Gene))
))

# Plot proportional Venn Diagram   
plot(fit,fills=cvi_cols,quantities = TRUE,labels = list(font = 15),shape = "ellipse")


########################################################################
# KEGG pathway analysis
#######################################################################

# Select differentially expressed genes in macrophages treated with Ucenprubart
eg.CD200R = bitr(d.e.macrophage.genes.CD200R$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Select differentially expressed genes in macrophages treated with Dupilumab
eg.Dupi = bitr(d.e.macrophage.genes.Dupi$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Select background dataset for macropghages
all.macrophage.genes.CD200R<-CD200R.cells.summary %>% filter (Celltype_CD200R == "Macrophages") %>% dplyr::select (Gene)
all.macrophage.genes.Dupi<-Dupi.cells.summary %>% filter (Celltype_Dupi == "Macrophages") %>% dplyr::select (Gene)
macrophage.universe<-unique(c(all.macrophage.genes.Dupi$Gene,all.macrophage.genes.CD200R$Gene))
universe = bitr(macrophage.universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Create list for clusterprofiles mapping and KEGG enrichments
macrophage.list<-list(Ucenprubart=eg.CD200R$ENTREZID,Dupilumab=eg.Dupi$ENTREZID)

# Perform enrichment analysis in KEGG
ck<-compareCluster(geneCluster = macrophage.list, fun = enrichKEGG,universe=universe$ENTREZID)

# Plot pathway results for differentially expressed genes in macrophages
dotplot(ck,showCategory=10)

##########################################################
### Gene Set Enrichment Analysis
### Pathway enrichment for hallmark genesets from MSiGDB GSEA
### Pathway Analysis using fGSEA
### Figure main text
#########################################################
library(pheatmap)
library(fgsea)

# Define Pathway Dataset from MsigDb (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2)
# Hallmark pathways were obtained from the MSIGdB 2024 release: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/
# See file: "h.all.v2024.1.Hs.symbols.gmt.txt" that containts geneset annotations for gene symbols

pathways <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt.txt")

# Define dataset to store GSEA enrichment information
Celltypes.GSEA<-data.frame()
cellytpes<-unique(human_data_integrated$celltypes.author.annotation)

for (i in cellytpes){
  
  # Prepare Dupi Data set 
  cell.type.enrichment<-Dupi.cells.summary[Dupi.cells.summary$Celltype_Dupi == i,]
  x<-cell.type.enrichment$avg_log2FC_Dupi
  names(x)<-cell.type.enrichment$Gene
  x<-x[order(x)]
  dupi.enrichment<-x[!duplicated(names(x))]
  
  # Prepare CD200R Data set 
  cell.type.enrichment<-CD200R.cells.summary[CD200R.cells.summary$Celltype_CD200R == i,]
  x<-as.numeric(cell.type.enrichment$avg_log2FC_CD200R)
  names(x)<-cell.type.enrichment$Gene
  x<-x[order(x)]
  CD200R.enrichment<-x[!duplicated(names(x))]
  
  # Run GSEA enrichment analysis
  fgseaRes.CD200R<-fgseaRes <- fgsea(pathways = pathways, stats    = CD200R.enrichment,minSize  = 15, maxSize  = 500,eps= 0.0)
  fgseaRes.Dupi<-fgseaRes <- fgsea(pathways = pathways, stats    = dupi.enrichment,minSize  = 15, maxSize  = 500,eps= 0.0)
  fgseaRes.CD200R$Celltype<-paste(i,"CD200R",sep = "_")
  fgseaRes.CD200R$Treatment<-paste("CD200R")
  fgseaRes.CD200R$Cells<-i
  fgseaRes.Dupi$Celltype<-paste(i,"Dupilumab",sep = "_")
  fgseaRes.Dupi$Treatment<-paste("Dupilumab")
  fgseaRes.Dupi$Cells<-i
  x<-rbind(fgseaRes.CD200R,fgseaRes.Dupi)
  Celltypes.GSEA<-rbind(Celltypes.GSEA,x)
}

## Prepare to plot Overview Pathways in a heatmap
Celltypes.GSEA<-data.frame(Celltypes.GSEA[,-"leadingEdge"])
GSEA.Plot<-Celltypes.GSEA[Celltypes.GSEA$padj < 0.1,]

# Generate Plotting matrix for fGSEA heatmap     
pathway.plot<-matrix(nrow=length(na.omit(unique(GSEA.Plot$pathway))),ncol=length(unique(Celltypes.GSEA$Celltype)))
colnames(pathway.plot)<-unique(Celltypes.GSEA$Celltype)
rownames(pathway.plot)<-na.omit(unique(GSEA.Plot$pathway))

# Fill pathway plot with Normalized Enrichment Score values
for (i in 1:length(colnames(pathway.plot))){
  pathway.plot[GSEA.Plot[GSEA.Plot$Celltype %in% colnames(pathway.plot)[i],"pathway"],i]<-GSEA.Plot[GSEA.Plot$Celltype %in% colnames(pathway.plot)[i],"NES"]
}
is.na(pathway.plot)<-sapply(pathway.plot, is.infinite)
pathway.plot[pathway.plot == 0] <- NA

# Sort columns to match wit other plots
pathway.plot<-pathway.plot[,c(sort(colnames(pathway.plot)[grep("CD200R",colnames(pathway.plot))]),sort(colnames(pathway.plot)[grep("Dupi",colnames(pathway.plot))]))]

# Plot pathways using heatmap
pheatmap(pathway.plot,border_color = NA,cluster_rows = F,cluster_cols = F,labels_col = colnames(pathway.plot),na_col = "lightgrey",treeheight_row = 0)





######################################################################
### Correlation Plots for average log2 fold changes across treatments#
### Supplemental Figure
######################################################################


# Create a matrix of data frames harboring log2 fold change correlations for
# cells treated with Ucenprubart or Dupilumab
corr.plots.cells<-list()

for (i in unique(Dupi.cells.summary$Celltype_Dupi)){
  x<-Dupi.cells.summary[Dupi.cells.summary$Celltype_Dupi == i,c("avg_log2FC_Dupi","Gene")]
  colnames(x)[1]<-paste(i,"Dupi",sep="_")
  y<-CD200R.cells.summary[CD200R.cells.summary$Celltype_CD200R == i,c("avg_log2FC_CD200R","Gene")]
  colnames(y)[1]<-paste(i,"CD200R",sep="_")
  correlation<-merge(x,y,by="Gene")
  corr.plots.cells[[i]]<-correlation
}

#merge all data frames in list
corr.plots<-Reduce(function(x, y) merge(x, y, all=TRUE), corr.plots.cells)
corr.plots<-na.omit(corr.plots)
row.names(corr.plots)<-corr.plots[,1]
M<-cor(as.matrix(corr.plots[,-1]))

# Perform correlation analysis
cor.matrix2<-corr.test(as.matrix(corr.plots[,-1]),method="pearson")


# Order Celltypes
level_order <- c( 'CD4+_T-cells','CD8+_T-cells','Gamma-delta_T-cells','Tregs','NK_cells','Macrophages','Monocytes','DC','GMP','B_cells','Basophils')

# use this matrix for corrplot
correlation.matrix<-cor.matrix2$r
correlation.matrix<-correlation.matrix[rownames(correlation.matrix)[grep("CD200R",rownames(correlation.matrix))],colnames(correlation.matrix)[grep("Dupi",colnames(correlation.matrix))]]
correlation.matrix<-correlation.matrix[sort(rownames(correlation.matrix)),sort(colnames(correlation.matrix))]
correlation.matrix<-correlation.matrix[paste(level_order,"CD200R",sep="_"),paste(level_order,"Dupi",sep="_")]

# Add celltypes to correlation matrix
rownames(correlation.matrix)<-gsub("_CD200R","",rownames(correlation.matrix))
colnames(correlation.matrix)<-gsub("_Dupi","",colnames(correlation.matrix))

M<-corrplot(correlation.matrix ,method="number", sig.level = c( 0.05), insig = 'label_sig',tl.col='black')

###############################################################
# KEGG pathway enrichment CD4 T-cells for Dupilumab treatment
# Supplemental Figure
#################################################################
library(clusterProfiler)
library("org.Hs.eg.db")

# KEGG pathway enrichment in CD4 T-cells for Cells stimulated with Dupilumab
all.CD4.Tcell.genes.Dupi<-Dupi.cells.summary %>% 
  filter (Celltype_Dupi == "CD4+_T-cells", p_val_adj_Dupi < 0.05) %>% 
  dplyr::select (Gene)
eg.Dup.CD4.Tcells = bitr(all.CD4.Tcell.genes.Dupi$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Identify significantly enriched pathways in CD4+ T cells stimulated with Dupilumab
kk.Dupi <- enrichKEGG(gene         = eg.Dup.CD4.Tcells$ENTREZID,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)

CD4.t.cell.pathways.Dupi<-data.frame(kk.Dupi)
CD4.t.cell.pathways.Dupi$Treatment<-"Dupilumab"

# KEGG pathway enrichment in CD4 T-cells for Cells stimulated with Ucenprubart
all.CD4.Tcell.genes.CD200R<-CD200R.cells.summary %>% 
  filter (Celltype_CD200R == "CD4+_T-cells",p_val_adj_CD200R < 0.05)  %>% 
  dplyr::select (Gene)
eg.CD200R.CD4.Tcells = bitr(all.CD4.Tcell.genes.CD200R$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Identify  enriched pathways in CD4+ T cells stimulated with Ucenprubart
kk.CD200R <- enrichKEGG(gene       = eg.CD200R.CD4.Tcells$ENTREZID,
                        organism     = 'hsa',
                        pvalueCutoff = 1)

CD4.t.cell.pathways.CD200R<-data.frame(kk.CD200R)
CD4.t.cell.pathways.CD200R$Treatment<-"Ucenprubart"


### Combine Pathway annotations and plot pathways as barplot
pathways.not.significant.in.CD4.Ucenprubart.trt<-setdiff(CD4.t.cell.pathways.Dupi$Description,CD4.t.cell.pathways.CD200R$Description)
# Include mock table for pathways that are not enriched with Ucenprubart treatment
pathways.not.significant.in.CD4.Ucenprubart.trt<-data.frame(category=0,subcategory=0,ID=0,Description=pathways.not.significant.in.CD4.Ucenprubart.trt,GeneRatio=0,BgRatio=0,pvalue=0.99,p.adjust=0.99,qvalue=0.99,geneID=0,Count=0,Treatment="Ucenprubart")
combined.pathways<-rbind(rbind(CD4.t.cell.pathways.CD200R,CD4.t.cell.pathways.Dupi),pathways.not.significant.in.CD4.Ucenprubart.trt)

# Select pathways that are significantly enriched after Dupilumab treatment
pathways.significant.for.Dupi<- combined.pathways %>% filter(qvalue < 0.05 & Treatment == "Dupilumab") %>% dplyr::select (Description)
pathways.significant.for.Dupi<-combined.pathways[combined.pathways$Description %in% pathways.significant.for.Dupi$Description,]

# Plot pathway enrichment as barplot
ggplot(pathways.significant.for.Dupi, aes(fill=Treatment, y=-log10(qvalue), x=Description)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_classic() +
  labs(title="KEGG pathways: CD4 T cells",y="-log10(FDR p-value)", x = "Pathways") +
  geom_hline(yintercept=1.3)+
  theme(legend.position="right")+
  scale_fill_manual(values= cvi_cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



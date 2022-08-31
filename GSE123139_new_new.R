library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
a<-list.files("c:/Users/xjmik/Downloads/GSE123139/GSE123139/GSE123139_RAW/")
setwd("c:/Users/xjmik/Downloads/GSE123139/GSE123139/GSE123139_RAW/")
c<-read.table("c:/Users/xjmik/Downloads/GSE123139/GSE123139/GSE123139.txt",sep = "\t",header = TRUE)
pbmc.data<-read.table(a[1],sep = "\t",header = TRUE,row.names = 1)
pbmc.metadata<-data.frame(colnames(pbmc.data),rep(c[1,2],length(colnames(pbmc.data))))
colnames(pbmc.metadata)<-c("Barcode","Group")
for (i in 2:length(a)) {
  pbmc_new.data<-read.table(a[i],sep = "\t",header = TRUE,row.names = 1)
  pbmc_new.metadata<-data.frame(colnames(pbmc_new.data),rep(c[i,2],length(colnames(pbmc_new.data))))
  colnames(pbmc_new.metadata)<-c("Barcode","Group")
  b<-intersect(rownames(pbmc.data),rownames(pbmc_new.data))
  pbmc.data<-pbmc.data[b,]
  pbmc_new.data<-pbmc_new.data[b,]
  remove(b)
  pbmc.data<-cbind(pbmc.data,pbmc_new.data)
  pbmc.metadata<-rbind(pbmc.metadata,pbmc_new.metadata)
}
remove(pbmc_new.data,pbmc_new.metadata,i,a,c)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.metadata<-pbmc.metadata[colnames(pbmc.data),]

pbmc.metadata_Tumor<-pbmc.metadata[which(as.character(pbmc.metadata$Group) ==" Tumor"),]
pbmc.metadata_PBMC<-pbmc.metadata[which(as.character(pbmc.metadata$Group) == " PBMC"),]
pbmc_Tumor.data<-pbmc.data[,rownames(pbmc.metadata_Tumor)]
pbmc_PBMC.data<-pbmc.data[,rownames(pbmc.metadata_PBMC)]

pbmc_Tumor <- CreateSeuratObject(counts = pbmc_Tumor.data, project = "IMMUNE_pbmc_Tumor",meta.data = pbmc.metadata_Tumor,min.cells = 3,min.features = 200)
pbmc_Tumor$type <- "pbmc_Tumor"

pbmc_Tumor <- NormalizeData(pbmc_Tumor, verbose = FALSE)
pbmc_Tumor <- FindVariableFeatures(pbmc_Tumor, selection.method = "vst", nfeatures = 2000)

pbmc_PBMC <- CreateSeuratObject(counts = pbmc_PBMC.data, project = "IMMUNE_pbmc_PBMC",meta.data = pbmc.metadata_PBMC,min.cells = 3,min.features = 200)
pbmc_PBMC$type <- "pbmc_PBMC"

pbmc_PBMC <- NormalizeData(pbmc_PBMC, verbose = FALSE)
pbmc_PBMC <- FindVariableFeatures(pbmc_PBMC, selection.method = "vst", nfeatures = 2000)

pbmc_new1.data<-Read10X("c:/Users/xjmik/Downloads/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")
pbmc_new2.data<-Read10X("c:/Users/xjmik/Downloads/frozen_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")
pbmc_new3.data<-Read10X("c:/Users/xjmik/Downloads/frozen_pbmc_donor_b_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")
pbmc_new4.data<-Read10X("c:/Users/xjmik/Downloads/frozen_pbmc_donor_c_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/")

for (i in 1:length(colnames(pbmc_new1.data))) {
  colnames(pbmc_new1.data)[i] <- paste(colnames(pbmc_new1.data)[i],"pbmc_new1",i,sep = "-")  
}

for (i in 1:length(colnames(pbmc_new2.data))) {
  colnames(pbmc_new2.data)[i] <- paste(colnames(pbmc_new2.data)[i],"pbmc_new2",i,sep = "-")  
}

for (i in 1:length(colnames(pbmc_new3.data))) {
  colnames(pbmc_new3.data)[i] <- paste(colnames(pbmc_new3.data)[i],"pbmc_new3",i,sep = "-")  
}

for (i in 1:length(colnames(pbmc_new4.data))) {
  colnames(pbmc_new4.data)[i] <- paste(colnames(pbmc_new4.data)[i],"pbmc_new4",i,sep = "-")  
}

pbmc_new.data<-cbind(pbmc_new1.data,pbmc_new2.data,pbmc_new3.data,pbmc_new4.data)

pbmc_new <- CreateSeuratObject(counts = pbmc_new.data, project = "IMMUNE_pbmc_new",min.cells = 3,min.features = 200)
pbmc_new$type <- "pbmc_new"

pbmc_new <- NormalizeData(pbmc_new, verbose = FALSE)
pbmc_new <- FindVariableFeatures(pbmc_new, selection.method = "vst", nfeatures = 2000)


immune.anchors <- FindIntegrationAnchors(object.list = list(pbmc_PBMC,pbmc_Tumor,pbmc_new), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

FeaturePlot(immune.combined,features = "FOXP3")
VlnPlot(immune.combined,features = "FOXP3",sort = TRUE,pt.size = 0,assay = "RNA")
Treg<-subset(immune.combined,idents = "1")
Idents(Treg)<-Treg@meta.data$Group
library(ggpubr)
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()
FeaturePlot(Treg,features = "JMJD1C",cols = c("#00008B","#FF69B4"),min.cutoff = "q10",max.cutoff = "q95",pt.size = 2,order = TRUE,reduction = "umap",split.by = "Group")
Treg<-ScaleData(Treg)
Treg <- RunPCA(Treg, npcs = 30, verbose = FALSE)

Treg <- RunUMAP(Treg, reduction = "pca", dims = 1:20)
Treg <- FindNeighbors(Treg, reduction = "pca", dims = 1:20)
Treg <- FindClusters(Treg, resolution = 1)
FeaturePlot(Treg,features = "FOXP3")
DimPlot(Treg, reduction = "umap", label = TRUE)

VlnPlot(Treg,features = "FOXP3",sort = TRUE,pt.size = 0,assay = "RNA")
Treg_new<-subset(Treg,idents=c("6"))
Idents(Treg_new)<-Treg_new@meta.data$type
VlnPlot(Treg_new,features = "JMJD1C",sort = TRUE,pt.size = 0)
library(ggpubr)
VlnPlot(Treg_new,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()
Treg_new_new<-subset(Treg_new, ident=c("pbmc_Tumor","pbmc_new"))
VlnPlot(Treg_new_new,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()

FeaturePlot(Treg_new,features = "JMJD1C",split.by = "Group",order = TRUE)
JMJD1C<-FetchData(Treg_new_new,vars = "JMJD1C")
JMJD1C_metadata<-data.frame(rownames(Treg_new_new@meta.data),Treg_new_new@meta.data$type)
JMJD1C_new<-cbind(JMJD1C,JMJD1C_metadata[,2])
write.table(JMJD1C_new,file = "GSE121638.txt",sep = "\t")

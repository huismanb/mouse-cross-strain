# Three-strain mouse data: preprocess, assign HTOs, QC

library(tidyverse)
library(Seurat)
library(Matrix)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(future)
library(ggplot2)

plan("sequential")
date=""

cmap <- c("#D7B5A6","#4E79A7","#D37295","#A0CBE8",
          "#F28E2B","#59A14F","#FFBE7D","#8CD17D","#B6992D","#F1CE63",
          "#D4A6C8","#86BCB6","#E15759","#FF9D9A","#499894","#BAB0AC",
          "#FABFD2","#B07AA1","#9D7660")

#############
# Initialize the Seurat object with the raw (non-normalized data)

initial_fn <- function(dat,proj_name){
  hto_adt_remove_from_analysis<-c("HT4","HT7","HT10","HT11","HT13","HT14","HT15")
  hashtag<-dat[[2]]
  hashtag <- hashtag[!rownames(hashtag) %in% hto_adt_remove_from_analysis,]
  
  gene<-dat[[1]]
  seurat_table <- CreateSeuratObject(counts = gene, project = proj_name, min.cells=3, min.features = 0)
  
  seurat_table[['HTO']] <- CreateAssayObject(counts=hashtag)
  seurat_table <- NormalizeData(seurat_table, assay = "HTO", normalization.method = "CLR")
  seurat_table <- ScaleData(seurat_table, assay = "HTO")
  seurat_table <- HTODemux(seurat_table, assay = "HTO", positive.quantile = 0.99)
  
  return(seurat_table)
}

#custom HTO cutoffs
custom_hto_assignnment <- function(seurat_table, cutoff_list){
  keep <- vector(mode="logical", length=ncol(seurat_table)) #initialize empty vector
  for( i in 1:ncol(seurat_table) ) {
    a <- seurat_table@assays$HTO@data[1,i] > cutoff_list[1] 
    b <- seurat_table@assays$HTO@data[2,i] > cutoff_list[2]
    c <- seurat_table@assays$HTO@data[3,i] > cutoff_list[3]
    d <- seurat_table@assays$HTO@data[4,i] > cutoff_list[4]
    e <- seurat_table@assays$HTO@data[5,i] > cutoff_list[5]
    f <- seurat_table@assays$HTO@data[6,i] > cutoff_list[6]
    g <- seurat_table@assays$HTO@data[7,i] > cutoff_list[7]
    h <- seurat_table@assays$HTO@data[8,i] > cutoff_list[8]
    if ( sum(a+b+c+d+e+f+g+h)!=1 ) {
      keep[i] = F
    } else {
      keep[i] = T
    }
  }
  seurat_table <- seurat_table[,keep]
  
  #Assign HTO (hash.ident) to cells if above cutoff
  ident <- vector(mode="logical", length=ncol(seurat_table))
  for( i in 1:ncol(seurat_table) ) {
    if( seurat_table@assays$HTO@data[1,i] > cutoff_list[1] ) {
      ident[i] <- "b61"
    } else if ( seurat_table@assays$HTO@data[2,i] > cutoff_list[2] ) {
      ident[i] <- "b62"
    } else if ( seurat_table@assays$HTO@data[3,i] > cutoff_list[3] ) {
      ident[i] <- "balbc1"
    } else if ( seurat_table@assays$HTO@data[4,i] > cutoff_list[4] ) {
      ident[i] <- "balbc2"
    } else if ( seurat_table@assays$HTO@data[5,i] > cutoff_list[5] ) {
      ident[i] <- "Nod1"
    } else if ( seurat_table@assays$HTO@data[6,i] > cutoff_list[6] ) {
      ident[i] <- "Nod2"
    } else if ( seurat_table@assays$HTO@data[7,i] > cutoff_list[7] ) {
      ident[i] <- "test1"
    } else if ( seurat_table@assays$HTO@data[8,i] > cutoff_list[8] ) {
      ident[i] <- "test2"
    }
    else {
      ident[i] <- "na"
    }
  }
  names(ident) <- colnames(seurat_table)
  seurat_table@meta.data$hash.ident <- factor(ident)
  
  seurat_table$geno=''
  seurat_table$geno[seurat_table@meta.data[['hash.ident']] %in% c('b61','b62')] <- 'B6'
  seurat_table$geno[seurat_table@meta.data[['hash.ident']] %in% c('balbc1','balbc2')] <- 'BALBC'
  seurat_table$geno[seurat_table@meta.data[['hash.ident']] %in% c('Nod1','Nod2')] <- 'NOD'
  seurat_table$geno[seurat_table@meta.data[['hash.ident']] %in% c('test1','test2')] <- 'test'
  seurat_table$geno <- factor(seurat_table$geno, levels=c("B6","BALBC","NOD","test"))
  
  return(seurat_table)
}

dat = Read10X(data.dir='/DN/')
seurat_table_dn = initial_fn(dat,'dn')

RidgePlot(seurat_table_dn, assay = "HTO", features = rownames(seurat_table_dn[["HTO"]])[1:8], ncol = 4)

#HTO cutoffs:
cutoff_dn = c(2, 2, 1.5 , 1.5 , 2, 2, 1.5 , 1.5)
seurat_table_dn <- custom_hto_assignnment(seurat_table_dn, cutoff_dn)

#Show ridgeplot once added custom cutoffs
Idents(seurat_table_dn)<-"hash.ident"

RidgePlot(seurat_table_dn, assay = "HTO", features = rownames(seurat_table_dn[["HTO"]])[1:8], ncol = 4)

seurat_table=seurat_table_dn

##############################
###### filter data
seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^mt-")

#QC
x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=9000, color="red" ) + #line where cutoff is
  geom_hline( yintercept=300, color="red" ) + #line where cutoff is
  geom_violin( fill="steelblue" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  geom_violin( fill="pink" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=12.5, color="red" ) +
  geom_violin( fill="orchid4" ) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape=NA) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_prefilter_",date,".pdf" ), height=3, width=8 )
plot_grid( x,y,z, ncol=3 )
dev.off()
##### 
seurat_table <- subset(seurat_table, subset = percent.mt < 12.5)
seurat_table <- subset(seurat_table, subset = nFeature_RNA > 300) 
seurat_table <- subset(seurat_table, subset = nFeature_RNA < 9000)
#####
#normalize data
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf( paste0("figures/variable_feature_",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()


#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

sink( paste0("text_outputs/top_PC_genes_1-30_",date,".txt"), append=F )
print(seurat_table[["pca"]], dims = 1:30, nfeatures = 5)
sink()

pdf( paste0("figures/PC_dim_loading_",date,".pdf" ), height=40, width=12 )
VizDimLoadings(seurat_table, dims = 1:30, reduction = "pca")
dev.off()

pdf( paste0("figures/pca_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction = "pca", group.by = "orig.ident")
dev.off()

pdf( paste0("figures/PC_dim_heatmap_",date,".pdf" ), height=40, width=12 )
DimHeatmap(seurat_table, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

pdf( paste0("figures/PC_elbow_postfilter_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap
seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
seurat_table <- FindClusters(seurat_table, resolution = 2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)

pdf( paste0("figures/umap_prefiltering_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

pdf( paste0("figures/umap_prefiltering_origident_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1,group.by='orig.ident')
dev.off()

pdf( paste0("figures/umap_prefiltering_geno_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1,group.by='geno')
dev.off()

saveRDS(seurat_table, file = "Rdata/cross-strain_seurat-table_prefilter.rds")

##### 
seurat_table <- JoinLayers(seurat_table)
#####
#Determine identities of cells then remove contaminating cells 
keep <- seurat_table@active.ident!=26 & seurat_table@active.ident!=25 & 
  seurat_table@active.ident!=14 & seurat_table@active.ident!=19 & seurat_table@active.ident!=28
seurat_table <- seurat_table[,keep]

#Remove subsets that are not of interest
Idents(seurat_table)<-"geno"
seurat_table <- subset(seurat_table,idents = c("B6","BALBC","NOD"))

Idents(seurat_table)<-"orig.ident"
seurat_table <- subset(seurat_table,idents = c("dn"))

############### Repeat normalization, find var features, etc
#normalize data
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)
top10 <- head(VariableFeatures(seurat_table), 10)

#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

#####
# cluster and umap
seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
seurat_table <- FindClusters(seurat_table, resolution = 2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)

saveRDS(seurat_table, file = "Rdata/cross-strain_seurat-table.rds")

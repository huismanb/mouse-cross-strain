# Three-strain mouse data: assign labels, analyze, generate figures

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

cmap <- c("#D7B5A6",
          "#4E79A7","#D37295","#A0CBE8",
          "#F28E2B","#59A14F","#FFBE7D","#8CD17D","#B6992D","#F1CE63",
          "#D4A6C8","#86BCB6","#E15759","#FF9D9A","#499894","#BAB0AC",
          "#FABFD2","#B07AA1","#9D7660"
)

#############
seurat_table <- readRDS( "Rdata/cross-strain_seurat-table.rds")

############
# Heatmap
markers <- FindAllMarkers(seurat_table, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10
heatmap_features = top10$gene

breakslist = seq(1,3,by=0.01)
pdf( paste0("figures/markers_heatmap_purple_dn-only_3strain_",date,".pdf" ), height=15, width=7 )
DoHeatmap(subset(seurat_table, downsample=100), features = heatmap_features,angle=90, disp.min=0.5, size=4) + theme(axis.text.y.left = element_text(face = "italic"))+ 
  scale_fill_gradientn(colors = colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist))) 
dev.off()

## DotPlot
pdf( paste0("figures/dotplot_TFs_dn-only_3strain_",date,".pdf"), height=6, width=10 )
DotPlot(seurat_table, features=(c('Myog','Neurod1','Aire',
                                  "Grhl1",
                                  "Foxi1",
                                  "Pou2f3",
                                  "Sox8","Spib","Hnf4g","Hnf4a","Foxj1",
                                  "Foxn1","Trp63","Trp73","Foxi2","Ptf1a","Foxa2","Foxa3","Foxa1","Spdef","Gata3"
)),
cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 ) + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(face = "italic"))+
  scale_y_discrete(limits=rev) #reverse y-direction
dev.off()

#########
#Label clusters
# increase # of clusters called for improved granularity
seurat_table <- FindClusters(seurat_table, resolution = 3)
#seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)

DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1, cols=cmap)

keep <- seurat_table@active.ident==6 | seurat_table@active.ident==11 | seurat_table@active.ident==12 | seurat_table@active.ident==3  | 
  seurat_table@active.ident==8 | seurat_table@active.ident==14 
seurat_table@active.ident[keep] <- 3 

seurat_table <- RenameIdents( object=seurat_table,
                              "0"="Ciliated",
                              "2"="Goblet/Lung basal",
                              "10"="Enterohetapto",
                              "13"="MCell",
                              "16"="Muscle",
                              "4"="Neuroendocrine",
                              "7"="Skin basal",
                              "5"="Skin keratinized",
                              "3"="Tuft (neuro)",
                              "1"="Tuft (sensory)",
                              "9"="Tuft (neuro)", 
                              "15"="Tuft (neuro)" 
                              
)

seurat_table@active.ident <- factor(seurat_table@active.ident,
                                    levels=c('Muscle', 'Tuft (neuro)',
                                             'Tuft (sensory)','MCell', 'Neuroendocrine',
                                             'Enterohetapto',
                                             'Goblet/Lung basal', 'Skin basal', 'Skin keratinized',
                                             'Ciliated'
                                    ))

pdf( paste0("figures/umap_split_", date,".pdf"), height=5, width=18 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1, split.by='geno', cols=cmap ) 
dev.off()

pdf( paste0("figures/umap_split_nolabel_", date,".pdf"), height=5, width=18 )
DimPlot(seurat_table, reduction="umap", label=F, repel=T, pt.size=2, split.by='geno', cols=cmap ) & NoLegend() & NoAxes() 
dev.off()

#################################
#plot umap as a density plot
library(MASS)
library(BuenColors)

unique(seurat_table@meta.data$geno) #should be 3 subgroups (B6, BALBC, NOD)

h=3
n=100

for(group in c('B6','BALBC','NOD')){
  density_ingroup = kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$geno %in% c(group),"umap_1"], 
                           y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$geno %in% c(group),"umap_2"],
                           h=h,
                           n=n ) 
  
  xy = xy.coords( x=seurat_table@reductions$umap@cell.embeddings[,"umap_1"],
                  y=seurat_table@reductions$umap@cell.embeddings[,"umap_2"] )
  select = is.finite(xy$x) & is.finite(xy$y)
  x = cbind(xy$x, xy$y)[select, ]
  mkBreaks = function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin = cut(x[, 1], mkBreaks(density_ingroup$x), labels = FALSE)
  ybin = cut(x[, 2], mkBreaks(density_ingroup$y), labels = FALSE)
  dens = density_ingroup$z[cbind(xbin, ybin)]
  
  p = ggplot( as.data.frame( seurat_table@reductions$umap@cell.embeddings ),
              aes( x=seurat_table@reductions$umap@cell.embeddings[,"umap_1"],
                   y=seurat_table@reductions$umap@cell.embeddings[,"umap_2"] ) ) +
    geom_point( aes( color=oob_squish(dens, range=c(0.00,0.03)) ), size=2 ) +
    theme_bw() +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    scale_color_gradientn( colors=jdb_palette("solar_extra"), na.value = "grey50") +
    
    labs( color=paste0("Density: ",group )) +
    xlab("UMAP1") +
    ylab("UMAP2")
  
  p1=p & NoLegend() 
  
  pdf( paste0("figures/density_",group,"_legend_h",h,'_n',n,'_', date,".pdf"), height=5, width=7 )
  print(p)
  dev.off()
  
  pdf( paste0("figures/density_",group,"_nolegend_h",h,'_n',n,'_', date,".pdf"), height=5, width=6 )
  print(p1)
  dev.off()
}


##############################
#Barplot
cluster_dist = t(table(seurat_table@active.ident, seurat_table@meta.data$hash.ident))
cluster_dist
barplot(cluster_dist,beside=TRUE,ylab='Number of cells',las=2, col=rep(cmap, each=6))

cluster_dist_norm = cluster_dist/rowSums(cluster_dist)*100
cluster_dist_norm
barplot(cluster_dist_norm,beside=TRUE,ylab='% of mimetics in mouse',las=2, col=rep(cmap, each=6))

#individual donor calculations
cluster_df = as.data.frame(cluster_dist_norm)
names(cluster_df)[1] = 'donor'
names(cluster_df)[2] = 'celltype'
cluster_df

cluster_df$strain = '_'
cluster_df[cluster_df$donor=='b61','strain'] = 'B6'
cluster_df[cluster_df$donor=='b62','strain'] = 'B6'
cluster_df[cluster_df$donor=='balbc1','strain'] = 'BALBC'
cluster_df[cluster_df$donor=='balbc2','strain'] = 'BALBC'
cluster_df[cluster_df$donor=='Nod1','strain'] = 'NOD'
cluster_df[cluster_df$donor=='Nod2','strain'] = 'NOD'
cluster_df

avg_strain = aggregate(cluster_df$Freq, list(cluster_df$celltype, cluster_df$strain), FUN=mean)
avg_strain

pdf( paste0("figures/cluster_distribution_",date,".pdf" ), height=2.5, width=12 )
ggplot(cluster_df, aes(x = factor(celltype), y = Freq, fill = strain)) +
  geom_col(data=avg_strain, aes(x = factor(Group.1), y = x, fill = factor(Group.2)), color = "black", linewidth = 1, width = .8, position = "dodge") +
  geom_point(position = position_jitterdodge(jitter.width = 0), size=1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("grey", "white",'#494fc2')) + ylab('Percentage of mimetic cells') + 
  theme(axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"))
dev.off()  

##############################
# HEATMAP 
markers <- FindAllMarkers(seurat_table, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10
heatmap_features = top10$gene

breakslist = seq(1,3,by=0.01)
pdf( paste0("figures/markers_heatmap_purple_",date,".pdf" ), height=15, width=7 )
DoHeatmap(subset(seurat_table, downsample=100), features = heatmap_features,angle=90, disp.min=0.5, size=4, group.colors=cmap) + theme(axis.text.y.left = element_text(face = "italic"))+ #NoLegend() + #size changes cluster label size
  scale_fill_gradientn(colors = colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist))) 
dev.off()

#select features
heatmap_features = c('Acta1','Ckm','Myl1','Mymx','Cox8b','Dlk1','Tnnt3',
                     'Nrgn','Sh2d6','Alox5ap','Dclk1','Luzp2','Cntnap5a',
                     'Bex4','Scn2a','Dcc','Tas2r105','Nrxn3','Gnat3',
                     'Gp2','Prg2','Cldn13','Hamp','Ccl20','Ccl9',
                     'Snap25','Cacna2d1','Chga','Chgb','Stxbp5l',
                     'Lypd8','Nos2','Fabp9','Cd70','Igf1',
                     'Rptn','Aqp5','Aqp4','Krt5','Spink5',
                     'Ecm1','Cyp2b19','Lingo2',
                     'Csta1','Ivl','Krtdap','Krt1','Lor','Klk5','Klk6','Them5',
                     'Ccdc153','Sntn','Spag17','Tekt1','Rp1')

breakslist = seq(1,3,by=0.01)
pdf( paste0("figures/markers_heatmap_purple_select_",date,".pdf" ), height=8, width=7 )
DoHeatmap(subset(seurat_table, downsample=50), features = heatmap_features,angle=90, disp.min=0.5, size=4, group.colors=cmap) + theme(axis.text.y.left = element_text(face = "italic"))+ #NoLegend() + #size changes cluster label size
  scale_fill_gradientn(colors = colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist))) 
dev.off()

##############################
## overlay mouse signatures from Michelson et al Cell paper
setwd("/signatures_mouse")
idx <- list.files()[grepl("_sig_", list.files())]
sig_list <- list()
for (i in idx ) {
  sig <- read.delim(i,header=T, sep="\t")
  sig <- sig[,1]
  sig_name <- substr(i, 10, nchar(i)-30)
  sig_list[[sig_name]] <- sig
  seurat_table <- AddModuleScore(seurat_table, features=as.data.frame(sig), name=sig_name)
}

pdf( paste0("figures/overlay_mouse_signatures_",date,".pdf" ), height=15, width=15 )
FeaturePlot(seurat_table, features=c(
  "Tuft11", "Tuft21",
  "Muscle1","Ciliated1", 
  "Neuroendocrine1","Ptf1a+ ductal1", "Mcell1", 
  "Gut1", "Goblet1", "Lung, basal1",   
  "Skin, basal1", "Skin, keratinized1",
  "Ionocyte1"
), order=T, min.cutoff="q75",pt.size=1, cols=c("lightgray","#3F007D"))
dev.off()



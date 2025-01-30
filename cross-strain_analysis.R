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

cmap <- c("#D7B5A6","#4E79A7","#D37295","#A0CBE8",
          "#F28E2B","#59A14F","#FFBE7D","#8CD17D","#B6992D","#F1CE63",
          "#D4A6C8","#86BCB6","#E15759","#FF9D9A","#499894","#BAB0AC",
          "#FABFD2","#B07AA1","#9D7660")

#############
seurat_table <- readRDS( "Rdata/cross-strain_seurat-table.rds")

############
#Label clusters
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
#Plot UMAP as a density plot
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

source('Codes/Functions.R')
Initialize()
# can access the clusters this way too: sCVdata_list$res.0.2@Clusters

output_filename <- "refinedFiles/refined_seurat_protze.RData"
load(output_filename)
Idents(seur) <- paste0('cluster_', seur@meta.data$SCT_snn_res.0.2)
table(Idents(seur))

tsne.clust.df <- data.frame(cbind(getEmb(seur, "tsne"), 
                                  clusters=as.character(Idents(seur))))
tsne.clust.df$tSNE_1 <- as.numeric(as.character(tsne.clust.df$tSNE_1))
tsne.clust.df$tSNE_2 <- as.numeric(as.character(tsne.clust.df$tSNE_2))
head(tsne.clust.df)
ggplot(tsne.clust.df, aes(x=tSNE_1, y=tSNE_2, color= clusters))+
  geom_point()+scale_color_manual(values=colorPalatte)+theme_bw()

SA_NODE = 'clusters_2'
TRANSITIONAL = 'clusters_5'


## mit included , expression pattern of genes
seur.mit <- readRDS('Seur_objects/seur_mit_included_clustered.rds')
Idents(seur.mit) <- paste0('cluster_', seur.mit@meta.data$SCT_snn_res.0.2)
GENE_NAME = 'TBX3'
gene.expr <- GetAssayData(object= seur.mit, assay='SCT')[rownames(seur.mit)==GENE_NAME,]
tsne.gene.df <- data.frame(cbind(getEmb(seur.mit,"tsne"), gene.expr))
Plot.tsne.gene.expr(tsne.gene.df, GENE_NAME)



## labeling the Y shape region manually
# x > -5, 10
# y > 5, 25
tsne.clust.df.mit <- data.frame(cbind(getEmb(seur.mit, "tsne"), 
                                      clusters=as.character(Idents(seur.mit))))
tsne.clust.df.mit$tSNE_1 <- as.numeric(as.character(tsne.clust.df.mit$tSNE_1))
tsne.clust.df.mit$tSNE_2 <- as.numeric(as.character(tsne.clust.df.mit$tSNE_2))
head(tsne.clust.df.mit)
ggplot(tsne.clust.df.mit, aes(x=tSNE_1, y=tSNE_2, color= clusters))+
  geom_point()+scale_color_manual(values=colorPalatte)+theme_bw()

tsne.clust.df.mit$label = ifelse(tsne.clust.df.mit$tSNE_1 > (-5) & tsne.clust.df.mit$tSNE_1 < 10 &
                                   tsne.clust.df.mit$tSNE_2 > 5 & tsne.clust.df.mit$tSNE_2 < 25, 
                                 'Y_shape', 'None')
ggplot(tsne.clust.df.mit, aes(x=tSNE_1, y=tSNE_2, color= label))+
  geom_point()+scale_color_manual(values=c('gray', 'red3'))+theme_bw()

## Finding the barcodes of the Y-shaped points
table(tsne.clust.df.mit$label)
Y_shape_barcodes <- rownames(tsne.clust.df.mit[tsne.clust.df.mit$label=='Y_shape',])


## labeling them in the mit-included plot
tsne.clust.df$label <- ifelse(rownames(tsne.clust.df) %in% Y_shape_barcodes, 'Y_shape', 'None')
ggplot(tsne.clust.df, aes(x=tSNE_1, y=tSNE_2, color= label))+
  geom_point()+scale_color_manual(values=c('gray', 'red3'))+theme_bw()


seur.mit.manClust <- seur.mit
Idents(seur.mit.manClust) <- tsne.clust.df.mit$label
Y_shape_markers <- FindMarkers(seur.mit.manClust, ident.1 = 'Y_shape', ident.2 = 'None')
write.csv(Y_shape_markers, '~/Desktop/Y_shape_markers.csv')



## list of needed markers
Markers <- c('TNNT2', 'SIRPA', 'TBX3', 'SHOX2', 'NKX2-5','HCN4',
             "TBX18", 'CD34', 'SCN5A','GJA5','NPPA','MYL2','NFATC1',
             'NRG1','PECAM1','TCF21','VIM','DDR2','THY1','WT1','CSF1R','MYH11')

pdf('~/Desktop/markers_expression_mit.pdf')
for (i in 1:length(Markers)){
  print(i)
  GENE_NAME = Markers[i]
  gene.expr <- GetAssayData(object= seur.mit, assay='SCT')[rownames(seur.mit)==GENE_NAME,]
  tsne.gene.df <- data.frame(cbind(getEmb(seur.mit,"tsne"), gene.expr))
  print(Plot.tsne.gene.expr(tsne.gene.df, GENE_NAME))
}
dev.off()

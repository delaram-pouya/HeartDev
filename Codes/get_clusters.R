
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
INPUT_NAME = args[1]
INPUT_FILE = args[2]

source('Codes/Functions.R')
Initialize()

# INPUT_NAME = 'in_vivo'
# INPUT_FILE = '2.seur_dimRed_in_vivo_mito_0.593_lib_1500.rds'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
PC_NUMBER = 25

ClUSTER_RESOLUTION = 0.5
if(INPUT_NAME == 'in_vivo') {ClUSTER_RESOLUTION = 0.2
}else if(INPUT_NAME == 'in_vitro') ClUSTER_RESOLUTION = 1.2

  
seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
seur <- FindNeighbors(seur, dims = 1:PC_NUMBER)
seur <- FindClusters(seur, resolution = ClUSTER_RESOLUTION)
#seur$seurat_clusters <- paste0('cluster_', seur$seurat_clusters)
table(seur$seurat_clusters)


pdf(paste0('Results/',INPUT_NAME,'/clusters/clusters_',OUTPUT_NAME,'.pdf'))

tsne_df <- getTsneDF(seur)
ggplot(tsne_df,aes(x= tSNE_1, y= tSNE_2, color=clusters))+geom_point()+
  theme_bw()+ggtitle(paste0('tSNE ', OUTPUT_NAME, ' ,res=', ClUSTER_RESOLUTION))

umap_df <- getUmapDF(seur)
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=clusters))+geom_point()+
  theme_bw()+ggtitle(paste0('UMAP ', OUTPUT_NAME, ' ,res=', ClUSTER_RESOLUTION))

pca_df <- getPcaDF(seur)
ggplot(pca_df,aes(x= PC_1, y= PC_2, color=clusters))+geom_point()+
  theme_bw()+ggtitle(paste0('PCA ', OUTPUT_NAME, ' ,res=', ClUSTER_RESOLUTION))

dev.off()



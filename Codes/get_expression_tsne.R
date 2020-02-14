## Run this script as: 
# Rscript Codes/get_expression_tsne.R 'in_vivo' '2.seur_dimRed_in_vivo_mito_0.593_lib_1500.rds'


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()

# INPUT_NAME = 'in_vivo'
# INPUT_FILE = '2.seur_dimRed_in_vivo_mito_0.593_lib_1500.rds'
INPUT_NAME = args[1] 
INPUT_FILE = args[2]
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))


### importing expression matrix and list of markers
candidateGenes_mapped_df <- readRDS('Data/cardiac_markers_ensembl_mapped_table.rds')
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$ensembl_gene_id))
candidateGenes_mapped <- lapply(candidateGenes_mapped, function(x) x[x %in% rownames(seur)])


seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- as.matrix(seur[['RNA']]@data)


pdf(paste0("Results/",INPUT_NAME,'/markers/marker_expression_',OUTPUT_NAME,'.pdf'))

for(cell_index in 1:length(candidateGenes_mapped)){
  for(marker_index in 1:length(candidateGenes_mapped[[cell_index]])){
    
    print(paste0(cell_index,' ' ,marker_index))
    
    aMarker <- candidateGenes_mapped[[cell_index]][marker_index]
    aMarker_ensemble_id_index <- candidateGenes_mapped_df[[cell_index]]$ensembl_gene_id == aMarker
    aMarker_symbol = candidateGenes_mapped_df[[cell_index]]$hgnc_symbol[aMarker_ensemble_id_index]
    aMarker_expression <- exprMatrix[aMarker,]
    
    tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
    
    p= Plot.tsne.gene.expr(tsne_df, paste0(names(candidateGenes_mapped)[cell_index],' - ', aMarker_symbol))
    print(p)
  }
}

dev.off()







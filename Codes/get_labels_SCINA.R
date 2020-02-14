## Run this script as: 
# Rscript Codes/get_labels_SCINA.R 'in_vivo' '2.seur_dimRed_in_vivo_mito_0.593_lib_1500.rds'

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


seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- as.matrix(seur[['RNA']]@data)

### importing expression matrix and list of markers
candidateGenes_mapped_df <- readRDS('Data/cardiac_markers_ensembl_mapped_table.rds')
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$ensembl_gene_id))
candidateGenes_mapped <- lapply(candidateGenes_mapped, function(x) x[x %in% rownames(seur)])



SCINA_res = SCINA(exprMatrix, candidateGenes_mapped, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.99, sensitivity_cutoff = 1, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')

rownames(SCINA_res$probabilities)
table(SCINA_res$cell_labels)
saveRDS(SCINA_res, paste0('Results/',INPUT_NAME,'/SCINA/','SCINA_',OUTPUT_NAME,'.rds'))





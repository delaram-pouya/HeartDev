## Run this script as: 
# Rscript Codes/get_labels_AUCell.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

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
candidateGenes_not_included <- lapply(candidateGenes_mapped, function(x) x[!x %in% rownames(seur)])
candidateGenes_mapped <- lapply(candidateGenes_mapped, function(x) x[x %in% rownames(seur)])



#### Staring the analysis for signiture finding and labeling using AUCell
AUCell_dir = paste0("Results/",INPUT_NAME,"/AUCell/")


gmtFile <- paste0('Data/cardiac_markers_ensembl.gmt')
geneSets <- getGmt(gmtFile)

all_markers <- as.character(unlist(candidateGenes_mapped))
all_markers[!all_markers %in% rownames(seur)]

geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))



## build gene expression ranking for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=detectCores()-1, plotStats=TRUE)
saveRDS(cells_rankings, paste0(AUCell_dir,"cells_rankings_",OUTPUT_NAME,".rds" ))

# Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # the top 4 percent are being considered: can change with: aucMaxRank
saveRDS(cells_AUC, file=paste0(AUCell_dir, "cells_AUC_",OUTPUT_NAME,".rds"))

## converting it to a data frame structure
cells_AUC_df = data.frame(getAUC(cells_AUC))

## Normalizing to consider probabilistic interpretation 
cells_AUC_df_norm <- mapply(`/`, cells_AUC_df,  colSums(cells_AUC_df))  ## normalizing each column to have a probabilistic point of view
rownames(cells_AUC_df_norm) <- rownames(cells_AUC_df)
rough_cell_type_assignment <- rownames(cells_AUC_df_norm)[apply(cells_AUC_df_norm,2,which.max)]
cbind(table(rough_cell_type_assignment)) ## 93% are the hepatocytes!!

## visualzing rough estimation
cell_tsne <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1], tSNE_2 = getEmb(seur, 'tsne')[,2], cell_type=rough_cell_type_assignment)
ggplot(cell_tsne, aes(x=tSNE_1, y=tSNE_2, color=cell_type))+geom_point(alpha=0.5)+theme_bw()


par(mfrow=c(4,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
saveRDS(cells_assignment, file=paste0(AUCell_dir, "cells_assignment_",OUTPUT_NAME,".rds"))

Cell_type_assigned <- sapply(1:length(cells_assignment), function(i) cells_assignment[[i]][['assignment']], simplify = F)
names(Cell_type_assigned) <- names(cells_assignment)



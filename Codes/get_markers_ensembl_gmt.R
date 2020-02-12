source('Codes/Functions.R')
Initialize()

cardiac_marker_gene_set_df = read.table('Data/cardiac_markers.gmt', fill=T)
rownames(cardiac_marker_gene_set_df) <- cardiac_marker_gene_set_df[,1]

cardiac_marker_gene_set_df <- cardiac_marker_gene_set_df[,-c(1,2)]
cardiac_marker_gene_set <- lapply(1:nrow(cardiac_marker_gene_set_df), 
                                function(i) {
                                  cell_type_markers <- as.character(cardiac_marker_gene_set_df[i,])
                                  cell_type_markers[cell_type_markers!='']} )

names(cardiac_marker_gene_set) <- rownames(cardiac_marker_gene_set_df)


####  Convert gene list to ensembl id   
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)
cardiac_markers_ensembl_df <- lapply(cardiac_marker_gene_set, function( a_gene_set){getBM(filters="hgnc_symbol", 
                                                                              attributes= c('hgnc_symbol',"ensembl_gene_id"),
                                                                              values=a_gene_set, mart= ensembl)})
cardiac_markers_ensembl <- lapply(cardiac_markers_ensembel_df, function(x) x$ensembl_gene_id)
saveRDS(cardiac_markers_ensembl_df, 'Data/cardiac_markers_ensembl_mapped_table.rds')


sink("Data/cardiac_markers_ensembl.gmt")
for(i in 1:length(cardiac_markers_ensembl)){
  marker <- cardiac_markers_ensembl[[i]]
  cell_type_name <- names(cardiac_markers_ensembl)[i]
  
  cat(cell_type_name)
  cat('\t')
  cat(cell_type_name)
  cat('\t')
  for (j in 1:length(marker)){
    cat (marker[j])
    cat('\t')
  }
  cat('\n')
}
sink()



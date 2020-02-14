#!/bin/bash

cd ~/HeartDev/
for input_name in in_vivo in_vitro; do
  for input_file in objects/${input_name}/2.seur_dimRed_* ; do 
    echo ${input_name}
    echo ${input_file##*/}
    Rscript Codes/get_clusters.R ${input_name} ${input_file##*/}
    Rscript Codes/get_expression_tsne.R ${input_name} ${input_file##*/}
  done
done




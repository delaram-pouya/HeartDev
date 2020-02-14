#!/bin/bash

# INPUT_NAME = 'in_vivo'
# INPUT_FILE = '2.seur_dimRed_in_vivo_mito_0.593_lib_1500.rds'
cd ~/HeartDev/
for input_name in in_vivo in_vitro; do
  for input_file in objects/${input_name}/2.seur_dimRed_* ; do 
    echo ${input_name}
    echo ${input_file##*/}
    Rscript Codes/get_labels_AUCell.R ${input_name} ${input_file##*/}
    Rscript Codes/get_labels_SCINA.R ${input_name} ${input_file##*/}
    Rscript Codes/get_integrated_labels.R ${input_name} ${input_file##*/}
  done
done



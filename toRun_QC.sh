#!/bin/bash

# for the in-vitro data, lets use MAD= 2, 4, 6
# for in vivo 
# 1. remove all the MT genes and ignore them
# 2. use mit_cut_off = 0.5, 1, 2

Rscript Codes/Initialize.R 'in_vitro' 2 1500
Rscript Codes/Initialize.R 'in_vitro' 4 1500
Rscript Codes/Initialize.R 'in_vitro' 6 1500
Rscript Codes/Initialize.R 'in_vitro' 4 2000

Rscript Codes/Initialize.R 'in_vivo' 0.5 3000
Rscript Codes/Initialize.R 'in_vivo' 1 3000
Rscript Codes/Initialize.R 'in_vivo' 2 3000

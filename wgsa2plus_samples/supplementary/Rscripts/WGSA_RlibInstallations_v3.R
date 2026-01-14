#! /usr/bin/env Rscript


# for magicians
install.packages(c("devtools", "docopts","remotes")) #, "RColorBrewer","htmltools"))

#transform data
install.packages(c("easypackages","broom", "tidyverse", "data.table" ,"reshape2", "scales"))# "grid" is depreciated
#c("stringr", "readr", "dplyr", "tidyr", "tibble", "ggplot2") # part of tidyverse

#biological packages
install.packages(c("fossil","vegan")) ##, "survival", "ape", "limma"))
remotes::install_github("kasperskytte/ampvis2")
## info fo fossil/chao1: https://www.rdocumentation.org/packages/fossil/versions/0.4.0/topics/chao1

# tools for Bioc-magicians (old/new version packages coalitions, expect problems, its a mess)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") # there might be problems with newer versions and older are not supported
BiocManager::install("biomformat")
BiocManager::install("S4Vectors")
# BiocManager::install("dada2")
# BiocManager::install("DECIPHER")
# BiocManager::install("Biostrings")
# BiocManager::install("phyloseq")
# #source('http://bioconductor.org/biocLite.R')
# #biocLite('phyloseq')





#------------------- some  session info
#> R.Version()
#$platform
# [1] "x86_64-pc-linux-gnu"
#$version.string
#[1] "R version 4.1.0 (2021-05-18)"
#$nickname
#[1] "Camp Pontanezen"


# > library(BiocManager)
# Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.0 (2021-05-18)

# > BiocManager::valid()
# 
# * sessionInfo()
# 
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/local/intel/compilers_and_libraries_2020.2.254/linux/mkl/lib/intel64_lin/libmkl_rt.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BiocManager_1.30.16
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.1.0 tools_4.1.0   
# 
# Bioconductor version '3.14'

# * sessionInfo()
# 
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)




# installations for dada2 pipeline. 
# Need libraries:
# library(ggplot2)
# library(readr)
# library(reshape2)
# library(tidyr)
# library(dada2)
# library(biomformat)
# library(DECIPHER)
# library(S4Vectors)
# library(docopt)
# 
#all libraries are OK installed through R studio on local machine, and will be used on commandline R
#check .libPath()

# install.packages("tidyverse")
# install.packages("ggplot2") #part of tidyverse
# install.packages("readr") #part of tidyverse
# install.packages("tidyr") #part of tidyverse
install.packages("reshape2") #install through RStudio
install.packages("docopt") #install through RStudio



if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("S4Vectors") #install through RStudio

BiocManager::install("dada2") #install through RStudio
BiocManager::install("biomformat", force =T) #installs older version with dada2, but needs force-install through RStudio 
# (get the source latest of RCulr)
BiocManager::install("DECIPHER") #install through RStudio

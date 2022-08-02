#/usr/bin/env Rscript

# metadata needs to have columns:
# SampleID column
# Treatment (if div=T, > 2 factors; if div=F, 1 factor is ok)
# Metadata1 (if div=T, > 2 factors)

wdir=${1}
mappingFile=${2}
scriptCall=~/Documents/MyApps/RScripts/16Svis_pipeline/16S-Vis-Pip-v3.Rmd
reportName=16S_rdrTestReport # default is whatever RmD is titled
reportDIR=16S_rdrTest #default is ${1}
outputDir=~/Documents/MyApps/RScripts/16Svis_pipeline/16S_rdrTest 		#default is  ${1}/16Svis_plots 
title=YL27_transplant_to_OMM12_mice # cannot accept spaces or '' or ""
sampDEPTH=10000
asvDepth=10
diversity=T

## Creates a tempScript to run with Rscript
tempscript=`mktemp ~/Documents/MyApps/RScripts/16Svis_pipeline/tempscript.txt` || exit 1
echo "library(rmarkdown); 
rmarkdown::render('"${scriptCall}"', 'html_document', '"${reportName}"', '"${reportDIR}"')" > ${tempscript}

## run it as Rscript but without the -e flag
Rscript ${tempscript} \
 --wdir ${wdir}\
 --map ${mappingFile}\
 --outdir ${outputDir}\
 --title ${title}\
 --sampDepth ${sampDEPTH}\
 --asvDepth ${asvDepth}\
 --div ${diversity}

rm ${tempscript}


### Rscript -e ... does not work with docopt options
# Rscript -e "library(rmarkdown); rmarkdown::render('"${scriptCall}"', 'html_document')"\
#  --wdir ${wdir}\
#  --map ${mappingFile}\
#  --outdir ${reportDIR}\
#  --title ${title}\
#  --sampDepth ${sampDEPTH}\
#  --asvDepth ${asvDepth}
 
#2> ${logfile}
 # --diversity ${diversity}
# -------------- to contunue with
		# choose if to perform alpha, beta diversity calculations and stats
# logfile=16Svis_rlog.txt



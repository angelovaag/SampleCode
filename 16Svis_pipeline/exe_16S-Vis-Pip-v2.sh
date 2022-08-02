#/usr/bin/env Rscript

logfile=16S_vis_Rlog.txt
wdir=${1}
mappingFile=${2}
outdir=16Svis_plots
title="YL27 transplant to OMM12 mice"
sampDEPTH=10000
asvDepth=10
diversity=F # choose if to perform alpha, beta diversity calculations and stats

Rscript ~/Documents/MyApps/RScripts/16S_vis_pipeline/16S_Vis_Pip-v2.R\
 --wdir      ${wdir}\
 --map       ${mappingFile}\
 --outdir    ${outdir}\
 --title     ${title}\
 --sampDepth ${sampDEPTH}\
 --asvDepth  ${asvDepth} \
 --diversity ${diversity} \

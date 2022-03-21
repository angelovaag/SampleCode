#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --job-name mergTabs
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x.txt
#SBATCH --gres=lscratch:2 

#if [ ! -d "slogs" ]; then mkdir slogs; fi
#module load R
#supScriptPATH=/home/angelovaag/MyScripts/Nephele_WmGS_v2/v2.1_opt

supScriptPATH=~/Documents/Nephele_pipelines/WmGS_v2/tests/edits_4Nephele/step10_v3 #"< ************** NEPHELE/hard/path/to/R Scripts *********************** >"

echo $1
TAXorPWY=${1}                                 # Chose "PWY" or "TAX" script.
wdir=$(pwd)                                   #"< ************** NEPHELE/hard/path/to/projectDIR *********************** >"
indir=${1}profiles/merged_tables              #"< ************** NEPHELE/hard/path/to/${1}profiles/merged_tables (no / at end of path) *********************** >"
mfile=${wdir}/mapping.txt                     #"< ************** NEPHELE/hard/path/to/mapping #file without ^# ************************* >"
outdir=${1}profiles/DiversityPlots            #"< ************** NEPHELE/hard/path/to/${1}profiles/DiversityPlots (no / at end of path) *********************** >"

if [[ ! -d ${outdir} ]]; then mkdir -p ${outdir}; fi
sed 's/^#//g' ${mfile} > ${wdir}/mapping1.txt && mv ${wdir}/mapping1.txt ${wdir}/mapping.txt

echo $mfile
echo "---> Rscript start"
Rscript ${supScriptPATH}/${1}_DiversityPlots_v3.R \
  --wdir      ${wdir}  \
  --indir     ${indir} \
  --mfile     ${mfile} \
  --outdir    ${outdir}  


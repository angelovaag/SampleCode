#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --job-name mergTabs
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x.txt
#SBATCH --gres=lscratch:2 

if [ ! -d "slogs" ]; then mkdir slogs; fi
module load R
suppPATH=/home/angelovaag/MyScripts/Nephele_WmGS/v2.2_opt/Rscripts
# suppPATH=~/Documents/Nephele_pipelines/WmGS/tests

echo $1
TAXorPWY=${1}                           # Chose "PWY" or "TAX" script. For now restrict to TAX
wdir=$(pwd)                             # project working directory
mfile=mapping.txt                       # mapping file having no 1s/^#

if [ ${TAXorPWY} == "TAX" ]; then
  indir=TAXprofiles/ReadsTAX/merged_tables        # path to merged_tables (no / at end of path)
  outdir=TAXprofiles/ReadsTAX/DivPlots            # path for ouputdir: Diversity plots (no / at end of path)
elif [ ${TAXorPWY} == "PWY" ]; then
  if [ -d PWYprofiles/keggPWY/merged_tablesMP ]; then
    indir=PWYprofiles/keggPWY/merged_tablesMP        # path to merged_tables (no / at end of path)
    outdir=PWYprofiles/keggPWY/DivPlotsMP      # path for ouputdir: Diversity plots (no / at end of path)
  elif [ -d PWYprofiles/keggPWY/merged_tablesH3 ]; then
    indir=PWYprofiles/keggPWY/merged_tablesH3        # path to merged_tables (no / at end of path)
    outdir=PWYprofiles/keggPWY/DivPlotsH3     # path for ouputdir: Diversity plots (no / at end of path)
  elif [ -d PWYprofiles/metacycPWY/merged_tables ]; then
    indir=PWYprofiles/keggPWY/merged_tablesMP       # path to merged_tables (no / at end of path)
    outdir=PWYprofiles/keggPWY/DivPlotsMP      # path for ouputdir: Diversity plots (no / at end of path)
  fi
fi

if [[ ! -d ${outdir} ]]; then mkdir -p ${outdir}; fi
sed -i 's/^#//g' ${mfile} 



#cat $mfile
echo "---> Rscript start"
Rscript ${suppPATH}/DiversityPlots_${TAXorPWY}_v3.R \
  --wdir      ${wdir}  \
  --indir     ${indir} \
  --mfile     ${mfile} \
  --outdir    ${outdir} 

## removed only at completion of pipeline
# if [ ${TAXorPWY} == "TAX" ]; then unlink TAXprofiles/ReadsTAX; fi 
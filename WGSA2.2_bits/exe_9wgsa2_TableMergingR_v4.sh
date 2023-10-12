#!/bin/bash -e
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --job-name mergTabs
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x.txt
#SBATCH --gres=lscratch:2 
if [ ! -d "slogs" ]; then mkdir slogs; fi
module load R
module load biom-format
suppPATH=/home/angelovaag/MyScripts/Nephele_WmGS/v2.2_opt/Rscripts
# suppPATH=~/Documents/Nephele_pipelines/WmGS/tests

## *********** Settings
TAXorPWY=${1}                                # Chose "PWY" or "TAX" script
wdir=$(pwd)                                  # wdir will be where script is started from


if [[ "${TAXorPWY}" == "TAX" ]]; then 
  ## create link to TAXfolder for R script
  # inDIR=$(echo $(ls -d TAXprofiles/ReadsTAX_* ) ) 
  # linkNAME=TAXprofiles/ReadsTAX
  # ln --relative -s ${inDIR} ${linkNAME}
  ## should have been created in 3wgsa2

  binDIR=TAXprofiles/ReadsTAX/bin; echo "${binDIR}"
  genesDIR=NA; echo ${genesDIR} # genesDIR == NA for TAX mode
  prefix=wTAXid_4krona
  outdir=TAXprofiles/ReadsTAX/merged_tables

elif [[ "${TAXorPWY}" == "PWY" ]]; then
  if [ -d "PWYprofiles/keggPWY" ]; then 
    binDIR=PWYprofiles/keggPWY/pwybin.MP
    genesDIR=PWYprofiles/keggPWY/genebin
    prefix=ko2gg.MP
    outdir=${TAXorPWY}profiles/keggPWY/merged_tablesMP     # path & name an ouput dir
  else 
    binDIR=PWYprofiles/metacycPWY/pwybin.MP
    genesDIR=PWYprofiles/metacycPWY/genebin          # path/to/krona/files
    prefix=ec2mc.MP
    outdir=${TAXorPWY}profiles/metacycPWY/merged_tablesMP     # path & name an ouput dir
  fi
fi
# sleep 5s
# \*****

## *********** II.  Code on merging Tables for all samples
### 1)  Code for dealing with no-annotations samples that would cause an error
find ${binDIR}/*${prefix}.txt -type f -empty -delete   #delete empty krona files from failed samples
ls -1 ${binDIR} |xargs basename -a -s _${prefix}.txt  > SccList.txt #re-list non-empty existing krona files, in the correct order
list=SccList.txt

### 2) Code for R script that merges 
Rscript ${suppPATH}/TableMerging_${TAXorPWY}_v4.R \
  --binDIR ${binDIR} \
  --sccList ${list} \
  --outdir ${outdir} \
  --genesDIR ${genesDIR} && \
  echo "--- R script finished succesfully"

# if [[ "${TAXorPWY}" == "TAX" ]]; then unlink ${linkNAME} ; fi
# \***

## ********* III. Code on merging KO and EC gene tables. Using the SccList from step II.1
######## The capacity was added to the R script and there is no additional unix code

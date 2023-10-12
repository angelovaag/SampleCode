#!/bin/bash -e
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --job-name biomTabs
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x.txt
#SBATCH --gres=lscratch:2 
if [ ! -d "slogs" ]; then mkdir slogs; fi
module load biom-format
suppPATH=/home/angelovaag/MyScripts/Nephele_WmGS/v2.2_opt

## *********** Settings
TAXorPWY=${1}                         # Chose "PWY" or "TAX" script
wdir=$(pwd)                           # wdir will be where script is started from
mfile=mapping.txt

if [[ "${TAXorPWY}" == "TAX" ]]; then 
  binDIR=TAXprofiles/ReadsTAX/bin
  biomSuffix=Lineage
  tabtype='OTU table'                               # quotes included
  biomsDIR=TAXprofiles/ReadsTAX/bioms
  genesDIR=NA #"genesDIR == NA for TAX mode
  prefix=wTAXid_4krona
  outdir=${1}profiles/ReadsTAX/merged_tables     # path & name an ouput dir

elif [[ "${TAXorPWY}" == "PWY" ]]; then
  biomSuffix=allTiers
  tabtype='Pathway table'                           # quotes included
  if [ -d "PWYprofiles/keggPWY/pwybin.MP" ]; then 
    binDIR=PWYprofiles/keggPWY/pwybin.MP
    biomsDIR=PWYprofiles/keggPWY/bioms.MP
    genesDIR=PWYprofiles/keggPWY/genebin
    prefix=ko2gg.MP
    outdir=${1}profiles/keggPWY/merged_tablesMP     # path & name an ouput dir
  
  elif [ -d "PWYprofiles/keggPWY/pwybin.H3" ]; then 
    binDIR=PWYprofiles/keggPWY/pwybin.H3
    biomsDIR=PWYprofiles/keggPWY/bioms.H3
    genesDIR=PWYprofiles/keggPWY/genebin
    prefix=ko2gg.H3
    outdir=${1}profiles/keggPWY/merged_tablesH3 
  
  elif [ -d "PWYprofiles/metacycPWY/pwybin.MP" ]; then 
    binDIR=PWYprofiles/metacycPWY/pwybin.MP
    biomsDIR=PWYprofiles/metacycPWY/bioms.MP
    genesDIR=PWYprofiles/metacycPWY/genebin          # path/to/krona/files
    prefix=ec2mc.MP
    outdir=${1}profiles/metacycPWY/merged_tablesMP     # path & name an ouput dir
  fi
fi
# sleep 5s
# \*****

# ********** I. Code for creating bioms for each sample
## Creating biom files per sample
echo "---- Creating biom files per sample"
if [ ! -d "${biomsDIR}" ]; then mkdir -p ${biomsDIR}; fi

### Give the Biom table a "Taxonomy" column, so MicrobiomeMB can recognize its 'taxonomy'
ls -1 ${binDIR} |xargs basename -a -s "_${prefix}.txt" > SccList.txt

cat SccList.txt |xargs -I {} sh -c \
 "sed 's/\t/;/g' ${binDIR}/{}*.txt | sed 's/;/\t/1' |nl -n ln |sed '1s/^/\t/' | \
  sed '1s/1\s\+\t#/\t/' | sed 's/\s\+\t/\t/g' | \
  sed -e '1s/kingdom.*/Taxonomy/' |sed -e '1s/Level.*/Taxonomy/' > ${biomsDIR}/{}_4biom.txt"

# echo "check _4biom.txts now"
# sleep 10s

### Creating the individual sample bioms (pre table merging)
cat SccList.txt |xargs -I {} sh -c \
 "biom convert -i ${biomsDIR}/{}_4biom.txt -o ${biomsDIR}/{}_json.biom \
        --to-json --table-type=\"${tabtype}\" --process-obs-metadata taxonomy"
rm ${biomsDIR}/*_4biom.txt
# \*****


# ******* II. Code for creating biom from merged table (post table merging)
echo "---- Creating merged_tables biom file"
sed '1s/Lineage/Taxonomy/' ${outdir}/merged_Counts+${biomSuffix}.txt | sed '1s/allTiers/Taxonomy/' > ${outdir}/tmp_4biom.txt
biom convert -i ${outdir}/tmp_4biom.txt -m <(sed '1s/^/#/' ${mfile}) \
        -o ${outdir}/merged_Counts+${biomSuffix}_json.biom \
        --to-json --table-type="${tabtype}" --process-obs-metadata taxonomy &&\
rm ${outdir}/tmp_4biom.txt

# \***

#!/bin/bash -e
#$ -N minion                             # set job name
#$ -cwd                                  # work from current directory
#$ -j y                                  # use 1 log file
#$ -o minion_profiling_log.txt           # set log file name
#$ -l h_vmem=14G                         # set memory requirements for script (per core) [ default: h_mem=8G ]
#$ -m be                                 # email notification preferences
##$ -M <email>                           # provide email , remove 2' #
#$ -pe threaded 16                       # thread / core requirements for script

NSLOTS=16
module load kraken
# module load python
module load krona

list=${1%%.*}
DBname=${2}
suffix=${3}     # _asmb or _raw or _something-to-clarify-folder-name (optional)
indir=raw
outdir=profiles_${DBname}${suffix}

if [[ $DBname == "plus" ]]; then
export KRAKEN_DB=/hpcdata/bio_data/kraken_db/standardPLUSpf+nema/
elif [[ $DBname == "plus_wTmu" ]]; then
export KRAKEN_DB=/hpcdata/bio_data/kraken_db/standardPLUSpvf+TmuMGX
elif [[ $DBname == "GTdb" ]]; then
export KRAKEN_DB=/hpcdata/bcbb/poorani/dbs/kraken2/mydb/gtdb_r95_plusEu
elif [[ $DBname == "MGBCdb" ]]; then
export KRAKEN_DB=/hpcdata/bio_data/kraken_db/MGBCdb/MGBC-26640_KrakenBracken
elif [[ $DBname == "TrichDB" ]]; then
export KRAKEN_DB=/hpcdata/bio_data/kraken_db/TrichDB/
elif [[ $DBname == "decontamDB" ]]; then
export KRAKEN_DB=/hpcdata/bio_data/kraken_db/decontam_db/decontam_human+mouse_db/
fi

ext=F # or T

     if [ ! -d "${outdir}" ]; then mkdir ${outdir}; fi

     params=()
     if [[ ${ext} == "T" ]]; then
          params+=(--gzip-compressed)
     fi

for i in $(cat ${list}.txt); do
     echo "--------------- TAX classification for file ${i} against database $DBname"
       kraken2 --use-names "${params[@]}" --threads $NSLOTS \
               --confidence 0.05 --db $KRAKEN_DB \
               --output ${outdir}/${i}_kr2raw.txt \
               --report ${outdir}/${i}_kr2REPORT_$DBname.txt \
               --classified-out ${outdir}/${i}_classified_$DBname.fasta \
               --unclassified-out ${outdir}/${i}_unclassified_$DBname.fasta \
               ${indir}/${i}*.fasta

     grep "C" ${outdir}/${i}_kr2raw.txt |cut -f1-3 > ${outdir}/classified_headers.txt
                   
     python /hpcdata/bio_data/kraken_db/kreport2krona.py -r ${outdir}/${i}_kr2REPORT*.txt \
               -o ${outdir}/${i}_4krona.txt
done
ktImportText ${outdir}/*_4krona.txt -o ${outdir}/sample_profiles_$DBname.html && \
          rm ${outdir}/*_4krona.txt # ${outdir}/*_kr2raw.txt



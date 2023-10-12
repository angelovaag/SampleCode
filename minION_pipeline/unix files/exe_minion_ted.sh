#!/bin/bash -e


module load fastp
module load kraken
module load krona
# module load seqtk

NSLOTS=16
list=${1%%.*}
indir=raw
outdir=${indir}/ted
outdir2=${outdir}/contaminated
ext=${2} #use .gz or leave empty

export KRAKEN_DB=/hpcdata/bio_data/kraken_db/decontam_db/decontam_human+mouse_db/
e=10                                                              # threshold for whole-read average quality
l=250                                                             # threshold for min read length
lq=20                                                             # threshold for 5' end quality trimming  
rq=15                                                             # threshold of 3' end quality trimming
W=4                                                               # trimming window size (default: 4)

if [[ ! -d ${outdir2} ]]; then mkdir -p ${outdir2}; fi

params=()
if [[ ${ext} == ".gz" ]]; then params+=(--gzip-compressed) && ext=.gz; fi

for f in $(cat ${list}.txt); do
     echo "---------------
          trimming & filtering with fastp for" ${f}

      fastp -i ${indir}/${f}*.fastq${ext} \
                      -h ${outdir}/${logs}/${f}_log.html -j ${outdir}/${logs}/${f}_log.json \
                -o ${outdir}/${f}_TE.fastq${ext} -y -c  -A --trim_poly_x -e ${e} -l ${l} \
                -5 --cut_front_mean_quality ${lq} -3 --cut_tail_mean_quality ${rq}\
                -W ${W} -w 16 
     echo "---------------
          decontamination with Kraken2_db for" ${f}
      kraken2 --use-names "${params[@]}" --threads $NSLOTS \
               --confidence 0.05 --db $KRAKEN_DB \
               --output ${outdir2}/${f}_kr2raw_contamination.txt \
               --report ${outdir2}/${f}_kr2REPORT_contamination.txt \
               --classified-out  ${outdir2}/${f}_contamination.fasta \
               --unclassified-out ${outdir}/${f}_TED.fasta \
               ${outdir}/${f}_TE.fastq${ext}

     grep "C" ${outdir2}/${f}_kr2raw_contamination.txt |cut -f1-3 > ${outdir2}/contaminated_headers.txt
                   
     python /hpcdata/bio_data/kraken_db/kreport2krona.py -r ${outdir2}/${f}_kr2REPORT*.txt \
               -o ${outdir}/${f}_4krona.txt
done
ktImportText ${outdir}/*_4krona.txt -o ${outdir}/profile_TED.html && \
          rm ${outdir}/*_4krona.txt ${outdir2}/*_kr2raw_*  #${outdir}/${f}_TE.fastq${ext}



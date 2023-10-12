#!/bin/bash -e
#$ -N fastp	                # set job name
#$ -cwd                     # work from current directory
#$ -j y                     # use 1 log file
#$ -o fastp_trimming.txt    # set log file name
#$ -l h_vmem=8G             # set memory requirements for script (per core) [ default: h_mem=8G ]
#$ -m be                    # email notification preferences
#$ -M angelovaag@nih.gov    # provide email 
#$ -pe threaded 16          # thread / core requirements for script

module load fastp

list=${1%%.*}
indir=raw
outdir=${indir}/tfe
#trimming parameters for minION
e=10                                                              # threshold for whole-read average quality
l=250                                                             # threshold for min read length
lq=20                                                            # threshold for 5' end quality trimming  
rq=15                                                             # threshold of 3' end quality trimming
W=4                                                               # trimming window size (default: 4)

if [[ ! -d ${outdir} ]]; then mkdir -p ${outdir}; fi

for f in $(cat ${list}.txt); do
      echo "---------------
          trimming & filtering with fastp for" ${f}

          fastp -i ${indir}/${f}*.fastq.gz \
          		  -h ${outdir}/${logs}/${f}_log.html -j ${outdir}/${logs}/${f}_log.json \
                -o ${outdir}/${f}_tfe.fastq.gz -y -c  -A --trim_poly_x -e ${e} -l ${l} \
                -5 --cut_front_mean_quality ${lq} -3 --cut_tail_mean_quality ${rq}\
                -W ${W} -w 16 
done
             

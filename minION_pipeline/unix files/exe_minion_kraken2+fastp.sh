#!/bin/bash -e
#$ -N fastp	                # set job name
#$ -cwd                     # work from current directory
#$ -j y                     # use 1 log file
#$ -o fastp_trimming.txt    # set log file name
#$ -l h_vmem=14G             # set memory requirements for script (per core) [ default: h_mem=8G ]
#$ -m be                    # email notification preferences
#$ -M angelovaag@nih.gov    # provide email 
#$ -pe threaded 16          # thread / core requirements for script

NSLOTS=32
module load fastp
module load python
module load krona
module load pigz

indir=raw
outdir1=TFQreads
outdir2=kraken2
#trimming parameters for minION
l=250                                                             # threshold for min read length
W=4                                                               # trimming window size (default: 4)
export KRAKEN_DB=/hpcdata/bio_data/kraken_db/standardPLUSpf+nema/


 if [ ! -d "${outdir1}" ]; then 
      mkdir ${outdir1}
      mkdir ${outdir2}
  fi


for f in $(ls ${indir}/*fa*.gz ${indir}/*.fa ${indir}/*.fna |cut -f1 -d "." ); do
#     	echo "---------------   trimming adaptors with PoreChop for" ${f}
     	i=$(basename ${f})
#	python ~/MyTools/PoreChop_minion/Porechop/porechop-runner.py -i ${f}.fastq.gz \
		# -o ${outdir1}/${i}_tA.fastq.gz -t $NSLOTS
 		
  #       echo "---------------  trimming quality & filtering with fastp for" ${f}
  #    fastp   -i ${outdir1}/${i}_tA.fastq.gz -W ${W} -w 16 -y -c  -A \
  #       	-h ${outdir1}/${logs}/${i}_log.html -j ${outdir1}/${logs}/${i}_log.json \
  #         -o ${outdir1}/${i}_tfq.fastq.gz --trim_poly_x -l ${l} && \
  #         rm ${outdir1}/${i}_tA.fastq.gz
 		
           echo "--------------- classification with kraken for " ${i} 
     kraken2 --use-names --gzip-compressed --threads $NSLOTS \
          --confidence 0.0 --db $KRAKEN_DB \
          --output ${outdir2}/${i}_kr2raw.txt \
          --report ${outdir2}/${i}_kr2REPORT.txt \
          ${f}.fastq.gz
          #${outdir1}/${i}_tfq.fastq.gz

     python ~/MyScripts/krakenScripts/kreport2krona.py -r ${outdir2}/${i}_kr2REPORT.txt \
          -o ${outdir2}/${i}_4krona.txt
done       
     ktImportText ${outdir2}/*_4krona.txt -o ${outdir2}/sample_profiles.html



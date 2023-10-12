#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=2g
#SBATCH --job-name rawQC
#SBATCH --mail-type END
#SBATCH -o slurm_rawQC-log.txt
#SBATCH --time=8:00:00

module load fastqc
module load multiqc

#run this script in the folder with all the fastq.gz files

indir=raw
outdir=rawQC
fqc=fastQC


#create general output dir
if [ ! -d "${indir}/${outdir}" ]; then
mkdir ${indir}/${outdir}/ 
mkdir ${indir}/${outdir}/${fqc}/
fi


#run fastqc
for i in $(ls ${indir}/*fastq.gz |cut -f1 -d "." ); do
	b=$(basename ${i%%})
	#ls ${indir}/${outdir}/${fqc}/${b}*html
	
	 if [ -s ${indir}/${outdir}/${fqc}/${b}*html  ]; then
	  	echo "------- skipping sample ${b}*html, for being done before"
	 else
	 	echo "------- performing QC for ${b}"
	 	fastqc ${i}.fastq.gz -o ${indir}/${outdir}/${fqc} -t $SLURM_CPUS_PER_TASK
	 fi

 done

multiqc ${indir}/${outdir}/${fqc} -o ${indir}/${outdir} --title rawQC

rm ${indir}/${outdir}/${fqc}/*_fastqc.zip



#will not be using prinseq-lite as it uses uncompressed .fastq files
#poorani says she has a batch script for this:
#https://github.niaid.nih.gov/subramanianp4/runQC
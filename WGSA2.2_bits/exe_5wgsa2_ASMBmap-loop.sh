#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=96g
#SBATCH --job-name mapping
#SBATCH --output slogs/slog_%x_%j.txt
#SBATCH --mail-type END
#SBATCH --time=5-00:00:00
#SBATCH --gres=lscratch:800
#SBATCH --partition=norm  #### quick/norm, quick walltime must be < 4hr 
#### takes ~1hr per sample. Cannot get it to efficiently use its memory
#lscratch_mon --log slogs/lscratchUSAGE_bt2.log --interval 1 & #super important to add "&" at end!!!!
if [ ! -d "slogs" ]; then mkdir slogs; fi

module load bowtie/2
module load samtools
module load metabat
module load bbtools

list=${1%%.*} 								# list of indiv_vIDs/Sample names
asmbDIR=asmbMetaSpades		 				#asmbMetaSpades 	
tedPATH=TEDreads 		
infile=final.assembly  	
tmpDIR=/lscratch/$SLURM_JOB_ID 				# for large tmp files

for i in $(cat ${list}.txt); do
f=${asmbDIR}/${i}_asmb
 	printf "\n---------> Starting work on sample ${i}, in path ${f}"
 	# mapping using bt2 
	 	if [ -s ${f}/${infile}.bam.bai ]  && [ -s ${f}/${infile}_scaffCoverage.txt ]; then
	 		echo "---------> mapping has been done to assembly. Skipping ${i}"
	 	else
	 		if [ ! -s ${tmpDIR}/${i}.db.rev.1.bt2 ]; then
	 		echo "	-----> bt2 DB preparation for" ${d}
		 	bowtie2-build ${f}/${infile}.fasta ${tmpDIR}/${i}.db --threads $((SLURM_CPUS_PER_TASK - 2 )) #give bowtie 2 "breathing room" threads cuz it sometimes uses it
		 	else
		 	echo " -----> bt2 DB for ${i} found. Proceeding to mapping"
		 	fi

		echo "---------> mapping of ${i} to bt2 DB of ${d}"
		printf "\nTEDreads mapping stats:\n"  >> ${f}/${infile}_stats.txt				 # cuz some info will be added from next step
		bowtie2 --phred33 --sensitive-local --no-unal --seed 4 \
				-1 ${tedPATH}/${i}_R1_ted.fastq.gz   -2 ${tedPATH}/${i}_R2_ted.fastq.gz \
				-x ${tmpDIR}/${i}.db 			  	 -S ${tmpDIR}/${i}.sam \
				-p $((SLURM_CPUS_PER_TASK - 2 ))     2>> ${f}/${infile}_stats.txt		# here is the info added 

	#sam->bam, de-replicating alignments, indexing final bam file, removing tmp files
		echo "---------> convert, sort, collate, fix-mate, dedup & index alignment file for " ${i}
		samtools sort  	 		${tmpDIR}/${i}.sam 			-o ${tmpDIR}/${i}.bam 		  -@ $SLURM_CPUS_PER_TASK	&& rm ${tmpDIR}/${i}.sam
		samtools collate 		${tmpDIR}/${i}.bam 			-o ${tmpDIR}/${i}_colated.bam -@ $SLURM_CPUS_PER_TASK	&& rm ${tmpDIR}/${i}.bam
		samtools fixmate     -m	${tmpDIR}/${i}_colated.bam 	   ${tmpDIR}/${i}_fixmate.bam -@ $SLURM_CPUS_PER_TASK 	&& rm ${tmpDIR}/${i}_colated.bam
		samtools sort 			${tmpDIR}/${i}_fixmate.bam  -o ${tmpDIR}/${i}_fixmate_sorted.bam 					&& rm ${tmpDIR}/${i}_fixmate.bam
		samtools markdup  -r -s ${tmpDIR}/${i}_fixmate_sorted.bam \
							 -f ${tmpDIR}/${i}_tmp_markdup.txt         ${f}/${infile}.bam -@ $SLURM_CPUS_PER_TASK 	&& rm ${tmpDIR}/${i}_fixmate_sorted.bam
		samtools index 		 -b ${f}/${infile}.bam 										  -@ $SLURM_CPUS_PER_TASK 
	#adding de-replication stats to final.assembly stats
		cat ${f}/${infile}_stats.txt <(printf "\nRead alignment de-replication stats:\n") \
			${tmpDIR}/${i}_tmp_markdup.txt > ${tmpDIR}/${i}_stats.txt &&
		mv ${tmpDIR}/${i}_stats.txt ${f}/${infile}_stats.txt

		
	#getting scaffold coverage: it is not for TAX abund estimations.
		echo "---------> getting scaffoldCov for ${i} with bbtools::pilieup"
		bbtools pileup ${f}/${infile}.bam    out=${f}/tmp_scaffcov0.txt overwrite=true
		sort ${f}/tmp_scaffcov0.txt   >    ${f}/tmp_scaffcov1.txt 														##**new output name**
	
	#getting scaffold depth (per base coverage). This one is for TAX abund extimations. Used for xbin, metabat autometa & more
		echo "---------> creating the baseCov for ${f} with metabat2::jgi_depths script"
		jgi_summarize_bam_contig_depths ${f}/${infile}.bam --outputDepth ${f}/${infile}_depths.txt 
		sort ${f}/${infile}_depths.txt | sed '1 s/^/#/' |cut -f 1-3 > ${f}/tmp_basecov.txt
	
	#getting read counts
		echo "---------> getting readCounts or ${i} with samtools::idxstats"
		samtools idxstats ${f}/${infile}.bam  > ${f}/tmp_idx.txt
		cut -f 1-3 ${f}/tmp_idx.txt | sed -e '1s/^/#contigName\tLength\tReadsCount\n/' | \
				sort | grep -v "*" > ${f}/tmp_readcount.txt

	##one file to rule them all ( no #TPM conversion, since no sense for for scaffolds (they are not genes)):
		echo "---------> joining read counts & scaff ave coverage & scaff depths into one file to rule them all >> scaffCov"
		paste ${f}/tmp_readcount.txt 	${f}/tmp_scaffcov1.txt  ${f}/tmp_basecov.txt | \
		cut -f 1,2,3,8-14,17 -d $'\t' > ${f}/tmp_scaffcov2.txt 															##**new output name**
		
		## calculating RPK per scaffold 																				##**new code **
		echo "---------> calculating RPK >> scaffCov"
		awk 'BEGIN{FS=OFS="\t"} NR>1 {if ($2>0 && $3>0) $12=sprintf("%0.2f", $3*1e3/$2); else $12=0; print}' ${f}/tmp_scaffcov2.txt | \
	 	sed -e '1s/^/#NODE\tlen\tReads\t%ScaffCovered\t%bpCovered\tplusReads\tminusReads\t%GC\tMedFold\tstDev\tAveDepth\tRPK\n/' > ${f}/${infile}_scaffCoverage.txt
		rm ${f}/tmp_*.txt

		fi


		echo "---------> exe_asmbBowtie2.sh of ${l} done "

done
cat ${list}.txt |xargs -I {} bash -c  \
"if [ ! -s ${asmbDIR}/{}_asmb/${infile}.bam ]; then echo 'mapping failed for {}' >> failed_samples.txt; fi "
echo "exe_asmbBowtie2.sh of ${1} done"




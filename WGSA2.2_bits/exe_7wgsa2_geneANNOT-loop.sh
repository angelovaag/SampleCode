#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=96g
#SBATCH --job-name eggnog
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x_%j.txt
#SBATCH --gres=lscratch:800
#SBATCH --time=3-00:00:00
if [ ! -d "slogs" ]; then mkdir slogs; fi

export TMPDIR=/lscratch/$SLURM_JOB_ID
lscratch_mon --log slogs/lscratchUSAGE.log --interval 1 & #super important to add "&" at end!!!!
module load eggnog-mapper
module load seqtk

list=${1%%.*} 												#list of indiv_vIDs/SampleName
 infile=PREDgenes
 outdir=annotations
 prefix=annots
 ABUNfile=PREDgenes_ABUNtab.txt
dirPATH=asmbMetaSpades/{}_asmb/genes/${outdir}
suppPATH=/home/angelovaag/MyScripts/Nephele_WmGS/v2.2_opt
suppPATH2=${suppPATH}/supplementary
str=NODE_ 													# or TRINITY_ for trinity assemblies (not implemented yet)

### copy database to /dev/shm (biowulf convenience only)
	shmDIR=/dev/shm/eggnogDATA
	mkdir -p ${shmDIR}
    echo 'copying eggnog database to /dev/shm'
    free -h --wide
    cp $EGGNOG_DATA_DIR/eggnog.db ${shmDIR}/
    cp $EGGNOG_DATA_DIR/eggnog_proteins.dmnd ${shmDIR}/
    free -h --wide

for sampName in $(cat ${list}.txt ); do
	path=asmbMetaSpades/${sampName}_asmb/genes
	if [ ! -d ${path}/${outdir}/ ]; then mkdir -p ${path}/${outdir}/; fi

	if [[ -s ${path}/${outdir}/${prefix}.KEGGmap.txt ]]; then
		echo "-------> skipping annotation for ${sampName}, because done"
	else
		echo "-------> functional annotation for ${sampName}"
			
		emapper.py -i ${path}/${infile}.faa -o ${prefix} --output_dir ${path}/${outdir}/ --cpu $SLURM_CPUS_PER_TASK \
					--itype proteins --block_size 4 --index_chunks 2 --data_dir ${shmDIR} \
					--override --temp_dir /lscratch/$SLURM_JOB_ID  --scratch_dir /lscratch/$SLURM_JOB_ID #{--dbmem, --resume,--override}
			
		#grep ec number annotations
		grep "${str}" ${path}/${outdir}/${prefix}.emapper.annotations | cut -f 1,11 | sed -e 's/\t-//g' | \
				grep -e $'\t' | awk '{n=split($2,s,",");for (i=1;i<=n;i++) {$2=s[i];print}}' |sed -e 's/\s/\t/g' \
				> ${path}/${outdir}/${prefix}.ec.txt
		#grep cog annotations
		grep "${str}" ${path}/${outdir}/${prefix}.emapper.annotations | cut -f 1,5 \
				> ${path}/${outdir}/${prefix}.COG.txt

		#grep KO numbers
		grep "${str}" ${path}/${outdir}/${prefix}.emapper.annotations |cut -f1,12 |grep "ko:" | \
				 sed 's/ko://g' | awk '{n=split($2,s,",");for (i=1;i<=n;i++) {$2=s[i];print}}' |sed -e 's/\s/\t/g' \
				 > ${path}/${outdir}/${prefix}.ko.txt

		#grep KEGG maps
		grep "${str}" ${path}/${outdir}/${prefix}.emapper.annotations |cut -f1,13 |grep "ko" | sed 's/,map.*//g' | \
				sed 's/ko/map/g' | awk '{n=split($2,s,",");for (i=1;i<=n;i++) {$2=s[i];print}}' |sed -e 's/\s/\t/g'  \
				 > ${path}/${outdir}/${prefix}.KEGGmap.txt
		mv ${path}/${outdir}/${prefix}.emapper.annotations ${path}/${outdir}/${prefix}.emapper.annotations.txt
	fi 

		## making annotGenes.faa and .fna
	 seqtk subseq ${path}/${infile}.faa ${path}/${outdir}/${prefix}.ko.txt > ${path}/${outdir}/${prefix}.ko.genes.faa
	 seqtk subseq ${path}/${infile}.fna ${path}/${outdir}/${prefix}.ko.txt > ${path}/${outdir}/${prefix}.ko.genes.fna
	 seqtk subseq ${path}/${infile}.faa ${path}/${outdir}/${prefix}.ec.txt > ${path}/${outdir}/${prefix}.ec.genes.faa
	 seqtk subseq ${path}/${infile}.fna ${path}/${outdir}/${prefix}.ec.txt > ${path}/${outdir}/${prefix}.ec.genes.fna

	 ## calculating KO-based abundances for each instance of gene (iTPM) & making instance-based baundance table
	 echo "---------> calclulating KO-based iTPM for ${sampName}"
	 awk -v OFS="\t" 'NR==FNR { id[$1]=$0; next } ($1 in id){ print $2, id[$1]}'  \
	 ${path}/${ABUNfile} ${path}/${outdir}/${prefix}.ko.txt |\
	 sed '1s/^/#KO\tnodeID\tlen\treads\tcov\tiRPK\tiTPM\n/' > ${path}/${outdir}/ANNOTgenes_ABUNtab.ko.txt

	 		## calculating geneTPM (TPM for each unique gene across all its instances) in a sample 
 	 echo "---------> calculating KO-based geneTPM for ${sampName}"
	 cut -f1,6 ${path}/${outdir}/ANNOTgenes_ABUNtab.ko.txt | sed $'s/\#KO.*//g' |\
	 awk -F'\t' -f ${suppPATH2}/calc_geneTPM.awk > ${path}/${outdir}/tmp_geneTPMtab.ko.txt 
	 awk -F'\t' -f ${suppPATH2}/join2files_f2.awk ${suppPATH2}/KO_list.txt ${path}/${outdir}/tmp_geneTPMtab.ko.txt | \
	 	sed '1s/^/#KO\tgeneTPM\tgeneNAME\n/'> ${path}/${outdir}/geneTPMtab.ko.txt #cp in PWYdir/genebin, next script
	 # rm ${path}/${outdir}/tmp_geneTPMtab.ko.txt


	## calculating EC-based abundances for each instance of gene (iTPM) & making instance-based baundance table
	 echo "---------> calclulating EC-based iTPM for ${sampName}"
	 awk -v OFS="\t" 'NR==FNR { id[$1]=$0; next } ($1 in id){ print $2, id[$1]}'  \
	 ${path}/${ABUNfile} ${path}/${outdir}/${prefix}.ec.txt |\
	 sed '1s/^/#EC\tnodeID\tlen\treads\tcov\tiRPK\tiTPM\n/' > ${path}/${outdir}/ANNOTgenes_ABUNtab.ec.txt
 	 
	 		## calculating geneTPM (TPM for each unique gene across all its instances) in a sample 
 	echo "---------> calculating EC-based geneTPM for ${sampName}"
	cut -f1,6 ${path}/${outdir}/ANNOTgenes_ABUNtab.ec.txt | sed $'s/\#EC.*//g' |\
 	awk -F'\t' -f ${suppPATH2}/calc_geneTPM.awk > ${path}/${outdir}/tmp_geneTPMtab.ec.txt 
	awk -F'\t' -f ${suppPATH2}/join2files_f2.awk ${suppPATH2}/EC_list.txt ${path}/${outdir}/tmp_geneTPMtab.ec.txt | \
	sed '1s/^/#EC\tgeneTPM\tgeneNAME\n/'> ${path}/${outdir}/geneTPMtab.ec.txt  #cp in PWYdir/genebin, next script
	
	rm ${path}/${outdir}/tmp_geneTPMtab.*.txt
done


cat ${list}.txt |xargs -I {} bash -c \
" if [ ! -s ${dirPATH}/${prefix}.emapper.annotations ]; then echo 'gene annotation failed for' {} >> failed_samples.txt; fi "

echo "exe_asmbEggNog.sh done"

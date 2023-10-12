#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2g
#SBATCH --job-name MAGs
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=norm 
if [ ! -d "slogs" ]; then mkdir slogs; fi

module load metabat
module load m-tools
module load samtools
module load kronatools

list=${1%%.*}
magDIR=MAGs
magQA=magsQA
fasfix=final.assembly
magfix=mag
ftype=fa
plotsDIR=QCplots
inPATH=asmbMetaSpades/{}_asmb
plots=F 
krona=T

set -e

jobs=4
threads=$(( SLURM_CPUS_PER_TASK / ${jobs}))



cat ${list}.txt |xargs -I {} sh -c \
" if [ ! -f ${inPATH}/${fasfix}_depths.txt ]; then echo '---------> Where is my depths file? Remaking ...' && \
	jgi_summarize_bam_contig_depths ${inPATH}/${fasfix}.bam  --outputDepth ${inPATH}/${fasfix}_depths.txt; fi"

echo "---------> sorting scaffolds into bins (MetaBat2)"
parallel -a ${list}.txt -j $(( SLURM_CPUS_PER_TASK / ${threads} )) \
	metabat2 -i ${inPATH}/${fasfix}.fasta -o ${inPATH}/${magDIR}/mag -m 1500 --maxEdges 250 \
	--unbinned  -t ${threads} -a ${inPATH}/${fasfix}_depths.txt 

echo "---------> mags QA:: lineage analysis (CheckM) "
parallel -a ${list}.txt -j $(( SLURM_CPUS_PER_TASK / ${threads} )) \
	checkm lineage_wf --tab_table --nt --pplacer_threads ${threads}  \
	-t ${threads} -x ${ftype} ${inPATH}/${magDIR}/ ${inPATH}/${magDIR}/${magQA}/ 
	
echo " ---------> mags QA:: completeness & contamination levels (CheckM)"
parallel -a ${list}.txt -j $(( SLURM_CPUS_PER_TASK / ${threads} )) \
	checkm qa -o 2 --tab_table -t ${threads} -f ${inPATH}/${magDIR}/${magQA}/${magfix}_qa.txt \
	${inPATH}/${magDIR}/${magQA}/lineage.ms ${inPATH}/${magDIR}/${magQA}/

echo " ---------> mags QA:: TAX (CheckM)"
parallel -a ${list}.txt -j $(( SLURM_CPUS_PER_TASK / ${threads} )) \
	checkm tree_qa -o 2 --tab_table -f ${inPATH}/${magDIR}/${magQA}/${magfix}_tax.txt \
	${inPATH}/${magDIR}/${magQA} 

echo " ---------> mags QA:: coverage (CheckM)"
parallel -a ${list}.txt -j $(( SLURM_CPUS_PER_TASK / ${threads} )) \
	checkm coverage -x ${ftype} -t ${threads}   ${inPATH}/${magDIR}   \
	${inPATH}/${magDIR}/${magQA}/${magfix}_coverage.txt  ${inPATH}/${fasfix}.bam

echo "----------> bin QA:: profile assessments (CheckM)"
parallel -a ${list}.txt -j $(( SLURM_CPUS_PER_TASK / ${threads} )) \
	checkm profile -q ${inPATH}/${magDIR}/${magQA}/${magfix}_coverage.txt \
	--tab_table -f ${inPATH}/${magDIR}/${magQA}/${magfix}_profiles.txt

echo "----------> collecting data into 'one file to rule them all' (summary file)"
cat ${list}.txt |xargs -I {} sh -c \
"cut ${inPATH}/${magDIR}/${magQA}/${magfix}_tax.txt -f1,5 > ${inPATH}/${magDIR}/${magQA}/tmp1.txt && \
 cut ${inPATH}/${magDIR}/${magQA}/${magfix}_qa.txt -f1,6-11,13,15,17,19,23 > ${inPATH}/${magDIR}/${magQA}/tmp2.txt"
cat ${list}.txt |xargs -I {} sh -c \
"paste ${inPATH}/${magDIR}/${magQA}/tmp1.txt      ${inPATH}/${magDIR}/${magQA}/tmp2.txt \
	> ${inPATH}/${magDIR}/${magQA}/${magfix}_tmpsum.txt && rm ${inPATH}/${magDIR}/${magQA}/tmp* "
cat ${list}.txt | xargs -I {} sh -c \
"join ${inPATH}/${magDIR}/${magQA}/${magfix}_profiles.txt ${inPATH}/${magDIR}/${magQA}/${magfix}_tmpsum.txt -t $'\t' | \
	cut -f 1,3,4,6,7-19 |sed 's/final.assembly: //g' > ${inPATH}/${magDIR}/${magQA}/${magfix}_SUMMARY.txt && \
 rm ${inPATH}/${magDIR}/${magQA}/${magfix}_tmpsum.txt"

if [ ${plots} == T ]; then
echo "----------> creating bin summary plots (CheckM)"
parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	checkm coding_plot --image_type pdf  -x ${ftype}   \
	${inPATH}/${magDIR}/${magQA}/ ${inPATH}/${magDIR}/  \
	${inPATH}/${magDIR}/${magQA}/${plotsDIR}/    95

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	checkm marker_plot --image_type pdf  -x ${ftype} --dpi 400 \
	${inPATH}/${magDIR}/${magQA}/ ${inPATH}/${magDIR}/  \
	${inPATH}/${magDIR}/${magQA}/${plotsDIR}/ 

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	checkm nx_plot --image_type pdf  -x ${ftype} \
	${inPATH}/${magDIR}/ ${inPATH}/${magDIR}/${magQA}/${plotsDIR}/
fi

if [ ${krona} == T ]; then 
TAXbin=TAXprofiles/${magDIR}TAX_checkm
if [[ ! -d ${TAXbin} ]]; then mkdir -p ${TAXbin}/bin; fi
cat ${list}.txt | xargs -I {} sh -c \
"join ${inPATH}/${magDIR}/${magQA}/${magfix}_profiles.txt ${inPATH}/${magDIR}/${magQA}/${magfix}_tax.txt -t $'\t' | \
	cut -f6,10,11 | sed -e 's/..__/\t/g' |sed -e 's/\tunresolved//g' | sed -e '1 s/^/#/' \
	 > ${TAXbin}/bin/{}_4krona.txt"
ktImportText ${TAXbin}/bin/*_4krona.txt -o ${TAXbin}/TAXplots_${magDIR}.html
fi

echo "---> MetaGenome Assembled Genomes (draftGenomes) analysis done. Nazdrave"
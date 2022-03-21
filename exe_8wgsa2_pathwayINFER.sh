#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=16g
#SBATCH --job-name pathInfer
#SBATCH --mail-type END
#SBATCH --output slogs/slog_%x.txt
#SBATCH --time=8:00:00
#SBATCH --partition=norm #quick  #### quick, but walltime needs to be 4hr limited
if [ ! -d "slogs" ]; then mkdir slogs; fi

module load kronatools
module load python/2.7 #pretty sure it needs 2.7, not >3
#also genes2KronaTable.py has some sort of interference with eggnog-mapper module

list=${1%%.*}
mapping=${2} 
TPMfile=pregenes_ABUNtab.txt
outdir=pathways						
supScriptPATH=/home/angelovaag/MyScripts/Nephele_WmGS_v2/v2.1_opt
inPATH=asmbMetaSpades/{}_asmb/genes

	if [[ "${mapping}" == "" ]] || [[ "${mapping}" == "KEGG" ]]; then
		echo "------> mapping KO numbers to ${mapping} pathways"
		export 	 map=${supScriptPATH}/MinPath-master/data/KEGG_myfiles_2021/KEGG-ko2keggpath.txt
		export   hrr=${supScriptPATH}/MinPath-master/data/KEGG_myfiles_2021/KEGGBrite-hierarchy-me.txt
			  prefix=ko2gg
			  ANNOTSfile=annotations/annots.ko.txt
			  mapping=KEGG
			  PWYdir=PWYprofiles/keggPWY/							#bucket bin for PWY krona file collection
	else
		echo "------> mapping EC numbers to ${mapping} pathways"
			#my EC/GO->MetaCyc mapping & hierarchy files for MinPath & genes2pathway (created Feb 2021)
			#Code for EC ---> MetaCyc pathway transformations (EC#s from EggNog-mapper or Metaprokka)
		export	map=${supScriptPATH}/MinPath-master/data/MetaCyc_myfiles_2021/MetaCyc-ec2path-me.txt 
		export  hrr=${supScriptPATH}/MinPath-master/data/MetaCyc_myfiles_2021/MetaCyc-hierarchy-me.txt
			 prefix=ec2mc
			 ANNOTSfile=annotations/annots.ec.txt
			 PWYdir=PWYprofiles/metacycPWY/							#bucket bin for PWY krona file collection
	fi


cat ${list}.txt |xargs -I {} bash -c \
"if [ ! -d ${inPATH}/${outdir}/ ]; then mkdir ${inPATH}/${outdir}/ ; fi"  # if dir exists, exit code=1 => cannot use &&

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
		cp ${inPATH}/${ANNOTSfile} ${inPATH}/${outdir}/tmp_annots.txt '&&' \
		cut -f 1,5 ${inPATH}/${TPMfile} '>' ${inPATH}/${outdir}/tmp_TPM.txt 

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	cd ${inPATH}/${outdir}/ '&&' \
	python ${supScriptPATH}/MinPath-master/MinPath.py -any tmp_annots.txt -map ${map}  \
		-report  tmp_report.txt -details ${prefix}.details.txt '>>' ${prefix}.minpath.log.txt 

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	cd ${inPATH}/${outdir}/ '&&' \
 	grep 'minpath\ 1' tmp_report.txt '>' ${prefix}.report.txt 

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	cd ${inPATH}/${outdir}/ '&&'  \
	python ${supScriptPATH}/genes2KronaTable.py -i tmp_annots.txt \
		-m ${map} -H ${hrr} -n {} -c tmp_TPM.txt \
		-l ${prefix}.report.txt \
		-o tmp_4kr.txt 

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	cd ${inPATH}/${outdir}/ '&&' \
	sed '1s/^/#/' tmp_4kr.txt '>' {}_${prefix}_4krona.txt

parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	rm ${inPATH}/${outdir}/tmp_*.txt #\


if [[ ! -d "${PWYdir}" ]];then mkdir -p ${PWYdir}/bin ; fi
parallel -a ${list}.txt -j $SLURM_CPUS_PER_TASK \
	cp ${inPATH}/${outdir}/{}_${prefix}_4krona.txt ${PWYdir}/bin
ktImportText ${PWYdir}/*${prefix}_4krona.txt -o ${PWYdir}/PWYplots_${mapping}.html


cat ${list}.txt |xargs -I {} sh -c \
"if [[ ! -s ${inPATH}/${outdir}/{}*_4krona.txt ]]; then echo {} >> failed_samples.txt; fi "


echo "Pathway Inference with MinPATH done"

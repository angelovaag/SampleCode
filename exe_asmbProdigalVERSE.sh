#!/bin/bash
#SBATCH --cpus-per-task=42
#SBATCH --mem-per-cpu=2g
#SBATCH --job-name asmbPrdgl
#SBATCH --mail-type END
#SBATCH --output slurm_prodigallog.txt
#SBATCH --time=8:00:00
#SBATCH --partition=norm  #### quick/norm, quick walltime must be < 4hr 

module load prodigal

indir=asmbMetaSpades
prefix=assembly.genes # ${prefix}
type=gff
infile=final.assembly #needs both .bam and .fasta
outdir=predictedGenes
#ted_path=trimNdecon
readlen=150



for f in $(ls -d ${indir}/*_asmb/ |cut -f1-3 -d "_" ); do
path=${f}_asmb

if [ ! -d "${path}/${outdir}"  ]; then
	mkdir ${path}/${outdir}
fi

echo "
	----------------------
	running prodigal for " ${f} "full assembly"
prodigal -i ${path}/${infile}.fasta -o ${path}/${outdir}/${prefix}.${type} -f ${type} -q \
	-a ${path}/${outdir}/${prefix}.faa -d ${path}/${outdir}/${prefix}.fna -p meta

echo "----------------------
	converting gff3 to gtf for" ${f}
~/MyScripts/Nephele_WmGS_pipeline/prodigalgff2gtf.sh ${path}/${outdir}/${prefix}.gff > ${path}/${outdir}/${prefix}.gtf

#running VERSE to count reads per gene (input's bam!)
echo "----------------------
	getting gene counts with VERSE" ${f}
~/MyTools/verse -a ${path}/${outdir}/${prefix}.gtf -t 'CDS' -g gene_id -z 0 -s 0 \
		-o ${path}/${outdir}/${prefix}_counts     ${path}/${infile}.bam -T $SLURM_CPUS_PER_TASK

#running TPM normalization from the counts
echo "----------------------
	getting gene coverage & TPM normalization for sample" $(basename ${f})
cut -f4,5,9 ${path}/${outdir}/${prefix}.gtf | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | \
		tr ' ' '\t' | sed '1s/^/gene\tlength\n/' > ${path}/${outdir}/${prefix}_lengths.txt
join ${path}/${outdir}/${prefix}_lengths.txt ${path}/${outdir}/${prefix}_counts.CDS.txt -t $'\t' >  ${path}/${outdir}/${prefix}_counts.txt

grep "_" ${path}/${outdir}/${prefix}_counts.txt | awk -v OFS="\t" '{$5 = sprintf("%0.0f", $3*150/$2)}1' |sort | cut -f 1,2,3,5 | \
 		sed -e '1s/^/#name\tlength\tNumbAlignedReads\trelGENE_coverage\n/' > ${path}/${outdir}/${prefix}_gene.coverage.txt

awk 'NR==FNR{sum+= $4; next} FNR==1{print $0,"relGene_TPM"; next} {printf("%s %0.0f\n",$0,$4*1000000/sum)}' OFS='\t' \
  		${path}/${outdir}/${prefix}_gene.coverage.txt ${path}/${outdir}/${prefix}_gene.coverage.txt |sed -e 's/\s/\t/g' \
  		> ${path}/${outdir}/${prefix}_gene.TPM.txt


done
echo "prodigal.sh done"



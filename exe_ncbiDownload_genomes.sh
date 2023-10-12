#!/bin/bash
# DLs any genome from ncbi

# From https://www.ncbi.nlm.nih.gov/datasets/genome, 
# find assembly of interest and grab the "NCBI RefSeq assembly", thats $1

acc=${1}  #  GCF_000001405.40 or GCF_003339765.1: the assembly accession # (NCBI RefSeq assembly)
TAXid=${2}
dir=$(pwd)
annotation_type=GENOME_FASTA
#other annotation_types: GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT



curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${acc}/download?include_annotation_type=${annotation_type}&filename=${acc}.zip" -H "Accept: application/zip"

unzip ${acc}
mv ncbi_dataset/data/${acc}/${acc}*_genomic.fna ${dir}
rm -rf ncbi_dataset/ ${acc}.zip README.md

if [[ ! ${TAXid} == ""  ]]; then
	# mv ${acc}*_genomic.fna ${acc}_tid${TAXid}.fna
	rename "_genomic.fna" "_txid${TAXid}.fna" ${acc}*_genomic.fna
fi
echo "---> downloaded genome is named:" ${acc}*_txid${TAXid}.fna


sed -i "s/^>/>kraken:taxid\|${TAXid}\|/g"   ${acc}*_${TAXid}.fna
echo "---> added kraken header with taxID" ${TAXid}

# and the ftp site is: 
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz 
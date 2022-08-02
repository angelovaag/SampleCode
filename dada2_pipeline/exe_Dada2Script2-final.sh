#/usr/bin/env Rscript
# local script for L444
#Usage: edit these parameters to what is needed, run script as executable, "> log.txt"
logfile=Rlog.txt
wdir=${1}
indirName=${2}
outdirName=dada2

#module load R ### -> should not need loading with /usr/bin/env Rscript
# For the belkaid lab, the leftTrim 20bp, is unnecessary as the V4 primer is already trimmed
Rscript ~/Documents/MyApps/RScripts/dada2_pipeline/Dada2Script2-final.R\
 --wdir  ${wdir}\
 --path ${wdir}/${indirName}\
 --outdir ${wdir}/${outdirName}/\
 --filtpath ${wdir}/${outdirName}/flt\
 --fwdptn _R1_\
 --revptn _R2_\
 --taxDB ~/Documents/DB-rRNA/SILVA/silva_nr_v138_train_set-DADA2.fa.gz\
 --taxDB_species ~/Documents/DB-rRNA/SILVA/silva_species_assignment_v138-DADA2.fa.gz\
 --decipherDB ~/Documents/DB-rRNA/SILVA/SILVA_SSU_r138_2019.RData\
 --trimLeft 0\
 --nthreads 8\
 --nrds 1e5\
 --trimRight 0\
 --truncLen 0\
 --maxEE 5\
 --trunQ 4\
 --minLen 50

# def trimLeft = 20, but not for MicrobiomeCore V4 


#mkdir -p ${path}/${outdirName}/filtered/
#mv ${path}/${outdirName}/*fastq.gz ${path}/${outdirName}/filtered/

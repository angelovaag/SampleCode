#/usr/bin/env Rscript
# local script for L444
#Usage: edit these parameters to what is needed, run script as executable, "> log.txt"
logfile=dada2_Rlog.txt
projDir=$(pwd)
    echo $projDir
RAWsDir=${1}
outDir=${2} # dada2_V4seq_trimLeft0
            if [[ ${outDir} == "" ]]; then outDir=dada2_out; fi
sampName_segments=${3} #number of underscored segments in sampleNames
            if [[ ${sampName_segments} == "" ]]; then sampName_segments=1; fi
LeftTrim=0 #[default is 20, unless Belkaid V4seq]

DB=SILVA
if [[ ${DB} == "SILVA" ]]; then
taxDB=~/Documents/DB-rRNA/SILVA/SILVA138.1/silva_nr99_v138.1_train_set-DADA2.fa.gz
taxDB_spp=~/Documents/DB-rRNA/SILVA/SILVA138.1/silva_species_assignment_v138.1-DADA2.fa.gz
fi
if [[ ${DB} == "HOMD" ]]; then
taxDB=~/Documents/DB-rRNA/HOMD/v15.23_Apr2023/HOMD_16S_rRNA_RefSeq_V15.23.dada2.train_set.fa
taxDB_spp=~/Documents/DB-rRNA/HOMD/v15.23_Apr2023/HOMD_16S_rRNA_RefSeq_V15.23.dada2.species_assignment.fa
fi

taxDB_DphR=NULL # this will omit decipher as it takes a long time and does not produce better assignment
# taxDB_DphR=~/Documents/DB-rRNA/SILVA/SILVA138.1/SILVA_SSU_r138_2019-DECIPHER.RData

# For V4 sequences from the Belkaid lab, the leftTrim = 0 bp, as V4 primer is already trimmed
# For any other sequencing, set leftTrim = 20

Rscript ~/Documents/MyApps/RScripts/dada2_pipeline/dada2_script-local_v2.R \
 --projDir  ${projDir}\
 --path ${projDir}/${RAWsDir}\
 --outdir ${projDir}/${outDir}/\
 --filtpath ${projDir}/${outDir}/flt \
 --fwdptn _R1_ --revptn _R2_  \
 --taxDB ${taxDB} --taxDB_species ${taxDB_spp} \
--decipherDB ${taxDB_DphR} \
 --nseg ${sampName_segments}  --nthreads 8 \
 --trimLeft ${LeftTrim} --trimRight 0  \
 --nrds 1e3 --truncLen 0 --maxEE 5 --trunQ 4 --minLen 50

# def trimLeft = 20, but not for MicrobiomeCore V4 

## Just in case for troubleshooting:
cp ~/Documents/MyApps/RScripts/dada2_pipeline/dada2_script-local_v2.R ${outDir}/Rfiles/
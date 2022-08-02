library(docopt)
"Usage: my_program.R [options]
--help           show this screen
--wdir PATH      give path to working dir, [default: '~/Documents/TestData/subsampled-10K/']
--path PATH      give path to demultiplexed fastq files [default: '~/Documents/TestData/subsampled-10K/']
--filtpath PATH  path of a sub-directory to create for filtered fastq files [default: '~/Documents/TestData/subsampled-10K/output']
--fwdptn PTN     type in pattern for recognition of R1-fastq.gz files. [default: _R1_subsampled.fastq.gz]
--revptn PTN     type in pattern for recogniztion of R2-fastq.gz files. [default: _R2_subsampled.fastq.gz]
--outdir NAME    path to & name of output directory to create [default: ~/Documents/TestData/output]
--taxDB FILE     filname of DB [default: '~/DB-rRNA/RDP/rdp_train_set_16-DADA2.fa.gz']
--taxDB_species FILE filneme of DB with species [default: '~/DB-rRNA/RDP/rdp_species_assignment_16-DADA2.fa.gz']
--decipherDB FILE path & filename of the DECIPHER RData file.[default: '~/Documents/DB-rRNA/RDP/RDP_v16-mod_March2018.RData']
--nthreads NUMB  number of threads [default: 4] half the ones from interractive sessio
--nrds NUMB      number of reads to read in into memory [default: 1e5]
--trimLeft NUMB  number of bases to trim on the left [default: 20]
--trimRight NUMB number of bases to trim on the right [default: 0]
--truncLen NUMB  number of bases to truncate reads to [default: 0]
--maxEE NUMB     discard reads with expected errors bellow this number [default: 5]
--trunQ NUMB     remove reads with quality below this number [default: 4]
--minLen NUMB    filter out reads of length below this number [default: 50]" -> doc

docopt(doc)
args <-docopt((doc))#, args=c(
#  "--wdir", "~/Documents/TestData/subsampled-10K", 
#  "--outdir", "~/Documents/TestData/subsampled-10K/output",
#  "--path", "~/Documents/TestData/subsampled-10K",
#  "--filtpath", "~/Documents/TestData/subsampled-10K/output", 
#  "--fwdptn", "_R1_subsampled.fastq.gz", 
#  "--revptn", "_R2_subsampled.fastq.gz", 
#  "--taxDB", "~/Documents/DB-rRNA/RDP/rdp_train_set_16-DADA2.fa.gz", 
#  "--taxDB_species", "~/Documents/DB-rRNA/RDP/rdp_species_assignment_16-DADA2.fa.gz",
#  "--decipherDB", "~/Documents/DB-rRNA/RDP/RDP_v16-mod_March2018.RData",
#  "--trimLeft", "20",
#  "--nthreads", "2",
#  "--nrds", "1e3",
#  "--trimRight", "0",
#  "--truncLen", "0",
#  "--maxEE", "5",
#  "--trunQ", "4",
#  "--minLen", "50"))
args


# set workdir}
setwd(args$wdir)

#create and set output directory
print(paste("creating an output directory at", args$wdir))
dir.create(args$outdir)
outdir=args$outdir


# load packages}
library(ggplot2)
library(readr)
library(reshape2)
library(tidyr)
library(dada2)
library(biomformat)
library(DECIPHER)
# load paths and define data}
#"import" data
print(paste("importing .fastq.gz files from", args$path))
path <- args$path
filtpath <- args$filtpath
#file-patterns:
fwd.ptn=args$fwdptn
rev.ptn=args$revptn

fwd <- list.files(path, pattern=fwd.ptn) # CHANGE if different file extensions
rev <- list.files(path, pattern=rev.ptn) # CHANGE if different file extensions

#r parameters}
#threads
nthread= as.numeric(args$nthreads) #DO NOT EXCEED HALF THE NUMBER OF CPUS PROVIDED FOR YOUR interractive session

#taxonomy
taxDB=args$taxDB
taxDB_species=args$taxDB_species
load(args$decipherDB) #for the decipher training set

#for primer trimming (def = 0). Truncation is performed AFTER trimming. Trimming is first
trimLeft = as.numeric(args$trimLeft)
trimRight = as.numeric(args$trimRight)

#for truncation, minLen and exp.Error filtering 
truncLen= as.numeric(args$truncLen) #truncate all reads after # of bases (0 means no truncation); or for PE, do c(#,#) for F&R trunc
maxEE= as.numeric(args$maxEE) #discard reads with expected errors below maxEE
truncQ= as.numeric(args$trunQ) #remove reads at first qual below truncQ
minLen= as.numeric(args$minLen)



#r inspect quality and make qc plots}
print("inspecting quality")
readqc.fwd=plotQualityProfile(file.path(path,fwd))
readqc.rev=plotQualityProfile(file.path(path,rev))

qplotheight <- min(2.5*ceiling(length(fwd)/4), 49)

pdf(paste0(outdir, "qualityProfile-fwd.pdf"), width=16, height=qplotheight); readqc.fwd; dev.off()
pdf(paste0(outdir, "qualityProfile-rev.pdf"), width=16, height=qplotheight); readqc.rev; dev.off()



#r filter and trim}
# Filtering
print("filtering and trimming")
filtered=filterAndTrim(file.path(path,fwd), file.path(filtpath,fwd), 
                       file.path(path,rev), file.path(filtpath,rev),
                       compress=TRUE, verbose=TRUE, rm.phix=TRUE,
                       truncLen=as.numeric(truncLen), maxEE=as.numeric(maxEE),
                       truncQ=as.numeric(truncQ), trimLeft = as.numeric(trimLeft),
                       minLen=as.numeric(minLen), multithread=as.numeric(nthread),
                       n=as.numeric(args$nrds)) 
#head(filtered)



#r link filename to sample name with function}
#for forward reads
flts.fwd <- list.files(filtpath, pattern = fwd.ptn, full.names = T)
flts.fwd
sample.names.fwd<-sapply(strsplit(basename(flts.fwd), "_"), 
                         function (x) paste(x[1:2], collapse="_"), simplify=F)

#for the reverse
flts.rev <- list.files(filtpath, pattern = rev.ptn, full.names = T)
flts.rev
sample.names.rev<-sapply(strsplit(basename(flts.rev), "_"), 
                         function (x) paste(x[1:2], collapse="_"), simplify=F)


#check if names for F and R are identical
if(!identical(sample.names.fwd, sample.names.rev)) stop("Forward and reverse files do not match.")
#if no errors, continue to link the sample names to the fastq file names
names(flts.fwd) <- sample.names.fwd; sample.names.fwd
names(flts.rev) <- sample.names.rev; sample.names.rev



#r learn the error rates, sample inference, merger of PE reads}
print("learning error rates")
set.seed(100)
# Learn forward error rates
errF <- learnErrors(flts.fwd, nbases=1e8, multithread=nthread)
# Learn reverse error rates
errR <- learnErrors(flts.rev, nbases=1e8, multithread=nthread)

#Plot error rates
errF.plot=plotErrors(errF, nominalQ = T); errF.plot
errR.plot=plotErrors(errR, nominalQ = T); errR.plot
pdf(paste0(outdir, "errRate_FWD.pdf")); errF.plot; dev.off()
pdf(paste0(outdir, "errRate_REV.pdf")); errR.plot; dev.off()


#r Sample inference and merging}
#Sample inference & merging
print("merging reads")
mergers<- vector("list", length(sample.names.fwd))
names(mergers)=sample.names.fwd
#mergers

for(sam in sample.names.fwd) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(flts.fwd[[sam]])
  ddF <- dada(derepF, err=errF, multithread=nthread)
  derepR <- derepFastq(flts.rev[[sam]])
  ddR <- dada(derepR, err=errR, multithread=nthread)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
#rm(derepF); rm(derepR)

#Construct seq table
print("contstrucing seq table")
seqtab<- makeSequenceTable(mergers)
saveRDS(seqtab, paste0(outdir,"seqtab.rds"))



#r remove chimeras, assign Taxonomy}
#merge multiple runs (if necessary)
#seqtab1<-readRDS("output/seqtab.rds") # to read back the chimeric ASV
#seqtab2<-readsRDS("ouptu/seqtab1.rds") #and so on
#seqtab<-mergeSequenceTable(seqtab1, seqtab2)

#remove chimeras
print("removing chimeras")
seqtab<-removeBimeraDenovo(seqtab, method="consensus", multithread=nthread)

#explore some output
print("explore some output")
  #head(ddF); head(ddR)
  #head(mergers[[1]])
# Inspect distribution of sequence lengths
  table(nchar(getSequences(seqtab)))
#track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered, sapply(mergers, getN), rowSums(seqtab))#sapply(ddF, getN), sapply(ddR, getN),
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "nonchim") # "denoisedF", "denoisedR",
rownames(track) <- sample.names.fwd
  #head(track)
write.table(track, paste0(outdir, "readtracker.txt"), sep="\t", col.names = NA)

#assign Taxonomy
print("assigning taxonomy")
tax<-assignTaxonomy(seqtab, taxDB, multithread = nthread); #head(tax)
tax<-addSpecies(tax, taxDB_species); #head(tax)

#assign Taxonomy with DECIPHER
print("and now with DECIPHER")
library(DECIPHER)
dna<-DNAStringSet(getSequences(seqtab)); head(dna) #creates a DNA string set
ids<-IdTaxa(dna, trainingSet, strand="both", processors = nthread, verbose=F)
ranks<-c("domain", "phylum", "class", "order", "family", "genus", "species")
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)
#head(taxid)


#export
print("export some data")
saveRDS(seqtab, paste0(outdir, "seqtab-nonChimASVs.rds"))
saveRDS(taxid, paste0(outdir, "tax.rds"))



#export final data}
print("export some more data")
##creating a COUNTS table (no taxonomy)
seqtab.print<-(seqtab)
colnames(seqtab.print)<-paste0("ASV", 1:ncol(seqtab))
head(seqtab.print)
write.table(t(seqtab.print), paste0(outdir, "ASV_counts.txt"), sep="\t", na="", col.names=NA, quote=F)


##creating tab-sep tax table (from assignTaxonomy)
tax.print<-tax
row.names(tax.print)<-paste0("ASV", 1:nrow(tax)) # <-NULL raw.names are the entire ASV sequence => replacing
head(tax.print)
write.table(tax.print, paste0(outdir, "ASV_taxtable.txt"), sep="\t", col.names=NA, quote = F, na="")

#ASV_fullTable (contains counts+taxonomies in 1)
asvtab=as.data.frame(cbind(t(seqtab.print), tax.print))
#head(asvtab)
write.table(asvtab, paste0(outdir, "ASV_fulltable.txt"), sep="\t", col.names=NA, quote=F)
write_biom(asvtab, paste0(outdir, "ASV_fulltable.biom"))

#creating tab-separated tax table (from DECIPHER)
taxid.print<-taxid
row.names(taxid.print)<-paste0("ASV", 1:nrow(taxid)) # <-NULL raw.names are the entire ASV sequence => replacing
#head(taxid.print)
write.table(as.data.frame(taxid.print), paste0(outdir, "ASV_taxIDtable2.txt"), sep="\t", col.names=NA, quote=F)


##creating a .fna of the ASVs
seq.print<-seqtab
seq.fna=paste(paste0(">ASV", 1:ncol(seq.print)), colnames(seq.print), sep="\n")
#head(seq.fna)
write.table(seq.fna, paste0(outdir, "ASV_seqs2.fna"), col.names=F, row.names = F, quote = F)


#r save workspace}
#save.image("Dada2Script1-Mrkdown.RData")
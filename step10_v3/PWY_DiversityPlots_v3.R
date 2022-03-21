#! /usr/bin/env Rscript
### load packages
## info fo fossil/chao1: https://www.rdocumentation.org/packages/fossil/versions/0.4.0/topics/chao1
# other packages not needed for this script
#y=c("ape", "broom","readr", "reshape2", "dplyr", "tidyr", "tidyverse",
#"Biostrings","devtools","tibble", "stringr", "htmltools")
library(easypackages) #can also install multiples with packages()
x<-c("data.table", "fossil", "ggplot2", "vegan", "ampvis2" , "tidyverse", "S4Vectors")
libraries(x);# libraries(y)
rm(x)

print("------> start of PWY Diversity Plots <------")
print("------> docopt settings ------")
library(docopt)
### prep arguments for docopt #
#(use only one space between --argument & PATH/NAME, ok to have multiple  after  or TABs)
doc <- 'Usage: DiversityPlots_v3.R [options]
  --help            show this screen.
  --wdir PATH       give path to working dir /project main dir  [default: pwd ].
  --indir PATH      give path to input files directory [default: profiles/merged_tables ].
  --mfile NAME      give name of mapping file [default: mapping.txt ].
  --outdir NAME     give path & name of output directory to create [default: pwd/profiles/DiversityPlots].' # -> doc
doc
args <- docopt(doc); 
print(args)
# ## manual args set:
# args <- docopt(doc, args = c(
#    "--wdir"   ,  "/Users/angelovaag/Documents/Nephele_pipelines/WmGS_v2/tests/filesFromDuc",
#    "--indir"  ,  "/Users/angelovaag/Documents/Nephele_pipelines/WmGS_v2/tests/filesFromDuc/merged_tables_PWY",
#    "--mfile"  ,  "/Users/angelovaag/Documents/Nephele_pipelines/WmGS_v2/tests/filesFromDuc/mapping.txt",
#    "--outdir" ,  "/Users/angelovaag/Documents/Nephele_pipelines/WmGS_v2/tests/filesFromDuc/PWYout")) #this is if the .sh does not provide arguments
#  print(args)
print("-----------> docopt settings------^^^^")

# set workdir
setwd(args$wdir)

# #load("PWY_DiversityPlots_v3.RData")
# # 
# # #=============== starting tab input ==============
path=list()
path$input = paste0(args$wdir)
path$tables= paste0(args$indir)
path$outdir=paste0(args$outdir); dir.create(path$outdir); 

# #importing metadata objects
tables=NA
tables$meta <- read.table(paste0(args$mfile), header = T, sep="\t", quote = "");
row.names(tables$meta)= tables$meta$SampleID # tables$meta=tables$meta %>% relocate(SampleID)
tables$meta=tables$meta[order(tables$meta$SampleID), ]
tables$meta; print("-----------> read in mapping file^^^")

### Beta diversity with Ampvis2
#https://madsalbertsen.github.io/ampvis2/articles/ampvis2.html
print("-----------> prepping AMPVIS2 object")
#tables$amp=cbind(tables$counts, tables$tax) # or 
tables$amp=read.table(paste0(path$tables, "/merged_Counts+PWY.txt"),  header = T, row.names = 1, sep = "\t", quote = "" )
tables$amp = tables$amp %>% dplyr::rename(Species = Tier4, Family = Tier3 , Class = Tier2, Phylum = Tier1) %>%
  dplyr::select(-c(starts_with("Tier") ) ) #bc AWS::ampvis2.7.17 makes issue with the extra columns from metaCyc DB
head(tables$amp); tail(tables$amp)
# tables$amp = tables$amp %>% S4Vectors::rename("Tier4" = "Species", "Tier3" = "Family", "Tier2" = "Class", "Tier1" = "Phylum")
# try({
#   print("--------> attempting removal of extra metaCyc columns (existing only for metaCYC pwys)")
#   # remtiers = paste0(rep("Tier"), seq(5,14,1))
#   # tables$amp<- dplyr::select(tables$amp, -c(any_of(remtiers) ) ) #bc AWS::ampvis2.7.17 makes issue with the extra columns from metaCyc DB
#   })

amp=amp_load(tables$amp, tables$meta)
# print("-----------> checking some AMPVIS2 data")
# head(amp$abund); tail(amp$abund)
# head(amp$tax); tail(amp$tax)
# head(amp$metadata)
# print("-----------> checking AMPVIS2 data ^^^")

print("-----------> making PCOA plot")
try({
pcoa_bray <- amp_ordinate(amp, filter_species = 0.01, type = "PCOA", sample_label_by = "SampleID", sample_label_size = 3,
                          distmeasure = "bray", sample_color_by = "SampleID", sample_point_size = 4, opacity=1.5,
                          detailed_output = TRUE, transform = "none", sample_shape_by = "TreatmentGroup")
#pcoa_bray$plot= pcoa_bray$plot + scale_shape_manual(values=c(17, 15,13,18, 4,5,6,7,8) )
pdf(paste0(path$outdir, "/PWY_BetaDiv_PCoA.pdf"), width=12, height=10); plot(pcoa_bray$plot); dev.off()
})

print("-----------> making nMDS plot")
try({
nmds_bray <- amp_ordinate(amp, filter_species = 0.01, type = "NMDS", sample_label_by = "SampleID", sample_label_size = 3,
                          distmeasure = "bray", sample_color_by = "SampleID", sample_point_size = 4, opacity = 1.5,
                          detailed_output = TRUE, transform = "none", sample_shape_by = "TreatmentGroup")
#nmds_bray$plot=nmds_bray$plot + scale_shape_manual(values=c(17, 15,13,18, 4,5,6,7,8))
pdf(paste0(path$outdir, "/PWY_BetaDiv_nMDS.pdf"), width=12, height=10); plot(nmds_bray$plot); dev.off()
})

print("-----------> making heatmap profile graph")
try({
amp_htmp=amp_heatmap(amp, group_by = "SampleID", facet_by = "TreatmentGroup",
                     tax_aggregate = "Species", tax_show=35, plot_values = T,  normalise=T,
                     plot_values_size = 3, showRemainingTaxa = F,  tax_add = "Class") + 
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8), legend.position="right")
pdf(paste0(path$outdir, "/PWY_Profile_Heatmap_top35PWYs.pdf"), width=14, height = 10); plot(amp_htmp); dev.off()
})


#cleanup
rm(amp_htmp, nmds_bray, pcoa_bray)


#save R image
save.image("PWY_DiversityPlots_v3.RData")

#Attempting DiffAbund: not all dataset will have replicates & groups rich enough for statistical analysis
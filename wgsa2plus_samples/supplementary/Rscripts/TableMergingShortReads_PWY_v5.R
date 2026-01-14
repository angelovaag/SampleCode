# load packages
library("easypackages") # can also install multiples with packages()
x<-c("tidyverse", "reshape2", "broom", "data.table", "tibble", "docopt")
suppressMessages(libraries(x)); rm(x)

print("------> start of PWY TABLE COLLATION <------")
# library(docopt)
print("------> docopt settings ------")
doc <- "Usage: my_program.R [options] 
 --help              show this screen
 --binDIR PATH       hard path to bin dir [default: 'wdir/TEDreadsTAX/bin']
 --sccList NAME      hard path and name of file listing successful sample names [default: 'SccList.txt']
 --outdir NAME       path to & name of output directory to create [default: 'wdir/TEDreadsTAX/merged_tables']
 --genesDIR PATH     hard path to bin dir for KO or EC genes [default: 'wdir/PWYprofiles/keggPWY/genebin']"
print(doc) %>% cat()

args <- docopt(doc)
print(args)
# for manual prep #this is if the .sh does not provide arguments
# args <-docopt(doc, args = c(
#  # "--wdir"   , "/Users/angelovaag/Documents/tests/table_merging_tests/ ",                             # "/path/to/WmGS_v2/tests/" ,
#  "--binDIR" , "/Users/angelovaag/Documents/tests/table_merging_test/binDIR_ec/",                    # "/path/to/PWYprofiles/bin/",
#  "--sccList", "/Users/angelovaag/Documents/tests/table_merging_test/scclist.txt",               # "/path/to/SccList.txt",
#  "--outdir" , "/Users/angelovaag/Documents/tests/table_merging_test/merged_tables_v5_ec/",          # "/path/to/PWYprofiles/merged_tables/" )
#  "--genesDIR",  "/Users/angelovaag/Documents/tests/table_merging_test/genebin/") )
# # args
# print("-----------> docopt settings------^^^^")
# load("PWY_TableMerging_v3.RData")
# 
# #=============== finished ==============
input=list() #instead of input=NA
input$path=args$binDIR    #"~/Documents/DataTest/CAMI_HiSeq_data_wgs/4PWYonomic_plots/"
input$files=list.files(input$path, full.names = T, pattern = ".txt") 
# print(paste("detecting input files from ", input$path) )
# input$files

input$scclist <- read.table(args$sccList, header = F, sep="\t") %>% 
            setNames("SampNames") %>%  arrange(SampNames) %>% 
            rename_with(~ str_remove(.,".txt"), contains(".txt")) %>% 
            rename_with(~ str_remove(.,"_4krona"), contains("_4krona"))
names(input$files) <- input$scclist$SampNames;
print("------> naming input files:")
input$files

#load data from files
print("------> load data from files")
tiers=paste0(rep("Tier"), seq(1,18,1)) #assuming max(metaCyc PWY Tiers) = 18 (its actually 14)
data=list()
data=lapply(input$file, function(x){
  y<-read.table(x, header = FALSE, sep = "\t", na.strings = c("NA", "", " "), fill = TRUE,
                quote = "", check.names = FALSE, col.names = c("Counts", tiers) ) %>% remove_rownames() %>%
    group_by(across(any_of(tiers))) %>%  relocate("Counts") %>% filter(!Counts %in% "X") %>%
    unite("allTiers", any_of(tiers), sep="; ", remove=T, na.rm=T) %>%
    relocate("allTiers", .after=last_col()) %>%
    mutate(PWYname = gsub(".*; ", "", allTiers) ) %>% 
    mutate(Counts = round(as.numeric(Counts), digits = 2))
})

# lastcol <- dim(data[[1]])[[2]]; lastrow<- dim(data[[1]])[[1]]
# data[[1]][1:10, c(1:lastcol)]
# data[[1]][134:144, c(1:lastcol)]


print("------> merging PWY data")
merged0=Reduce(function(x,y) merge(x,y, by=c("allTiers","PWYname"), all=T ), data)
# head(merged0)
merged <- merged0 %>% setNames(c("allTiers", "PWYname", names(data))) %>% 
  mutate(rn = PWYname) %>% column_to_rownames("rn") %>%
  separate(allTiers, into = tiers, sep = "; ", remove = F) %>%
  discard(~all(is.na(.x))) %>%
  mutate_if(is.numeric, ~replace(. , is.na(.), 0)) %>%
  relocate("allTiers", .after=last_col())
# View(merged)

# print("------> checking merged object:")
# head(merged); tail(merged); dim(merged);

output=list()
output$dir<- args$outdir
print(paste("------> exporting PWY data at:", output$dir))
dir.create(output$dir)

write.table(merged, file.path(output$dir, "merged_Counts+PWY+allTiers.txt"),
            sep="\t", col.names=NA, quote = F, na="")

noHier <- setdiff(names(merged), c("allTiers"))
write.table(merged[,noHier], file.path(output$dir, "merged_Counts+PWY.txt"),
            sep="\t", col.names=NA, quote = F, na="")

Lin <- select(merged, "allTiers")
write.table(Lin, file.path(output$dir, "merged_allTiers.txt"),
            sep="\t", col.names=NA, quote = F, na="")

noTAX <- setdiff(names(merged), tiers)
noTAX <- merged[, noTAX] %>% relocate(allTiers, .after=last_col() )
write.table(noTAX, file.path(output$dir, "merged_Counts+allTiers.txt"),
            sep="\t", col.names=TRUE, row.names = F, quote = F, na="")

write.table(merged[,names(data)], file.path(output$dir, "merged_Counts.txt"),
            sep="\t", col.names=NA, quote = F, na="")

noCounts<-setdiff(names(merged), c(names(data), "allTiers") )
write.table(merged[,noCounts], file.path(output$dir, "merged_PWY.txt"),
            sep="\t", col.names=NA, quote = F, na="")

# bim <-biomformat::make_biom(merged[,names(data)], observation_metadata = merged[,noCounts])
# biomformat::write_biom(bim, biom_file = paste0(output$dir, "merged_Counts+PWY_json.biom") )
# print("WARNING: please import biom with read_biom(), not with phyloseq, qiime or ampvis biom-import features")
# test biom file:
# t <-read_biom(paste0(output$dir, "merged_Counts+TAX_json.biom")) #imports
# t <- amp_import_biom(paste0(output$dir, "merged_Counts+TAX_json.biom")) #does not import
# rm(t)



if(args$genesDIR == "NA" ){ print("---> no geneTABLES to merge. Disabled.")}else{
            ######## Code for merging geneTPM tables (KO and EC)
          print("------> Starting work on genes files:")
          ingenes=list()
          ingenes$files=list.files(args$genesDIR, full.names = T, pattern = "TPM") 
          ingenes$names=read.table(args$sccList, header = F) %>% setNames("SampNames") %>% arrange(SampNames)
          names(ingenes$files) <- ingenes$names$SampNames ;


          #load data from files
          print("------> load GENE data from gene files")
          genedata=list()
          genedata=lapply(ingenes$files, function(x){
            y<-read.table(x, header=F, check.names = F, sep="\t", na.strings =c("NA",""," "), row.names = NULL,
                          fill=T, quote = "", col.names = c("annot", "geneTPM", "geneNAME") ) %>%
              group_by(across(any_of("annot"))) 
          })
          # names(genedata)

          print("------> merging genedata, setting names to samples, replacing NAs with 0s, and removing rowSum==0s")
          genemerged=Reduce(function(x,y)  merge(x,y, by=c("annot", "geneNAME"), all=T ), genedata)
          genemerged = genemerged %>% setNames(c("annot", "geneNAME", names(genedata)) ) %>% 
            mutate_if(is.numeric, ~replace(. , is.na(.), 0)) %>% 
            column_to_rownames("annot") %>% filter(rowSums(select_if(., is.numeric)) !=0)
          # print("------> checking genemerged object:")
          # head(genemerged); tail(genemerged); dim(genemerged);

          print("------> exporting GENE data")
          # output=list()
          # output$dir<- args$outdir
          # print(paste("------> creating output at:", output$dir))
          # dir.create(output$dir)

          write.table(genemerged, file.path(output$dir, "merged_geneTPMtable.txt"),
                      sep="\t", col.names=NA, quote = F, na="")
 }


# rm(noCounts, noHier)
#save R image
save.image(paste0(output$dir, "/TableMerging_v4.RData") )


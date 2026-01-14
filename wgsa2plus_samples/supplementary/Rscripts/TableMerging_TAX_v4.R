# load packages
library("easypackages") # can also install multiples with packages()
x<-c("tidyverse", "reshape2", "broom", "data.table", "tibble", "docopt", "stats")
suppressMessages(libraries(x))
rm(x)

# library(docopt)
#prep arguments for docopts
doc <- "Usage: my_program.R [options]
--help           show this screen
--binDIR PATH    hard path to bin dir [default: 'getwd()/TEDreadsTAX/bin']
--sccList NAME   hard path and name of file listing successful sample names [default: 'SccList.txt']
--outdir NAME    path to & name of output directory to create [default: 'getwd()/TEDreadsTAX/merged_tables']
--genesDIR PATH  hard path to bin dir for KO/EC genes [None for TAX mode, but flag is needed]"
doc
args <- docopt(doc)
args
# #for manual prep
#  args <-docopt(doc, args = c(
#   "--binDIR" , "/path/to/PWYprofiles/bin/",
#   "--sccList", "/path/to/SccList.txt",
#   "--outdir" , "/path/to/PWYprofiles/merged_tables/" )) #this is if the .sh does not provide arguments
# args


#load("TAX_TableMerging_v1.RData")
# 
# #=============== finished ==============
input=list() #instead of input=NA
input$path= args$binDIR 
input$files=list.files(input$path, full.names = T) #pattern = args$infix,
print("------> detecting input files:")
input$files

input$scclist <- read.table(args$sccList, header = F, sep="\t") %>% 
            setNames("SampNames") %>%  arrange(SampNames) %>% 
            rename_with(~ str_remove(.,".txt"), contains(".txt")) %>% 
            rename_with(~ str_remove(.,"_4krona"), contains("_4krona"))

names(input$files) <- input$scclist$SampNames;
print("------> naming input files:")
input$files

#load data from files
print("------> load data from files")
data=list()
tiers=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species", "TAXid")
data=lapply(input$files, function(x){
  y<-read.table(x, header=F, check.names = F, sep=c("\t"), na.strings =c("NA",""," "),
                fill=T, quote = "", col.names = c("Counts",tiers) ) %>%
    remove_rownames() %>%  group_by(across(all_of(tiers))) %>% mutate(Species=coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom)) %>%
    summarize(across(everything(), sum)) %>%  relocate("Counts") %>% mutate(Species=gsub("_sp._", "_", Species)) %>%
    mutate(Species=gsub("\\[|\\]|Candidatus_", "", Species)) %>%
    mutate(Species=gsub("(.*)", "\\1_sp", Species)) %>% mutate(Species=gsub("(.*_.*)_sp", "\\1", Species)) %>%
    unite("Lineage", all_of(tiers), sep="; ", remove=F) %>%
    relocate("Lineage", .after=last_col()) %>%
    return(y)
})
print("------> checking data loading objects:")
head(data[[1]]); tail(data[[1]])
names(data)


print("------> merging data")
merged=Reduce(function(x,y)
  merge(x,y, by=c(tiers, "Lineage"), all=T ), data)

colnames(merged)=c(tiers, "Lineage", names(data))
# head(merged)
row.names(merged)=paste0("TAX", row.names(merged))
#row.names(merged)=paste0("id", merged$TAXid)
merged = merged %>% relocate("Lineage", .after=last_col()) %>%
  mutate_if(is.numeric, ~replace(. , is.na(.), 0)) %>%
  select(-TAXid) %>% 
  filter(rowSums(select_if(., is.numeric))!=0); #dim(merged) 
print("------> checking merged object:")
head(merged %>% select(-Lineage)); dim(merged);


output=list()
output$dir<- args$outdir
dir.create(output$dir)
print(paste("------> exporting data at:", output$dir))

write.table(merged, file.path(output$dir, "merged_Counts+TAX+Lineage.txt"),
            sep="\t", col.names=NA, quote = F, na="")

noHier <- setdiff(names(merged), c("Lineage"))
write.table(merged[,noHier], file.path(output$dir, "merged_Counts+TAX.txt"),
            sep="\t", col.names=NA, quote = F, na="")

Lin <- select(merged, "Lineage")
write.table(Lin, file.path(output$dir, "merged_Lineage.txt"),
            sep="\t", col.names=NA, quote = F, na="")

noTAX <- setdiff(names(merged), tiers)
noTAX <- merged[, noTAX] %>% relocate(Lineage, .after=last_col() )
write.table(noTAX, file.path(output$dir, "merged_Counts+Lineage.txt"),
            sep="\t", col.names=NA, quote = F, na="")

write.table(merged[,names(data)], file.path(output$dir, "merged_Counts.txt"),
            sep="\t", col.names=NA, quote = F, na="")

noCounts<-setdiff(names(merged), c(names(data), "Lineage") )
write.table(merged[,noCounts], file.path(output$dir, "merged_TAX.txt") ,
            sep="\t", col.names=NA, quote = F, na="")

# bim <-biomformat::make_biom(data = merged[,names(data)], observation_metadata = merged[,noCounts], sample_metadata = names(data))
# biomformat::write_biom(bim, biom_file = paste0(output$dir, "merged_Counts+TAX_json.biom") )
# print("WARNING: please import biom with read_biom(), not with phyloseq, qiime or ampvis biom-import features")

# test biom file:
# library(biomformat)
# bim <- read_biom(paste0(output$dir, "merged_Counts+TAX_json.biom")) #imports
# bim <- amp_import_biom(paste0(output$dir, "merged_Counts+TAX_json.biom")) #does not import
 # library(rbiom)
 # bim=rbiom::read.biom(paste0(output$dir, "merged_Counts+TAX_json.biom")) # imports
 # write.biom(bim, file = paste0(output$dir, "merged_Counts+TAX_hdf5.biom"), format = "hdf5")
# rm(bim)

# rm(noCounts, noHier)
#save R image
save.image(paste0(output$dir, "/TableMerging_v4.RData") )


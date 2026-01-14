#' some random funcitons of mine

#' for creating a kable table
mytab <- function(df, cap=paste0(""), w = as.character(60), h = as.character(500), arrange_by=""){
  ktab <- df %>%  kable("html", caption=cap) %>% 
    kable_styling("striped", "condensed") %>% #, position="float_left" 
    scroll_box(width = paste0(w, "%"), height= paste0(h, "px") )
  return(ktab)
}

##' ----- convert amp2phy & phy2amp -------------------------------------------------------

amp2phy <- function(amp_obj){
  library(phyloseq)
  OTU = otu_table(as.matrix(amp_obj$abund) , taxa_are_rows = T)
  TAX = tax_table(as.matrix(amp_obj$tax) )
  MET = sample_data(amp_obj$metadata)
  
  PHY<- phyloseq(OTU, TAX, MET)
  return(PHY)
}

phy2amp <- function(phy_obj){
  library(ampvis2)
  abund <- otu_table(phy_obj) %>% as.matrix() %>% as.data.frame()
  tax   <- tax_table(phy_obj) %>% as.data.frame()
  meta  <- sample_data(phy_obj) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("rowNames") ## has to be this way
  
  amp <- amp_load(cbind(abund, tax), meta)
  return(amp)
}


##' ----- select a string from list of sublists ------------------------
##' when you have a list of sub-lists, extract from that list only the sublists
##' which contain the mentioned string? 
##' may be my least useful function

sel.from.list <- function(list, string){ 
  sublist <- lapply(list, \(x) { 
    if(string %in% x){ return(x) }   
  } )
  sublist <- sublist[!sapply(sublist, is.null)]
  return(sublist)
}

##' ----- head & tails ----------

heil <- function(mtx, n=10){
  print(head(mtx, n)); print(tail(mtx, n))
}


##' ---- source My Scripts ---

myScripts <- function(script = "Libraries_n_Colors.R"){
  source(paste0(myRpath, script))
}

##' ----- simple tolerance-around-value function or wrap-around-value
##' I use it to see if test-R2 values are within 10% range of the base-R2 values (r)
##' other ways I can call this is tolr, or wrap

rrange <- function(r, range = 0.1) {
  perc <- 100*range
  maxr <- r + range*r
  minr <- r - range*r
  # info <- paste("max ==", maxr, "\nmin ==", 
                # minr, "\ntol ==", perc, "%")
  # info <- paste0(perc, "% range: ", minr, " to ", maxr)
  info <- paste0("wrapOf ", perc, "% = ", minr, " to ", maxr)
  return(cat(info))
}

##' --- see data
##' when dfs are too big head tail or heil show only 
seeTailEnd <- function(df, n = 10){
  corr_val <- 1
  nrows <- dim(df)[1]
  nrmin <- nrows - n ## + corr_val
  ncols <- dim(df)[2]
  ncmin <- ncols - n + corr_val
  
  tailend <- df[nrmin:nrows, ncmin:ncols]
  return(tailend)
}


##' ---- calc fraction of select taxa  (ecological function fraction) ------------------
##' create a taxVector of representative taxa (Genera) with select ecological role (e.g. SCFA producers 
##' or pathogens). This function calculates the relative abundance of these taxa within each sample
##' Use alpha_plot separately, to plot their abundance between groups (use which_alphas)
##' 


calcTAXfrac <- function(ampObj, taxVector, taxRole="repTAXa", taxLevel = "Genus"){
  
  ampObj$fulltab <- ampObj$abund
  ampObj$fulltab[[taxRole]] <- ifelse(ampObj$tax[[taxLevel]] %in% taxVector, "yes", "no")
  table(ampObj$fulltab[[taxRole]])
  
  #3) calculate and add sums of the taxRole taxa
  ampObj$fulltab <- as.data.table(ampObj$fulltab)[, lapply(.SD, sum), by = taxRole ]
  # ampObj$fulltab[[taxRole]]
  
  #' calculate relAbu based on those sums, with 100% = (taxRole_yes + taxRole_no), for each sample
  ampObj$fulltab[, 2:ncol(ampObj$fulltab)] <- 100 * vegan::decostand(ampObj$fulltab[, -1, with = FALSE],
                                                                     method = "total", MARGIN = 2)
  # ampObj$fulltab %>% View
  
  ##' grab only the taxRole_yes and add to metadata table
  # ampObj$metadata$SampObjleID
  ampObj$metadata[[taxRole]] <- ampObj$fulltab[2,-1] %>% t() %>% as.numeric()
  # ampObj$metadata[[taxRole]] %>% class()
  
  # ampObj$fulltab <- NULL # rm(fulltab)
  # ampObj$metadata %>% View
  
  #' plot meta in boxplot
  # ampObj$metadata$ABX_var %<>% factor(levels=c("NA", "NONE", "Single", "Multi" ))
  # aplot <- alpha_plot(ampObj$metadata,fct1 = fct1, fct2=fct2, which_alphas = taxRole, x.angle = 0, 
  #                     sigplot = sigplot, y.lab = "sqrt(Relative Abundance)",  yScale = yScale , ylim = ylim,
  #                     points = T, colVector = colVector ) #+
  # # geom_signif(comparisons=combn(2, 2, simplify = F), test="wilcox.test", map_signif_level = F,
  # #                       step_increase = 0.05, color="gray60", tip_length = 0.01, vjust = 0.55)
  # all <- list(aplot, ampObj$metadata)
  # names(all) <- c("boxplot", "metadata")
  # return(all)

  return(ampObj$metadata)
}

##' sets:
taxVec <- list()
taxVec$PTGNtax <- c( ## known pathogenic genera (Bacteria; human gut)
              "Acinetobacter", "Aeromonas", "Anaplasma", "Vibrio", "Yersinia", 
              "Bacillus", "Bartonella", "Bordetella", "Brucella", "Burkholderia", 
              "Campylobacter", "Chlamydia", "Clostridium", "Corynebacterium", 
              "Coxiella", "Enterococcus", "Escherichia", "Francisella", "Haemophilus", 
              "Helicobacter", "Klebsiella", "Legionella", "Listeria", "Mycobacterium", 
              "Mycoplasma", "Neisseria", "Pseudomonas", "Rickettsia", "Salmonella", 
              "Shigella", "Staphylococcus", "Streptococcus")
taxVec$TLOCtax <- c( ## putative translocating genera (human gut)
              "Segatella", "Butyricimonas", "Bordetella", "Acinetobacter", 
              "Citrobacter", "Escherichia", "Klebsiella", "Shigella", "Salmonella", 
              "Collinsella", "Actinomyces", "Mycolicibacterium", "Corynebacterium", 
              "Veillonella", "Lachnospira", "Ihubacter", "Gemella", "Kurthia", 
              "Enterococcus", "Abiotrophia", "Aerococcus", "Lactococcus", "Streptococcus", 
              "Weissella", "Limosilactobacillus", "Lactobacillus", "Ligilactobacillus" )

## curated last Apr 2025
taxVec$BUTYtax <- c( ## butyrate producers (SCFA producers) 
  "Faecalibacterium", "Roseburia", "Butyricicoccus", "Eubacterium", 
  "Agathobacter", "Ruminococcus", "Coprococcus", "Butyrivibrio", 
  "Anaerostipes", "Oscillospira", "Blautia", "Lachnospira", 
  "Roseburia", "Anaerotruncus", "Paraclostridium", "Poryphyromonas",
  "Subdoligranulum", "Megasphaera", "Tannerella", "Instestinimonas")
taxVec$PROPtax <- c( ## propionate producers (SCFA producers)
  "Bacteroides", "Prevotella", "Dialister", "Coprococcus", "Veillonella", 
  "Blautia", "Megasphaera", "Roseburia", "Acidaminococcus", "Ruminococcus", 
  "Akkermansia", "Intestimonas", "Peptoniphilus", "Lactobacillus", 
  "Enterococcus", "Listeria", "Fusobacterium", "Phascolarctobacterium" )
taxVec$ACETtax <- c( ## acetate producers (SCFA producers)
  "Bacteroides", "Prevotella", "Bifidobacterium", "Escherichia", 
  "Akkermansia", "Ruminococcus", "Faecalibacterium", "Roseburia", 
  "Blautia", "Coprococcus", "Lachnospira", "Dorea", "Anaerostipes", 
  "Clostridium", "Eubacterium", "Butyrivibrio", "Veillonella", 
  "Megasphaera", "Collinsella", "Desulfovibrio",
  "Anaerotruncus", "Oscillospira", "Parabacteroides", "Subdoligranulum" )
taxVec$SCFAtax <- unique( c(taxVec$BUTYtax, taxVec$PROPtax,taxVec$ACETtax) )
##--

taxVec$PTGNeuk <- c( ## Euk genera pathogenic in **human & mouse** (elected from VEukPathDB; 
                     #https://tinyl.io/BFyf from Release 68,  May 7, 2024, with help of chatGPT)
              "Acanthamoeba", "Alternaria", "Apophysomyces", "Aspergillus", "Babesia", 
              "Balamuthia", "Blastomyces", "Candida", "Coccidioides", "Cryptococcus", 
              "Cryptosporidium", "Cyclospora", "Cystoisospora", "Dermacentor", "Emergomyces",
              "Encephalitozoon", "Entamoeba", "Exophiala", "Fonsecaea", "Fusarium", "Giardia",
              "Haemaphysalis", "Histoplasma", "Ixodes", "Leishmania", "Lichtheimia", "Madurella", 
              "Malassezia", "Microsporidium", "Microsporum", "Naegleria", "Neospora", "Nosema",
              "Paracoccidioides", "Pediculus", "Penicillium", "Phialophora", "Phlebotomus", 
              "Plasmodium", "Pneumocystis", "Rhipicephalus", "Rhizopus", "Sarcocystis", 
              "Scedosporium", "Sporothrix", "Stachybotrys", "Talaromyces", "Theileria", 
              "Toxoplasma", "Trichomonas", "Trichophyton", "Trichosporon", "Trypanosoma", 
              "Uncinocarpus") ## there is 225, but those are primary


pwyVec <- list() ## mCYC pathways related to SCFA production
pwyVec$BUTYpwy <- c("PWY-5677", "P163-PWY", "PWY-5676","CENTFERM-PWY") ## Fermentation to Butanoate, 
## butanoate degradation pathways are "PWY-5022","GLUDEG-II-PWY" (degradative); CENTFERM-PWY is both
pwyVec$PROPpwy <- c("PWY-5088", "PWY-7013", "P108-PWY", "PWY-5494") ## Fermentation to Propanoate 
## propanoate degradative  "PWY-8188",PROPFERM-PWY, PWY-8086
pwyVec$ACETpwy <- c( "PWY-5100", "PWY-8328", "P142-PWY", "P461-PWY", "P162-PWY", "P41-PWY",  ## Fermenation to Acetate
                     "P124-PWY", "PWY-804", "PWY-5535", "PWY-5600", "PWY-5536", "PWY-5537",  "P122-PWY",
                     "PWY-5096", "PWY0-1312", "P161-PWY", "PWY-5768", "ANAEROFRUCAT-PWY", "PWY-8190") 
## acetate degradive pwhys are: "PWY-5482", "PWY-5485", "PWY-5483", "P163-PWY", "PROPFERM-PWY", "PWY-5538", "P142-PWY", "PWY-8377"  
pwyVec$SCFApwy <- c("PWY-5676", "ANAEROFRUCAT-PWY", "P124-PWY", "PWY-7013",  ## "PWY-5022",#PROPFERM-PWY, "PWY-8188",  degradative
                    "P461-PWY", "P161-PWY", "P122-PWY", "CENTFERM-PWY",   
                    "GLUDEG-II-PWY", "PWY-8190", "P162-PWY", "P163-PWY", "PWY-5494", 
                    "PWY0-1312", "PWY-5535", "PWY-8328", "PWY-5536", "PWY-8377", 
                    "PWY-5677", "PWY-5096", "P142-PWY", "PWY-5537", "PWY-5600", 
                    "P41-PWY", "PWY-5100", "P108-PWY") ## n = 26, checked 


##' ----- get dominant (Katie) -----
##' this is katie's function, assuming x is a phy object
##' for how to use it go to 
##' ~/Documents/CoreMicroboime/SuchiHourgan_infants/2024projects/111NICU_project/reports/PreliminaryDataSummary.html

get_dominant <- function(phyObj, level="Species") {
  # if(level == "Species") {
  #   name <- c(tax_table(x)[apply(otu_table(x), 2, which.max), "Species"])
  #   val <- c(apply(otu_table(x), 2, max))
  # } else {
    glommed <- tax_glom(phyObj, level)
    name <- c(tax_table(glommed)[apply(otu_table(glommed), 2, which.max), level])
    val <- c(apply(otu_table(glommed), 2, max))
  # }
  ling <- tax_table(glommed) %>% as.data.frame() %>% filter(!!sym(level) %in% as.vector(name))
  return(list(name=name, val=val, tax=ling))
}






##' ------- merge OTU tables based on Lineage ------------------------
##' This should be the more complex version, like in AmFootballers 1-Metadata processing_v1.Rmd (chunk ln 150 or so)
##' 
##' 


##' -------- MERGE tables -------------------

## function for merging multi-sample otuMTXs (like merged_tables/.. from wgsa2)
## ways to use:
## if:
# > files[2:4]
# $otuTAB1
# [1] "/Users/angelovaag/Documents/Collaborations/LiscoMGX/wgsa2out_web/pt1/merged_tables/merged_Counts+TAX.txt"
# $otuTAB2
# [1] "/Users/angelovaag/Documents/Collaborations/LiscoMGX/wgsa2out_web/pt2/merged_tables/merged_Counts+TAX.txt"
# $otuTAB3
# [1] "/Users/angelovaag/Documents/Collaborations/LiscoMGX/wgsa2out_web/pt3/merged_tables/merged_Counts+TAX.txt"
## then:
# mergeTables(files[2:4])


mergeTables <- function(FileNamesList){
  
  data<- list()
  data=lapply(FileNamesList, function(x){
    y<-read.table(x, header=T, sep="\t", row.names = 1, 
                  check.names=F, quote="" ) %>% mutate(rowname=rownames(.))
    return(y)
  })
  
  
  
  taxRanks <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  otuMTX <- Reduce(function(x,y) merge(x,y, by=c("rowname", taxRanks), all=T), data) %>%
    group_by(rowname, Kingdom, Phylum, Class, Order, Family, Genus,Species) %>% 
    summarize(across(where(is.numeric) ,sum, na.rm=T)) %>% 
    column_to_rownames("rowname")
  
  
  return(otuMTX)
}





##' -------- EXTRACT object from RData -------------------

myRDataExtractor <- function(RDataFile, objectName) {
  #' Function for extracting an object from a .RData file created by R's save() command
  #' Inputs: RData file, object name
  E <- new.env()
  load(RDataFile, envir=E)
  return(get(objectName, envir=E, inherits=F))
}
# Rmkdown styles

# ---
# title: "IgAseq.PRE_RAG_diversity-all_v4"
# author: 
#   - name: "Angelina Angelova, PhD"
# affiliation: "BCBB/OCICB/OSMO/OD/NIAID/NIH"
#   - name: "[Metagenomics @NIAID](https://tinyl.io/97Be)"
# date: "`r Sys.Date()`"
# output: 
#   html_document: 
#   css: ~/Documents/MyApps/RScripts/myRscripts/CSSstyles.css
# code_folding: hide
# highlight: pygments
# df_print: kable
# theme: lumen
# editor_options: 
#   chunk_output_type: console
# markdown: 
#   wrap: 180
# ---


# Other font families I like: Arial Narrow, Optima, Palatino, Gill Sans , Comic Sans, Courier, Bradley Han
# Font families at: https://www.w3.org/Style/Examples/007/fonts.en.html 
# h1 {font-size: 34px;}
# h1.title {font-size: 38px;}
# h2 {font-size: 30px;}
# h3 {font-size: 24px;}
# h4 {font-size: 18px;}
# h5 {font-size: 16px;}
# h6 {font-size: 14px;}


#Load common libraries
library(easypackages) #can also install multiples with packages()
x<-c("tidyverse", "reshape2", "broom", "dplyr", "data.table", "readr",         # data organization
     "stringr", "stringi",  "tibble", "magrittr",                              # data organization
     "ggplot2","ggrepel", "kableExtra", "ggpubr", "plotly",                    # plotting & Visuals
     "devtools", "htmltools", "readxl", "parallel", "BiocParallel",            # dev tools
     "scales","RColorBrewer", "colorspace")                                    # colors                                        #stats
suppressMessages({libraries(x)}); rm(x)
# suppressMessages( lapply(x, library, character.only=T) ) ## too much outprint

#' loading microbiome & stats packages; info: https://tinyurl.com/y3syd5vn
y<-c("ape", "Biostrings", "fossil", "vegan", "car", "lmerTest",
     "ampvis2", "phyloseq", "mixOmics", "MetaLonDA", ## "microbiomeSeq", 
     "MMUPHin", "DESeq2", "ALDEx2", "Maaslin2", "EnhancedVolcano", 
     "ggsignif", "rstatix", "nestedRanksTest","pairwiseAdonis")
libraries(y); rm(y)

# Optional/supportive packages
# library(kableExtra) # for kable Tables
# library(tsne)    # different type of ordination
# library(plotly)  # plots ly?

#' path for myR functions
myRpath <- paste0("~/Documents/MyApps/RScripts/myRscripts/")
source(paste0(myRpath, "dataCorrections.R"))
source(paste0(myRpath, "dataClean.R"))
source(paste0(myRpath, "otherFunctions.R"))
source(paste0(myRpath, "Diversity.R"))
source(paste0(myRpath, "effectEval.R"))
source(paste0(myRpath, "myDiffAbuFunctions.R"))
# source(paste0(myRpath, "myDESeq.R"))
# source(paste0(myRpath, "myALDEX.R"))
# source(paste0(myRpath, "myMaAsLin.R"))



#' options: https://web.mit.edu/r/current/lib/R/library/base/html/options.html
#'  names(options()); getOptions("width")  or options()$width
#' getOption("max.print"); 
options(max.print = 200) # [def: 1000] # when print(mtx), how many max lines to print
options(warn = 1) ;# [def: 0 ] sometimes warnings are converted to errors by R. This sets a reduced threshold for warning-to-error conversion
#' If warn == 0 (default), warnings are stored until the top–level function returns. Signaled warnings 10+, R will only report their number, not actual message
#' If warn == 1, warnings are printed as they occur. 
#' If warn >= 2, all warnings are turned into errors.
options(width = 155) # [def: 80 ] the max # of characters printed on a line
options(digits = 5)  # the 0.+decimals displayed for numbers.
options(scipen = -1)  # [def: 0] penalty level for conversion sci number notation. Neg values: prefers sci notation. Pos values: prefers fixed notation.
#' fixed notation is preferred unless it is more than scipen digits. Check with:
# print(list(1e-7,1e-6,1e-5,1e-5,1e-4,1e-3,1e-2,1e-1,0,1,1e1,1e2,1e3,1e4,1e5,1e6,1e7))
Sys.setenv(PATH='/usr/bin:/bin:/usr/sbin:/sbin:/usr/X11/bin/:/usr/local/bin/')

#' multithreading
control.compute=list(save.memory=TRUE) ## for saving memory when only 16G is available
# library(parallel); library(BiocParallel)
threads=detectCores()-4
BiocParallel::register(MulticoreParam(threads))
# getOption("mc.cores") ## == NULL
options(mc.cores=threads)
setDTthreads(threads)
# getOption("mc.cores") ## == threads
#' ------------------------- Create pretty colors palletes & gradients
  #' Create a distinctive color pallete (a color_vector)
  # library(RColorBrewer)
  # display.brewer.all()
  #   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #   col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  #   n=15 #pick between 2 and 70
  #   set.seed(567)
  #   tol=sample(col_vector, n);   tol
  #   pie(rep(1,n), col=tol, labels = tol, cex=0.8)
  #   # col0=c("deepskyblue1", "palevioletred1", "seagreen1", "khaki", "maroon", "tomato", "plum1",  "chartreuse1", "#8DD3C7", "#FFF2AE" , "snow1", "chocolate2", tol)
  # pie(rep(1, length(col)), col=col, labels = col, cex=0.8)
  # # you can replace similar colors like this
  # tol=replace(col, col %in% c("#FFFF33", "#FFF2AE", "#B15928"), c("snow1", "deeppink3", "peachpuff"))
  # long_colors <-c(brewer.pal(12, "Set3"), brewer.pal(8, "Accent"))
#' === Create discrete colors vector of chosen colors
# colDscrt <- c("#F8766D", "#00BFC4", "#C77CFF", "#7CAE00", "#CD9600", "#00A9FF",  "#F0027F","seagreen1", "darkkhaki", 
#       "maroon","#CCEBC5", "tomato", "plum1",  "chartreuse1", "#8DD3C7", "#FFF2AE" , "chocolate2", "#999999",
#       "#FDB462", "snow1","#F0E442", "palevioletred1", "#6A3D9A", "#CBD5E8", "#FFD92F", "#D9D9D9", 
#       "deepskyblue1","#FFFF99","#80B1D3", "#B15928","#E31A1C", "#BC80BD", "#CCEBC5", "#FFED6F", "#666666",
#       "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",  "#7FC97F", "#BEAED4", "#386CB0",
#       "#FB8072","#BEBADA", "#FFFFB3", "#7FC97F", "#FBB4AE")
# library(scales)
# show_col(col)
# length(col)
#' === Create scaled colors 
# display.brewer.all() #https://ggplot2-book.org/scale-colour.html
# colScales <- c()
# # display.brewer.pal(5, "Reds");
# colScales$red <- brewer.pal(5 , "Reds");
# # display.brewer.pal(5, "Blues");
# colScales$blue <- brewer.pal(n=5, name = "Blues")  
# # display.brewer.pal(5, "Oranges");
# colScales$oranges <- brewer.pal(n=5, name = "Oranges")  
# # display.brewer.pal(5, "Greens");
# colScales$greens <- brewer.pal(n=5, name = "Greens")
# # display.brewer.pal(5, "Purples");
# colScales$purples <- brewer.pal(n=5, name = "Purples")
# colScales$pastels <- colorRampPalette(brewer.pal(12, "Set3"))(20);

myColors <- list()
  #' Create discrete color pallet
myColors$discrete <-  
  c("#F8766D", "#00BFC4", "#C77CFF", "#7CAE00", "#CD9600", "#00A9FF","#666666", "#F0027F","seagreen1",
    "maroon",  "#CCEBC5", "tomato", "plum1",  "chartreuse1", "#8DD3C7", "#FFF2AE" , "chocolate2", "#999999",
    "snow1",   "#F0E442", "palevioletred1", "#6A3D9A", "#CBD5E8", "#FFD92F", "#D9D9D9", "#FBB4AE",  "darkkhaki",
    "deepskyblue1","#FFFF99","#80B1D3", "#B15928","#E31A1C", "#BC80BD", "#CCEBC5", "#FFED6F", "#FCCDE5", 
    "#FDB462", "#B3DE69", "#7FC97F", "#BEAED4", "#386CB0",
    "#FB8072") #"#FFFFB3",  "#FDB462",
# names(myColors$discrete) <- 1:length(myColors$discrete)

# Generate custom pallets based on Hexadecimal codes: https://mycolor.space/
myColors$dscr8 <-  c("#0072B2","#D55E00","#CC79A7","#E69F00","#56B4E9","#009E73", "#F0E442","#999999")
myColors$dscrX <-  c("#00BFC4","#00ABD7","#2093DB","#7374C8","#9D509D", "#AA2C62","#9aad53",  "#f8766d", "#d39545", "#009E73")
myColors$salmon <- c("#FFC9BC", "#F8766D", "#C660A0", "#9762AE", "#5D64AF", "#0062A3","#9C2427")
myColors$pastel <- c("#8DD3C7", "#CFECBB", "#F4F4B9", "#CFCCCF", "#D1A7B9", "#F4867C",
                     "#C0979F", "#86B1CD", "#CEB28B", "#EDBC63", "#C2D567", "#CDD796",
                     "#F8CDDE", "#E9D3DE", "#C59CC5", "#C9DAC3", "#E1EBA0", "#FFED6F")

#' color name to hexadecimal code (function)
col2hex <- function(x, alpha = FALSE) {
  args <- as.data.frame(t(col2rgb(x, alpha = alpha)))
  args <- c(args, list(names = x, maxColorValue = 255))
  do.call(rgb, args)
}
# col2hex(colors(), alpha = TRUE)
# col2hex("orangered3")

    library(scales)
    show_col(myColors$discrete)

#' Create color gradients
myColors$reds <- brewer.pal(5 , "Reds");
# display.brewer.pal(5, "Blues");
myColors$blues <- brewer.pal(n=5, name = "Blues")  
# display.brewer.pal(5, "Oranges");
myColors$oranges <- brewer.pal(n=5, name = "Oranges")  
# display.brewer.pal(5, "Greens");
myColors$greens <- brewer.pal(n=5, name = "Greens")
# display.brewer.pal(5, "Purples");
myColors$purples <- brewer.pal(n=5, name = "Purples")
# display.brewer.pal(5, "Greys")
myColors$greys <- brewer.pal(n=5, name="Greys")
# display.brewer.pal(7, "Spectral")
myColors$spectral <- brewer.pal(n=7, name="Spectral")
# display.brewer.pal(7, "Greys")


# ---------- Create pch vector that I like (up to 19 values):
myPchR <- list()
myPchR$all <- c(16,15,17,18,   1,0,5,6,   3,4,8,   10,12,7,9,11,13,14,  21:25,  19,20)
myPchR$main <- c(16,15,17,18)
myPchR$fill <- c(21:25)
myPchR$empt <- c(1,0,2,6,5)
myPchR$crss <- c(3,4,8)
myPchR$strs <- c(10,12,7,9,11,13,14)
# ggpubr::show_point_shapes()



# ---------- setting knitr global options
knitr::opts_chunk$set(warning = FALSE, message = FALSE,  #out.width = '\\textwidth',
                      tidy.opts = list(width.cutoff=220), tidy = T, echo = TRUE,
                      fig.width = 8, fig.height = 6, width = 1200)#000)#, size="small")#
# df %>% kable("html") %>% kable_styling("striped") %>% scroll_box(width = "100%")
#' results defulat is "markup", but also try 'hold' (to output result at end of chunk exe) & "verbatim" to print exactly as in console

# common variables
              # c("***"= 0.001, "**"=0.01, "*"= 0.05, "." = 0.065, " "=1)
sig.cuts  <- c(0,      1e-04, 0.001, 0.01,  0.05, 0.8, 1)
sig.stars <- c(       '****', '***', '**' , '*' , '.', ' ')
taxRanks <- c("Kingdom", "Phylum", "Class", "Order","Family", "Genus","Species")

#' set 6threads to data.tables
# setDTthreads(6) #does not work
# getDTthreads(verbose=T) #if it still says 1, its problems for another day


## ----- load a ZymoTheoretical profile
ZymoTheo <- readRDS("~/Documents/myApps/RScripts/myRscripts/ZymoEvenTheoretical_amp.RDS")


## make a standard Zymo Mock amp
# ZymoSpecies <- c("Pseudomonas aeruginosa", "Escherichia coli", "Salmonella enterica", "Limosilactobacillus fermentum", 
#                  "Enterococcus faecalis", "Staphylococcus aureus", "Listeria monocytogenes", "Bacillus subtilis",  
#                  "Saccharomyces cerevisiae",  "Cryptococcus neoformans" )
# ZymoLineags <- abutax$tax %>% filter(Species %in% ZymoSpecies) %>% select(-SPP) %>% remove_rownames() %>%
#   mutate(OTU  = paste0("spp", 1:nrow(.)))  %>% mutate(rn = OTU) %>% column_to_rownames("rn")
# ZymoLineags
# ZymoAbunds <-  matrix(nrow = 10, ncol = 2) %>% as.data.frame() %>% setNames(c("SPP", "ZymoEven")) %>%
#   mutate(SPP = ZymoLineags$OTU) %>% mutate(ZymoEven = c(rep(12, 8), 2, 2) / 100) %>% column_to_rownames("SPP")
# ZymoAbunds
# 
# ZymoMeta <- matrix(nrow = 1, ncol = 2) %>% as.data.frame() %>% setNames(c("SampleID", "SampleType")) %>%
#   mutate(SampleID = "ZymoEven", SampleType = "ZymoTheoretical")
# 
# ZymoTheo <- amp_load(otutable = ZymoAbunds, taxonomy = ZymoLineags, metadata = ZymoMeta) %>% amp_agg_abund()
# 
# heat_plot(ZymoTheo, facet_fct = "SampleType", group_fct = "SampleID", title = NULL,subtitle = NULL)
# 
# saveRDS(ZymoTheo, "~/Documents/myApps/RScripts/myRscripts/ZymoEvenTheoretical_amp.RDS")

# To be used as:
# ampZymoMocks %<>% amp_agg_abund("Species")
# ampZymos0 <- amp_merge_ampvis2(ampZymoMocks, ZymoTheo, by_refseq = F) 
# ampZymos0$metadata %<>% mutate(SampleType = factor(SampleType, levels = c("ZymoMock", "ZymoPCR", "Theoretical")))
# heat_plot(ampZymos0,  facet_fct = "SampleType", group_fct = "SampleID", normalize = T, ntax = 15, x.angle = 45)
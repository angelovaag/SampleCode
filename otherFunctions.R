#' some random funcitons of mine

#' for creating a kable table
mytab <- function(df, cap=paste0(""), wdth = "60", hth = "500", arrange_by=""){
  ktab <- df %>%  kable("html", caption=cap) %>% 
    kable_styling("striped", "condensed") %>% 
    scroll_box(width = paste0(wdth, "%"), height= paste0(hth, "px") )
  return(ktab)
}

##' ----------------------- convert amp2phy & phy2amp -------------------------------------------------------

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
  abund <- as.matrix(otu_table(phy_obj))
  tax   <- tax_table(phy_obj) %>% as.data.frame()
  meta  <- sample_data(phy_obj) %>% as.data.frame()
  
  amp <- amp_load(cbind(abund, tax), metadata)
  return(amp)
}

##'/--------------------------------------------------------------------------

##' ----------- select a string from list of sublists ------------------------
##' when you have a list of sub-lists, extract from that list only the sublists
##' which contain the selected string? 
##' may be my least useful function

sel.from.list <- function(list, string){ 
  sublist <- lapply(list, \(x) { 
    if(string %in% x){ return(x) }   
  } )
  sublist <- sublist[!sapply(sublist, is.null)]
  return(sublist)
}
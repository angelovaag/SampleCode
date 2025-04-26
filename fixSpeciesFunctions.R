###' --------- FIX species ------------
#' functions fix empty fields of TAX$Species table to be represented as:
#' LowestAvailableTAXrank_speciesName
#' 
#' functions assume
#' - tax table is separated per rank already (each rank has in its own column)
#' - headers of ranks are "Kingdom,Phylum,...,Genus,Species
#' - symbol for missing Species names is provided as the object "mpt" ( e.g. "", "<NA>", "none"). Use no quotes for "NA"
#' 


#' One FixSpecies to rule them all! combine fixSpp1 & 2  
#' troubleshooting done, might still need more custom, but should generally be fine

fixSpecies<- function(taxMtx, mpt, addGenus=T){
  taxdf <- as.data.frame(taxMtx)
  # nr=c(70:90);  nc=c(4:7)
  # taxdf[nr, nc]
  
  
  taxdf %<>% mutate_all(~str_replace_all(., "Incertae Sedis|group|Candidatus|.__", "")) %>% 
    dplyr::mutate_all(~na_if(., mpt))
  # taxdf[nr, nc]
  
  if(addGenus==T){   
    #' if Genus needs adding to Species names (e.g output from dada2/SILVA db)
    taxdf %<>% mutate(Species=ifelse(is.na(Species), paste0(coalesce(Genus,Family,Order,Class,Phylum,Kingdom), " sp.") , paste(Genus, Species) ))
    # taxdf[nr, nc]
  }else{    
    #' if Species names already contain Genus string (pattern "Genus Species"), then 
    #' suffix & coalescence  is all that is needed to the empty Spp names
    #' (e.g. output form GTdb)
    taxdf %<>% mutate(Species=ifelse(is.na(Species), paste0(coalesce(Genus,Family,Order,Class,Phylum,Kingdom), " sp.") , Species ))
    # taxdf[nr, nc]
  }
  #' Rm placeholder strings
  taxdf %<>% mutate(Species=gsub(" group|Candidatus ", "",    Species)) %>% 
              mutate(Species=gsub("_sp"              , " sp.", Species)) %>%
              mutate(Species=gsub("_"                , " "   , Species))
  
  return(taxdf)
}


###' --------- BATCH correction ------------

batchCorr <- function(amp, corr="Batch", covars=NULL){
  library(MMUPHin)
  
        if(!is.null(covars) ){  
        amp$metadata[, covars] <- lapply(amp$metadata[, covars], as.factor) #multiple factors
        # amp$metadata[[covars]] %<>% factor()
        }
   amp$metadata[[corr]] %<>% factor() # 1 factors
   set.seed(124)
   bch_adj <- adjust_batch(feature_abd = amp$abund,
                        data= amp$metadata,
                        batch= corr,
                        covariates=covars)
  
  ampCORR <- amp_load(cbind(bch_adj$feature_abd_adj, amp$tax), amp$metadata )
  # return(bch_adj$feature_abd_adj)
  return(ampCORR)
}

###' ---------- Batch correction with sva::ComBat() 
#' usually performs worse than mmuphin 
#' #' library(sva)
#'
#' 1) input is abundance matrix 
#'    ampCLN$abund[1:5, 1:10] 
#' 2) needs normalization of the mtx
#'    ampNRM <- amp_rarefy(ampCLN, min(colSums(ampCLN$abund)))
#' 3) corrects only for 1 effect, use non-parametric (par.prior=F)
#'    corrCB <- ComBat(ampNRM$abund, batch = ampCLN$metadata$Batch, mean.only = T ,par.prior = F) %>% round(0)
#'    corrCB[1:5, 1:10]
#' 4) needs setting negative scores to 0
#'    corrCB[corrCB <0 ] <- 0  # some neg values will appear, need set to 0

batchComBat<- function(amp, corr="Batch", rareAT=min(colSums(amp$abund)) ){
  library(sva)
  ampRARE <- amp_rarefy(amp, rareAT)
  corrCB<- ComBat(ampRARE$abund, batch=amp$metadata[[corr]],
                  mean.only = T, par.prior = F) %>% round(0)
  corrCB[corrCB<0]<- 0
  
  ampCB <- amp_load(cbind(corrCB, amp$tax), amp$metadata )
  # ampCB <- amp
  # ampCB$abund <- corrCB
  return(ampCB)
}


###' ------------- ConQuR batch correction
###' https://github.com/wdl2459/ConQuR
###' https://github.com/wdl2459/ConQuR/blob/main/ConQuR_2.0.pdf





###' --------- BIAS correction pre DA analyses-------------------
###' info at: 
###; https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
###' and more at: https://github.com/FrederickHuangLin/ANCOMBC 
##' fix formula = fct$all
##' ANCOM-BC2 estimates and corrects both: 1) the sample-specific (sampling fraction), and 
##' 2) the taxon-specific (sequencing efficiency) biases.
##' It is important for pre- DA analyses, especially with ANCOM-BC

biasCorr <- function(amp, fix_factors=NULL , group= NULL , taxlvl=NULL , p.adj.met = "holm",
                     rand_formula = NULL, asvPREV=0.01){ 
  # taxlvl="Species" <- will agglomerate at spp for correction. 
        # txlvl should be NULL for no agglom to be performed performed (better but slow)
  # group = If the group of interest contains only two categories, leave it as NULL 
        # group should be only 1 factor with >2 categories
  # rand_formula should be  "(1|SubjectID)" or "(TimePoint|SubjID)" for DA, or NULL for BC
  #  prv_cut = 0.1 is default. that is 10% asv PREVALENCE. That is high!
  
  fix_formula <- paste(fix_factors, collapse=" + ")
  
  if(!is.null(fix_factors)){
  amp$metadata[, fix_factors] <- lapply(amp$metadata[, fix_factors], as.factor)
  # amp$metadata[[fix_factors]] %<>% factor() #does not work with multiple factors
  }

  library(phyloseq)
  phy=phyloseq(otu_table(amp$abund, taxa_are_rows = T), 
                    sample_data(amp$metadata), amp$tax)
  
  library(parallel)
  n_cl <- (detectCores()-4 ) ;  # setting for  12 cores. ( !>16!!! for the 8 core L444)
  # cl <- makeCluster(n_cl, type="FORK") # # Mac/Linux need to set as "FORK". Else "PSOCK"
  
  suppressWarnings(library(ANCOMBC))
  set.seed(124)
  anbc_out = ancombc2(data = phy, group=group, tax_level = taxlvl,
                      fix_formula = fix_formula, rand_formula=rand_formula,
                      p_adj_method = p.adj.met, prv_cut = asvPREV,
                      pseudo_sens = FALSE, verbose=T, n_cl= n_cl) 
  closeAllConnections() ## clean up parallel connections before starting
  # stopCluster(cl)
  # closeAllConnections()
  #' pseudo_sense is T by default but it takes a while, and produces no help
  
  ## we get a fraction with which to correct each sample's log(counts): anbc_out$samp_frac
  
  ## we have to replace NA in the frac with 0 
  samp_frac=anbc_out$samp_frac %>% replace(is.na(.), 0)
  
  # Add pesudo-count (1) to avoid taking the log of 0. And take log of counts
  #  more info at: https://github.com/FrederickHuangLin/ANCOM-BC-Code-Archive/blob/master/scripts/figure_5.Rmd#L105
  log_counts = log(anbc_out$feature_table + 1) # CODA adjustment
  
  # Adjust the log observed abundances
  log_biascor = t(t(log_counts) - samp_frac) # adjust log_counts with the samp_frac
  
  # Convert back to counts (correcting for the pseudo-count)
  biascor = round(exp(log_biascor)-1, 0)
  
  # Correct <0 coounts, which should be 0
  biascor[biascor <0] <- 0 
  
  # Prep returns
  outs <- list(anbc_out, biascor)
  names(outs) <- c("ancomBC_out", "corrMTX")
  return(outs)
}

#####' -------- MERGE tables -------------------

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

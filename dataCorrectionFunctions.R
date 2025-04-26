##' -------- normalize compositional data -----------------

normMyData <- function(countsMTX, norm_method = "RelAbu", pseudocount = 1){
  
  ##' --------  calculate log2(MTX) from amp_obj
  counts <- countsMTX
  tcounts <- t(countsMTX)
  
  # pseudocount <- min(counts[counts>0], na.rm = T)
  
  if (norm_method == "LogRel" ) {                       ## standard log2(relAbu)
    ylabel = "Log2(Relative Abundance)";  yinfo = paste("normalized by", ylabel); print(yinfo)
      relabu <- tcounts / rowSums(tcounts) 
    txfCounts <- log2( relabu + pseudocount )
    txfCounts %<>% t() %>% as.data.frame() 
    
  } else if (norm_method=="LogAbu") {                  ### log of raw abundance
    ylabel = "Log2(abundance)"; yinfo = paste("normalized by", ylabel); print(yinfo)
    txfCounts <- log2( (tcounts + pseudocount) )
    txfCounts %<>% t() %>% as.data.frame() 
    
  }else if (norm_method == "Umers"){   ### Umer's missinterpreted math == not correct
    #' the log2(abd) plot: log-ratios in RelAbu are captured
    ylabel = "compensated pseudocount abundance";  yinfo = paste("normalized by", ylabel); print(yinfo)
    txfCounts <- log((tcounts + pseudocount)/(rowSums(counts) + dim(tcounts)[2]))  ## dim() compensates for the 1s added: ?BS
    txfCounts %<>% t() %>% as.data.frame()
    
  }else if (norm_method == "TSS") {                      ## Total sum scaling (norm to 1, independent for taxa)
    ylabel = "Proportional abundance (TSS)";  yinfo = paste("normalized by", ylabel); print(yinfo)
    txfCounts <- tcounts / rowSums(tcounts)
    txfCounts %<>% t() %>% as.data.frame()
    
  }else if (norm_method == "CoDA") { ## CoDA transf (like TSS but taxa are interdependent; a pre-CLR state, really)
    ylabel = "Compositional abundance (CoDA)";  yinfo = paste("normalized by", ylabel); print(yinfo)
    txfCounts <- compositions::acomp(tcounts)
    txfCounts %<>% t() %>% as.data.frame() 
    
  }else if (norm_method == "CLR") { ## centered log ratio (uses the CoDA method internally)
    ylabel = "Centered Log Ratio abundance (CLR)";  yinfo = paste("normalized by", ylabel); print(yinfo)
    txfCounts <- compositions::clr(tcounts + pseudocount) ## it will auto-add 1, if counts has 0s
    txfCounts %<>% t() %>% as.data.frame()
    
  }else if (norm_method=="RelAbu") { ## standard RelAbu, to 100%
    # the abd plot: the RelAbu is captured
    ylabel = "Relative abundance (%)"; yinfo = paste("normalized by", ylabel); print(yinfo)
    # txfCounts <- 100 * vegan::decostand(tcounts, method = "total", MARGIN =1) # mrgn=2 if tcounts used
    txfCounts <- 100 * vegan::decostand(counts, method = "total", MARGIN = 2) 
    txfCounts %<>% as.data.frame()
    
  }else if (norm_method=="none"){
    ylabel = "Abundance"; yinfo = paste("normalized by", ylabel) ; print(yinfo) ;
    txfCounts <- counts
    txfCounts %<>% as.data.frame()
  }else{
    print("choices of norm_method are: 'none', 'RelAbu', 'LogRel', 'LogAbu', 'TSS', 'CoDA', 'CLR' or 'Umers' ")
      }
  
  outObj <- list(txfCounts, ylabel); names(outObj) <- c("normCounts", "ylabel") 
    return(outObj)
  
}


###' --------- BATCH correction ------------

batchMMUPHin <- function(amp, corr="Batch", covars=NULL){
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

##' ----------  sva::ComBat() Batch correction 
###' usually performs worse than mmuphin 
###' best of RNAseq (expression) data
#' #' library(sva)
#'
#' 1) input is abundance matrix The input data are assumed to be cleaned and normalized before batch effect removal.
#'    ampCLN$abund[1:5, 1:10] , expected to be cleaned and normalized (rarif ok) as input
#' 2) needs normalization of the mtx
#'    ampNRM <- amp_rarefy(ampCLN, min(colSums(ampCLN$abund)))
#' 3) corrects only for 1 effect, use non-parametric (par.prior=F)
#'    corrCB <- ComBat(ampNRM$abund, batch = ampCLN$metadata$Batch, mean.only = T ,par.prior = F) %>% round(0)
#'    corrCB[1:5, 1:10]
#' 4) needs setting negative scores to 0
#'    corrCB[corrCB <0 ] <- 0  # some neg values will appear, need set to 0

batchComBat<- function(amp, corr="Batch", rareAT=min(colSums(amp$abund)) ){
  # BiocManager::install("sva")
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



###' ------------- ConQuR batch correction (cannot get it to work )
###' https://github.com/wdl2459/ConQuR
###' https://wdl2459.github.io/ConQuR/ConQuR.Vignette.html
### devtools::install_github("wdl2459/ConQuR")

batchConQuR <- function(ampObj, corr, corr_refLvl, covars = NULL, pseudo_count = 1 ){
    # library(ConQuR)
    ##' conQuR expects normalized (non-integer) values, but simple log2(counts) or log3(relAbu) wont be enough
    ##'  Normalize Counts using TMM (edgeR; (log2 Counts Per Million transformation ))
    pseuMTX <- edgeR::DGEList(counts = ampObj$abund)                       ## 
    pseuMTX <- edgeR::calcNormFactors(pseuMTX)  # TMM normalization
    pseuMTX <- edgeR::cpm(pseuMTX, log = TRUE)  # Convert to log2-CPM
    
    metaMTX <- ampObj$metadata %>% as.matrix()
    
    corrMTX <- ConQuR::ConQuR(tax_tab     = pseuMTX,            # Feature table
                              batchid = metaMTX[[corr]],        # Batch variable
                              batch_ref = corr_refLvl,          # name of the reference level of batch
                              covariates = covars )             # covariates,  other than batch
    
    ## converting back to counts
    corrMTX <- 2^corrMTX - pseudo_count                                    ## convert back to count estimates
    corrMTX[corrMTX < 0] <- 0                                              ## replace negatives with 0
    corrMTX <- round(corrMTX, 0)                                           ## rounds to counts
    
    return(corrMTX)
}

###' ------------ limma:removeBatchEffect batch correction
###' https://rdrr.io/bioc/limma/man/removeBatchEffect.html
###' input is log2-transformed counts matrix that has been cleaned &/or rarified
###' I have set up the function for 1 factor to be corrected for (e.g. Batch)
###' the output is also log2-transformed matrix. 
###' I have created a conversion, converted counts are estimates

batchLimma<- function(ampObj, corr, pseudo_count = 1 ){
  # library(limma)
  pseuMTX <- log2( ampObj$abund + pseudo_count)                         ## needs log2 transfomraiton
  metaMTX <- ampObj$metadata[[corr]]
  corrMTX <- limma::removeBatchEffect(pseuMTX, metaMTX)                 ## returns fixed log-transformed matrix
  range(corrMTX) %>% print()                                            ## if range is from 0-~20, its a log2 transformation
  
  ## converting back to counts
  # maxVal <- max(corrMTX[is.finite(corrMTX)])
  # corrMTX[is.infinite(corrMTX) ] <- maxVal ## set up max value instead of Inf
  
  corrMTX <- 2^corrMTX - pseudo_count                                    ## convert back to count estimates
  corrMTX[corrMTX < 0] <- 0                                              ## replace negatives with 0
  corrMTX <- round(corrMTX, 0)                                           ## rounds to counts
  return(corrMTX)  
}


##' ----- batch correction explorations, do

###' 1) PERMANOVA, look for batch effect reduction
###' 2) PCoA plot, look for batch cluster reduction
###' 3) do boxplot of the abundance values (best if they are still log-transformed, but oh well)
###'      #  boxplot(cmp$abund, main = "mmuphinCorr", las =2) 
###'  ## see if the counts distribution looks similar to raw (good) or is way too squished
###' 4) check out range(cmp$abund), see if the range is similar to raw (should be)
###' 5) do exploreDepths on the batch-corrected. Is the rarecurve crazy now or looks beautiful


###' --------- BIAS correction pre DA analyses-------------------
###' info at: 
###; https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
###' and more at: https://github.com/FrederickHuangLin/ANCOMBC 
##' fix formula = fct$all
##' ANCOM-BC2 estimates and corrects both: 
##' 1) the sample-specific (sampling fraction), and 
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
  # BiocManager::install("microbiome") ## dependancy
  
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

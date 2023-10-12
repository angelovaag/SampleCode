#' -----------------------------------------------------------------------------------------------------------------------------
# --------------------------------------- my DiffAbu plotter function --------------------------------------
myDAplotter <- function(DA_obj, fcts,  rand=NULL, corr=NULL, transf="LOG", effZ=1, sig = 0.2, ntax=25, check=F, 
                        fc.color=NULL, fc.ncol=NULL, fc.nrow=NULL, colVector = myColors$discrete){
  
  DA_resMTX <- DA_obj$DA_resMTX
  amp_obj   <- DA_obj$amp_obj
  agg_at    <- DA_obj$agg_at
  DA_type   <- DA_obj$DA_type
  
  ##' flt &fmt DA_resMTX    
  if(DA_type=="MaAsLin2"){coef_col = "coef";           sig_col="qval";    metCols <- NULL } # "metadata"
  if(DA_type=="ALDEx2"  ){coef_col = "effect";         sig_col="p.adj";   metCols <- NULL }
  if(DA_type=="DESeq2"  ){coef_col = "log2FoldChange"; sig_col="padj";    metCols <- NULL }
  
  ### lets rename columns to keep them consistent
  colnames(DA_resMTX)[match(coef_col,colnames(DA_resMTX) ) ] <- "EffectSize"
  colnames(DA_resMTX)[match(sig_col, colnames(DA_resMTX))]   <- "Adj.Significance"
  # colnames(DA_resMTX)[match("metadata", colnames(DA_resMTX)) ] <- fcts[[1]]
  
  #### select significance thresholds
  taxCols <- c("TAX", "Taxon", "Class", "Phylm")
  DA.Cols <- c("EffectSize", "Adj.Significance") #EffectSize, qValue
  
  DA_resMTX.sig <- DA_resMTX %>%   filter(Adj.Significance < sig & abs(EffectSize) > effZ ) %>% arrange(desc(abs(EffectSize))) %>%
    relocate(TAX, .before=Taxon) %>%  dplyr::select(all_of(c(DA.Cols, taxCols))) #metCols,
  
  DA.sigTAX=row.names(DA_resMTX.sig) ##%>% sort()
  
  DA_resMTX.sig_kbl <-  unique(DA_resMTX.sig) %>% select(-c(TAX)) %>% 
    mutate(across(EffectSize, \(x) round(x, 2) ) ) %>%    #mutate(across(EffectSize, round, 2) ) %>% 
    mutate(Adj.Significance = sprintf("%.2e", Adj.Significance)) %>%  kable("html", digits=150) %>%  kable_classic_2(full_width=F) %>%
    kable_styling("striped") %>% scroll_box(width = "70%", height = "500px")
  
  info = c( DA_type ,"total differential TAXa ~", fcts[[1]],"(at ", agg_at, "level) : ",  dim(DA_resMTX.sig)[1], ".", '\n', 
            "Significance threshold < ", sig, '\n', 
            "Effect size threshold >", effZ, '\n')
  
  ### getting the flt & fmt DA res out (just in case here)
  res_obj <-  list(DA.sigTAX, DA_resMTX.sig, DA_resMTX.sig_kbl)
  objNames <- c("DA.sigTAXa", "DA.sig.tab", "DA.sig_kbl") 
  res_inf <-  list(info,  sig,   effZ,    ntax,    DA_type,   agg_at)
  infNames <- c("info", "sig", "effZ", "ntax", "DA_type", "agg_at")
  
  all <- do.call(c, list(res_obj, res_inf))
  names(all)<- c(objNames, infNames)
  
  ##' plotting volcanos
  vlc.title =  paste("TAXa changes ~", fcts[[1]])
  vlc.subtitle = DA_type
  vlc.caption = paste ("total = ", nrow(DA_resMTX), agg_at)
  
  library(EnhancedVolcano)
  volcano.plot <- EnhancedVolcano(DA_resMTX, x="EffectSize", y="Adj.Significance",  lab=DA_resMTX$Taxon,
                                  FCcutoff= effZ, pCutoff=sig, legendPosition = "bottom",
                                  caption = vlc.caption, captionLabSize = 10,
                                  ylab="-Log10(Adj.Significance)", titleLabSize = 12, title=vlc.title, xlab="EffectSize",
                                  subtitle = vlc.subtitle, legendLabSize=10 , axisLabSize = 10, labSize = 3); 
  # volcano.plot
  
  
  ##' --------  flt & fmt amp_obj
  
  ##' agglomerated amp used from DA_obj
        ##' agglomerate if needed 
        # if(agg_at!= "OTU"){
        #   print(paste("aggregating TAXa to", agg_at, "level"))
        #   source(paste0("~/Documents/MyApps/RScripts/myRscripts/dataClean.R"))
        #   amp_obj <- amp_agg_abund(amp_obj, agg_at = agg_at)
        #   # amp_obj
        #   # print(amp_obj$tax %>% head(12))
        # }
  
  counts <- amp_obj$abund
  tax    <- amp_obj$tax
  meta   <- amp_obj$metadata
  meta_sub <- meta %>% select(any_of(c(fcts, rand, corr))) %>% mutate_all(factor)
  
  if (transf=="LOG"){
    # the log2(abd) plot: log-ratios are captured
    txfCounts <- log((counts+1)/(rowSums(counts)+dim(counts)[2])) 
    txfCounts %<>% t() %>% as.data.frame() %>% select(DA_resMTX.sig$TAX) 
  }else{
    txfCounts <-  (counts+1)/(rowSums(counts)+dim(counts)[2]) 
    txfCounts %<>% t() %>% as.data.frame() %>% select(DA_resMTX.sig.sig$TAX)
  }
  
  
  ##' try DF4plot
  if(check==T){
    print("checkpoint DF4plot")
    print(DA_resMTX.sig %>% head) 
  }
  
  
  topSigTAX <- DA.sigTAX[1:ntax] %>% na.omit() %>% c()
  df4plot0<-NULL
  
  try({  
    for(i in topSigTAX){
      tmp<-NULL
      # meta_sub <- meta %>% select(any_of(c(fcts, rand, corr))) # is same as meta
      tmp<-data.frame(txfCounts[,i], meta_sub, rep(i),
                      rep(paste0(DA_resMTX.sig$Taxon[DA_resMTX.sig$TAX == i], " (", i,")" ) ),
                      # rep(paste0("EffectSize = ", sprintf("%.4g", top_sig_tax$effect[top_sig_tax$TAX== i]),
                      #            "; qval.BH = ", sprintf("%.1e", top_sig_tax$glm.eBH[top_sig_tax$TAX==i]) ) ),
                      rep( DA_resMTX.sig$EffectSize[DA_resMTX.sig$TAX== i]) , 
                      rep( DA_resMTX.sig$Adj.Significance[DA_resMTX.sig$TAX==i]) )
      if(is.null(df4plot0)){df4plot0<-tmp} else { df4plot0<-rbind(df4plot0,tmp)}
    }
    # print(tmp)
    
    #the value is log(relAbu) of TAXa in each Sample
    colnames(df4plot0)<-c("Value", c(fcts, rand, corr), "TAX" , "Taxon", "EffectSize", "Adj.Significance"); 
    df4plot0$Taxon %<>% factor( levels=unique(df4plot0$Taxon) )
  })
  
  
  
  ##' plotting FC/EffZ plots   
  # df4plot0$Taxon %>% unique()
  df4plot <- df4plot0 %>% mutate(SigVal=paste0("Eff.size=", formatC(EffectSize, digits=2), 
                                               "; pVal =", formatC(Adj.Significance, format="E", digits=2)) )
  
  # fc.color=NULL; fc.ncol=4; fc.nrow=NULL;
  if(is.null(fc.color)){fc.color = fcts[[1]] }
  fc.title = paste(DA_type, "top", ntax, "differential TAXa ( at",agg_at, "level) b/w", fcts[[1]], "groups")
  if (transf=="LOG"){ ylabel = "Log2(Relative abundance)" }else{ylabel = "Relative abundance" }
  
  ## info about the aes_string() replacement: https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
  foldchange.plot <- ggplot(unique(df4plot),  aes_string(fcts[[1]], "Value", color= fc.color )) + # #aes(vars(fcts[[1]]), "Value", color= fc.color )
    ylab(ylabel) +    labs(title=fc.title,   color=fc.color) +
    geom_boxplot(outlier.size = 0) +  xlab(fcts[[1]]) +
    geom_jitter(position = position_jitterdodge(), alpha=0.3) +
    facet_wrap( ~ Taxon + SigVal, scales="free_x", #labeller = label_wrap_gen(width=45),
                     ncol=fc.ncol, nrow=fc.nrow  ) + theme_bw() + scale_color_manual(values=colVector) +
    theme( axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5, size = 14),  legend.position = "top",
          strip.text.x = element_text(size = 12, colour = "black")) #, angle = 90
  foldchange.plot
  
  res_plt<- list(volcano.plot, foldchange.plot, df4plot)
  pltNames <- c("VolcanoPlot", "EffSizePlot", "DF4plot")
  
  all<- do.call(c, list(res_plt, res_obj, res_inf))
      names(all) <- c( pltNames, objNames, infNames)
  
  return(all)
}








#' -----------------------------------------------------------------------------------------------------------------------------
# ----------------------------- my DESeq2 function --------------------------------------
##' e.g. use: DA_objDSQ <- myDESeq2(ampHPV, fcts = fcts, agg_at = "OTU")
##' For the correction design: **add the biological factors to the fml(fcts),**
##' ** but the stat confounders (demographics, batches, etc) go in the correction factor**

myDESeq2<-function(amp_obj, fcts, corr=NULL, agg_at= "OTU", check=F, reffLvL=NULL){
  
  ##' agglomerate if needed
  if(agg_at!= "OTU"){
    print(paste("aggregating TAXa to", agg_at, "level"))
    source(paste0("~/Documents/MyApps/RScripts/myRscripts/dataClean.R"))
    amp_obj <- amp_agg_abund(amp_obj, agg_at = agg_at)
    # amp_obj
    # print(amp_obj$tax %>% head(12))
  }
  
  counts  <- amp_obj$abund
  tax     <- amp_obj$tax 
  meta    <- amp_obj$metadata
  meta_sub <-  meta %>% select(any_of(c(fcts, corr))) %>% mutate_all(factor)
  # if(is.null(ref_lvl)){reff=NULL}else{ reff <- paste0( ref_lvl , ",", fcts[[1]] ) }
  # reff
  
  #'------------------------ data PREP ---------------------------------------------
  countdata = round(as(counts, "matrix"), digits = 0) +1
  DA_type="DESeq2"
  #'------------------------ DDS calculation ------------------------
  if(is.null(corr)){
    # fml <- as.formula(paste(" ~",  fcts[[1]]) ) ## paste(fcts, collapse = "+") ) # !!!!!!!!!!!!!!!!!!!
    fml <- as.formula(paste(" ~", paste(rev(fcts), collapse="+")))
    redfml <- as.formula(paste("~1") )
  }else{
    fml <- as.formula(paste(" ~ ", paste(rev(corr), collapse = " + "), "  +  ", paste(rev(fcts), collapse="+")  ) ) # paste(rev(fcts), collapse = " + "))
    redfml <- as.formula(paste("~1 + ", paste(rev(corr), collapse  = "+") ) )
  }
  library(DESeq2)
  library(BiocParallel)
  threads=detectCores()-4
  register(MulticoreParam(threads))
  
  dds <- DESeqDataSetFromMatrix(countdata, meta_sub, fml)
  dsqRAW = DESeq(dds, test="LRT", reduced = redfml, parallel = T);
  # dispPlot <- DESeq2::plotDispEsts(dsqRAW)
  #' -------------------------- extracting DDS results ---------------------------------------------
  #' If results is run without specifying contrast or name (name is for continuous variables), \
  #' it will return the comparison of the last level of the last variable in the design formula \
  #' over the first level of this variable. (log2(lastGroup/firstGroup) => 
  #' **firstGroup listed in the levels of a factor, will be considered RefflvL (the Control/Negative/untreated))** 
  #' To set contrast manually, set contrast=c("condition","treated","untreated.Control.Negative")
  dsqRES = results(dsqRAW, cooksCutoff = F)#, contrast = reffLvL) 
  dsqRES <- dsqRES[order(dsqRES$padj),];
  
  # if(check==T){  print(head(res, 4)) }
  
  #' ----------------------- adding TAXnames to dsqRES ---------------------------------------------
  if(agg_at=="OTU"){lRank="Species"}else{lRank=agg_at}
  dsqRES_tax = cbind(as.data.frame(dsqRES) , # as.matrix(countdata[rownames(dsqRES), ]), 
                  Taxon = tax[[lRank]][match(rownames(dsqRES), rownames(tax) ) ],
                  Class =    tax$Class[match(rownames(dsqRES), rownames(tax) ) ] ,
                  Phylm =   tax$Phylum[match(rownames(dsqRES), rownames(tax) )] )
  if(check==T){
    print(fml); print(redfml)
    print(head(dsqRES, 4))
    y<- ncol(dsqRES_tax); y1=y-5
    print(dsqRES_tax[1:5, y1:y] )
  }
  
  ##' --- output formatting
  dsqRES_tax <- dsqRES_tax %>% arrange(desc(abs(log2FoldChange))) %>% 
    mutate(TAX=rownames(.)) %>% relocate(TAX, .before="Taxon")
  
         all<- list(dsqRES_tax, agg_at,  DA_type,   amp_obj,   fcts, corr,       dsqRAW, dsqRES)
  names(all) <- c("DA_resMTX", "agg_at", "DA_type", "amp_obj", "fcts", "corr", "dsqRAW", "dsqRES")
  
  return(all)
}









#' -----------------------------------------------------------------------------------------------------------------------------
# -------------------------------------- wrapper Maaslin function -------
##' e.g. use DA_objMSN<- wrpMaAsLin2(ampHPV, fcts, maasDIR=DIR, model="LM", agg_at = "Genus", rand = rand)

 
wrpMaAsLin2 <- function(amp_obj, fcts, rand = NULL, maasDIR= NULL, min_prevalence=0.1, min_abundance=0.0,
                        reff_lvl=NULL,  agg_at = "OTU", model = "all"){
  
  ##' agglomerate if needed
  if(agg_at!= "OTU"){
    print(paste("aggregating TAXa to", agg_at, "level"))
    source(paste0("~/Documents/MyApps/RScripts/myRscripts/dataClean.R"))
    amp_obj <- amp_agg_abund(amp_obj, agg_at = agg_at)
    # amp_obj
    # print(amp_obj$tax %>% head(12))
  }
  
  counts <- amp_obj$abund
  tax    <- amp_obj$tax
  meta   <- amp_obj$metadata
  meta %<>% select(any_of(c(fcts, rand))) %>% mutate_all( factor)
  
  
  library(Maaslin2)
  DA_type <- "MaAsLin2"
  
  if(model %in% c("LM", "CPLM","NEGBIN","ZINB") ){norm=NULL; transf=NULL }  # wrapper will auto do all combinations 
  if(model=="all"){ model=NULL; norm=NULL; transf=NULL  } # wrapper will auto do all combinations 
  
  
  #' ------------------ the wrapper code --------------------
  
  # print(list(model, norm, transf) )
  
  
  source(paste0(myRpath, "otherFunctions.R"))
  phy <- amp2phy(amp_obj)
  
  source(paste0(myRpath, "maaslin_wrapper_AA.R"))
  start.time <- Sys.time()
  suppressMessages( suppressWarnings( {
    maasOUT = multimodel_maaslin(phy=phy, fixed_effects=fcts, output_dir=maasDIR, random_effects= rand, 
                                 min_abundance = min_abundance ,  min_prevalence = min_prevalence, #needs % at >2 samples
                                 standardize=F, model=model, txf = transf, norm = norm)
  }) ) 
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time); print(time.taken)
  
  # print(maasOUT %>% head() )
  
  ##' --- output formatting

  ##' feature to TAX & rownames
  maasOUT <- maasOUT %>% filter(metadata %in% fcts[[1]]) %>% select(-all_of(c("metadata", "value", "name") ) ) %>%  
            mutate(TAX=feature) %>% relocate(TAX, .after=last_col() )
  ##' in cases fct[main] has 2+ categories, we get pair-wise comparisons to REFF_level => duplicated features => 
  ##' => problem with assigning row names => 
  ##' will extract only distinct(feature) with higher abs(coef)
  maasOUT <- maasOUT %>% arrange(desc(abs(coef)) ) %>% distinct(feature, .keep_all = TRUE) %>% 
              column_to_rownames("feature") %>% arrange(qval) 

  
  
  if(agg_at=="OTU"){lRank="Species"}else{lRank=agg_at}
  # maasOUTt = cbind(maasOUT,
  #                  Taxon  = tax[[lRank]][match(rownames(maasOUT), rownames(tax) ) ],
  #                  Class  =    tax$Class[match(rownames(maasOUT), rownames(tax) ) ],
  #                  Phylm  =   tax$Phylum[match(rownames(maasOUT), rownames(tax) ) ] )

  maasOUTt = cbind(maasOUT,
                   Taxon  = tax[[lRank]][match(maasOUT$TAX, rownames(tax) ) ],
                   Class  =    tax$Class[match(maasOUT$TAX, rownames(tax) ) ],
                   Phylm  =   tax$Phylum[match(maasOUT$TAX, rownames(tax) ) ] )
  
  all <- list(maasOUTt, agg_at,  DA_type, amp_obj, fcts, rand)
  names(all) <- c("DA_resMTX", "agg_at", "DA_type", "amp_obj", "fcts", "rand")
  
  return(all)
} 


#' ---------------------------------------------------------------------------------------------------------------------------
# --------------- my Maaslin funciton ----------------------------------------------------------

myMaAsLin2 <- function(amp_obj, fcts, rand = NULL, maasDIR= "maasOUT", minPrev=0.1, minAbu=0.0, 
                       reffLvL=NULL,  agg_at = "OTU", model = "LM", transf="LOG", norm="CLR"){
  
  ##' agglomerate if needed
  if(agg_at!= "OTU"){
    print(paste("aggregating TAXa to", agg_at, "level"))
    source(paste0("~/Documents/MyApps/RScripts/myRscripts/dataClean.R"))
    amp_obj <- amp_agg_abund(amp_obj, agg_at = agg_at)
    # amp_obj
    # print(amp_obj$tax %>% head(12))
  }
  
  counts <- amp_obj$abund
  tax    <- amp_obj$tax
  meta   <- amp_obj$metadata
  meta %<>% select(any_of(c(fcts, rand))) %>% mutate_all( factor)
  
  
  library(Maaslin2)
  DA_type <- "MaAsLin2"
  ####### ---- direct maaslin run ---
  ##' To use reffLvL do: revLvL = paste0('FCT[[1]], REFlvl'). e.g. **reffLvL=paste0("HPVstatus,HPVneg")**
  # met$Groups=paste(met[[fxd_effects[[1]]]], met[[fixed_effects[[2]]]], sep="-")
  # if(!is.null(reff_lvl)){ reff_lvl <- c(paste0(fcts[[1]], ",", reff_lvl ) ) }
  # threads=detectCores()-4
  # register(MulticoreParam(threads))
  
  print("using 1 model maaslin ")
  
  
  start.time <- Sys.time()
  # suppressMessages( suppressWarnings(
  maasOUTa = Maaslin2::Maaslin2(input_data = counts, input_metadata = meta, output = maasDIR,
                                fixed_effects = fcts,
                                random_effects = rand,
                                reference= reffLvL,
                                min_abundance = minAbu,
                                min_prevalence = minPrev, #neeeds at least 2 samples, so 10% is safest
                                normalization= norm, transform= transf, analysis_method = model,
                                plot_heatmap=F, plot_scatter=F, save_models = T, cores=1) ###keep cores at 1, fastest!!
  # ))
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time); print(time.taken)
  
  #' Add Poorani's AIC function (June 12)
      AIC_forcplm <- function(fit) {
        #' extracts aic value from slot of fit
        #' is of class cpglm or uses AIC() otherwise
        #' @param fit model fit
        #' @return AIC value
        if(inherits(fit, 'cpglm')) {
          return(fit$aic)
        } else{
          return(AIC(fit))
        }
      }
      
  AICs <- data.frame(AIC=unlist(lapply(maasOUTa[5]$fits, function(x) tryCatch(AIC_forcplm(x), error=function(e) NA))))
  maasOUT <- merge(maasOUTa$results, round(AICs, 2), by.x="feature", by.y=0)
  
  ##' --- output formatting
  ##' feature to TAX & rownames
  maasOUT <- maasOUT %>% filter(metadata %in% fcts[[1]]) %>% select(-all_of(c("metadata", "value", "name") ) ) %>%  
    mutate(TAX=feature) %>% relocate(TAX, .after=last_col() )
  ##' in cases fct[main] has 2+ categories, we get pair-wise comparisons to REFF_level => duplicated features => 
  ##' => problem with assigning row names => 
  ##' will extract only distinct(feature) with higher abs(coef)
  maasOUT <- maasOUT %>% arrange(desc(abs(coef)) ) %>% distinct(feature, .keep_all = TRUE) %>% 
    column_to_rownames("feature") %>% arrange(qval) 
  
  
  
  if(agg_at=="OTU"){lRank="Species"}else{lRank=agg_at}
  maasOUTt = cbind(maasOUT,
                   Taxon  = tax[[lRank]][match(rownames(maasOUT), rownames(tax) ) ],
                   Class  =    tax$Class[match(rownames(maasOUT), rownames(tax) ) ],
                   Phylm  =   tax$Phylum[match(rownames(maasOUT), rownames(tax) ) ] )
  
  all <- list(maasOUTt, agg_at,  DA_type, amp_obj, fcts, rand)
  names(all) <- c("DA_resMTX", "agg_at", "DA_type", "amp_obj", "fcts", "rand")
  
  return(all)
}  



#' -----------------------------------------------------------------------------------------------------------------------------
# ---------------------------------- ALDEX function (fancy) --------------------------------------
##' can do 1 factor with 2 levels in KW or GLM mode (GLM mode depending on sufficient contrast)
##' can do 1 factor with 2+ levels in GLM mode
##' can do 1 factor with >=2 levels with adjustments (corr)
##' is slow even in multi-thread mode, but KW is fastest
##' can agglomerate to TAX level chosen
##' fct1 should be pre-set to list CTRL level first in list of factors
##' works with myDAplotter

myALDEx2 <- function(amp_obj, fct1, corr= NULL, agg_at = "OTU", mode="glm", check=F, denom="all", threads=NULL){
  #ensure fct1 & corr levels are listed with CTRLlvl listed first
  
  ##' agglomerate if needed
  if(agg_at!= "OTU"){
    print(paste("aggregating TAXa to", agg_at, "level"))
    source(paste0("~/Documents/MyApps/RScripts/myRscripts/dataClean.R"))
    amp_obj <- amp_agg_abund(amp_obj, agg_at = agg_at)
    # amp_obj
    # print(amp_obj$tax %>% head(12))
  }
  
  meta   <- amp_obj$metadata # %>% mutate_at(., c(fct1, corr), factor) #factorl lvls should be pre-made, CTRL lvl list 1st
  counts <- amp_obj$abund
  tax    <- amp_obj$tax
  
  ###' ------------ settings
  ###' for KW mode                            ## only for 1 factor with 2 levels (quicker)
  if(mode == "kw") { #only for 1 factor with 2 levels
    fml <- as.formula(paste0("~", fct1))
    if(denom!="all"){denom="iqlr"}
    mm  <- meta[[fct1]]
    if(check==T){print(levels(mm))}
  }
  ###' ------------ settings
  ###' for GLM mode                           ## for 1 main factor with multiple levels & corrections
  if(mode== "glm"){
    if(is.null(corr) ){   # main factor with 2+, but not corrections
      fml <- as.formula(paste("~ ", fct1) ) 
    }else{                 # main factor with 2+ levels & correction
      fml <- as.formula(paste("~ ", fct1, " + ", paste(corr, collapse="+") )) 
    }
    
    denom = "all"   ## iqrl not applicable in GLM mode
    refLvL= "RefLvL"
    mm  <- model.matrix(fml, data=meta) #assumes refLvL listed first
  }
  

  DA_type="ALDEx2"
  library(ALDEx2)
  library(BiocParallel)
  if(is.null(threads)) {threads=detectCores()-4 }
  register(MulticoreParam(threads))
  
  xx <- aldex.clr(counts,           # the rows are SampleNames and Columns are TAXa
                  mm,               # the model matrix
                  mc.samples=128,   # monte carlo (random) samples to use for estimating distribution: 16 or [def: 128]
                  denom=denom,      # only all is accepted for glm
                  verbose=T, useMC = T)
  # row.names(xx@analysisData$DNA103)
  
  ###' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###' for KW mode  #only for 1 factor with 2 levels
  if(mode == "kw"){  
    x.kw  <- aldex.kw(xx, verbose = T)
    x.eff <- aldex.effect(xx, verbose = T) 
    x.all <- data.frame(x.kw, x.eff)
    x.all %<>% arrange(desc(abs(effect))) %>% dplyr::select(any_of(c("effect", "glm.eBH"))) %>% 
      mutate(TAX=rownames(.)) %>% dplyr::rename(p.adj = "glm.eBH")
    if(check==T){ print(x.all %>% head(20) ) }
  }
  ###' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###' for GLM mode
  if(mode == "glm"){
    print("starting aldex GLM")
    # Buggy function aldex.glm(). 
    # Common error: "missing value where TRUE/FALSE needed" @ "verbose[1] == TRUE"?
    # Has smth to do with no insufficient contrasts bw groups?  Try lowering agg_at level or seqDepth
    xx.glm <-aldex.glm(xx, mm) #do not input verbose (errors). Cannot useMC.
    
    # colnames(xx.glm)
    xx.glm1 <- xx.glm %>% dplyr::select(contains(fct1) )  %>%                 ## as w Maaslin we select only the res of main factor 
      dplyr::select(contains("pval.holm")) %>%                        ## only pval is impotant. No effect GLM
      select(-contains("Intercept"))%>% rename_all(~str_replace(., fct1, paste0(refLvL, "~") ) ) 
    if(check==T){ colnames(xx.glm1) %>% print() } 
    
    print("starting aldex effect")
    xx.eff <- aldex.glm.effect(xx, verbose = TRUE, useMC = TRUE); 
    # xx.eff$Group1HPVcarc #%>% conames()
    
    xx.eff1 =  data.frame(xx.eff) %>% dplyr::select(contains("effect")) %>% ## only effect is collected here
      dplyr::select(contains(fct1))  %>%                              ## as w Maaslin we select only the res of main factor
      rename_all(~str_replace(.,fct1, paste0(refLvL, "~")));
    
    if(check==T){  colnames(xx.eff1)  %>% print() }
    
    xx.all = data.frame(xx.glm1, xx.eff1) # we have already pre-select only the res of main factor
    if(check==T){ xx.all %>% head(5)%>% print() }
    
    #creating new column with the value.which = max(abs(effect)) & another column with value.which = min(qval)
    xx.all1 <- xx.all %>% mutate(effect  = apply(dplyr::select(., contains("effect")), 1, \(x)  x[which.max(abs(x))] ))  %>%
      mutate(p.adj = apply(dplyr::select(., contains(".holm")), 1, min)) %>%
      select(c("effect", "p.adj")) 
    
    x.all <- xx.all1 %>% arrange(desc(abs(effect))) %>% mutate(TAX=rownames(.)) 
    if(check==T){ xx.all1 %>% head(5) %>% print() }
  }  
  
  if(agg_at=="OTU"){lRank="Species"}else{lRank=agg_at}
  
  x.all = cbind(x.all,
                Taxon  = tax[[lRank]][match(rownames(x.all), rownames(tax) ) ],
                Class  =    tax$Class[match(rownames(x.all), rownames(tax) ) ],
                Phylm  =   tax$Phylum[match(rownames(x.all), rownames(tax) ) ] )
  
  
  all <-         list(x.all, agg_at,  DA_type, amp_obj, fct1, corr)
  names(all) <- c("DA_resMTX", "agg_at", "DA_type", "amp_obj", "fcts", "corr")
  
  return(all)      
  
}
#'------------


#' ---------------------------------------------------------------------------------------------------
# ---------------------------------- Compare lists function  ---------------------- 
##'  to compare output of sigABU lists
list_compare <-function(list1, list2, taxa_list, agg_at="Species"){
  comTax<-try({Reduce(intersect, list(list1, list2) ) })
  comNam<- taxa_list[[agg_at]][row.names(taxa_list) %in% comTax]
  common <- list(comTax, comNam)
  names(common) <- c("comTax", "comNam")
  return(common)
} 

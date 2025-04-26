#' ---------------------------------------------------------------------------------------------------------
# --------------------------------------- my DiffAbu plotter function --------------------------------------
myDAplotter <- function(DA_obj, fcts,  rand=NULL, corr=NULL, norm="LogRel", norm_add = 1,  
                        effZ=1, sig = 0.2, ntax=25, check=F, fltONsig=T,  mkLongiPlot = F, 
                        fc.color=NULL, fc.ncol=NULL, fc.nrow=NULL, colVector = myColors$discrete, 
                        kble_width="100%", kble_heigth="500px", volcLabSize = 3, x.angle = 0 ){
  
  DA_resMTX <- DA_obj$DA_resMTX
  agg_at    <- DA_obj$agg_at
  DA_type   <- DA_obj$DA_type
  amp_obj   <- DA_obj$amp_obj
  counts <- amp_obj$abund
  tax    <- amp_obj$tax
  meta   <- amp_obj$metadata %>% select(any_of(c(fcts, rand, corr))) %>% mutate_all(factor)
  
  ##' flt &fmt DA_resMTX    
  if(DA_type=="MaAsLin2"){coef_col = "coef";           sig_col="qval";    metCols <- NULL; fit = "glm" } 
  if(DA_type=="ALDEx2"  ){coef_col = "effect";         sig_col="p.adj";   metCols <- NULL; fit = "glm" }
  if(DA_type=="DESeq2"  ){coef_col = "log2FoldChange"; sig_col="padj";    metCols <- NULL; fit = "glm" }
  if(DA_type=="lmerDA"  ){coef_col = "EffSiZ";         sig_col="padj";    metCols <- NULL; fit = "glm" }
  if(DA_type=="bamDA"   ){coef_col = "estimate";       sig_col="p.adj";   metCols <- NULL; fit = "gam" }
  
  ### lets rename columns to keep them consistent
  colnames(DA_resMTX)[match(coef_col,colnames(DA_resMTX) ) ] <- "EffectSize"
  colnames(DA_resMTX)[match(sig_col, colnames(DA_resMTX))]   <- "Adj.Significance" #long name needed for spacing
  # colnames(DA_resMTX)[match("metadata", colnames(DA_resMTX)) ] <- fcts[[1]]
  
  #### select significance thresholds
  taxCols <- c("TAX", "Taxon", "Class", "Phylm")
  DA.Cols <- c("EffectSize", "Adj.Significance") #EffectSize, qValue
  
  if (fltONsig == T ){
    DA_resMTX.sig <- DA_resMTX %>%   filter(Adj.Significance < sig & abs(EffectSize) > effZ ) %>% 
      arrange(desc(abs(EffectSize))) %>%
      relocate(TAX, .before=Taxon) %>%  dplyr::select(all_of(c(DA.Cols, taxCols))) #metCols,
  }else{
    DA_resMTX.sig <- DA_resMTX %>%   filter(abs(EffectSize) > effZ ) %>% 
      arrange(desc(abs(EffectSize))) %>%
      relocate(TAX, .before=Taxon) %>%  dplyr::select(all_of(c(DA.Cols, taxCols))) #metCols,
  }
  
  DA.sigTAX=row.names(DA_resMTX.sig) ##%>% sort()
  
  DA_resMTX.sig_kbl <-  unique(DA_resMTX.sig) %>% select(-c(TAX)) %>% 
    mutate(across(EffectSize, \(x) round(x, 2) ) ) %>%    #mutate(across(EffectSize, round, 2) ) %>% 
    mutate(Adj.Significance = sprintf("%.2e", Adj.Significance)) %>%  kable("html", digits=150) %>% 
    kable_classic_2(full_width=F) %>%
    kable_styling("striped") %>% scroll_box(width = kble_width, height = kble_heigth)
  
  info = c( DA_type ,"total differential TAXa ~", fcts[[1]],"(at ", agg_at, "level) : ",  
            dim(DA_resMTX.sig)[1], ".", '\ \n', 
            "Significance threshold < ", sig, '\ \n', 
            "Effect size threshold >", effZ, '\ \n')
  
  ### getting the flt & fmt DA res out (just in case here)
  res_obj <-  list(DA.sigTAX, DA_resMTX.sig, DA_resMTX.sig_kbl)
  objNames <- c("DA.sigTAXa", "DA.sig.tab", "DA.sig_kbl") 
  res_inf <-  list(info,  sig,   effZ,    ntax,    DA_type,   agg_at)
  infNames <- c("info", "sig", "effZ", "ntax", "DA_type", "agg_at")
  
  all <- do.call(base::c, list(res_obj, res_inf)) # can never use c<-"variable" anymore :(
  names(all)<- c(objNames, infNames)
  
  ##' plotting volcanos
  vlc.title =  paste("TAXa changes ~", fcts[[1]])
  vlc.subtitle = DA_type
  vlc.caption = paste ("total = ", nrow(DA_resMTX), agg_at)
  
  library(EnhancedVolcano)
  if (fltONsig == T){
    volcano.plot <- EnhancedVolcano(DA_resMTX, x="EffectSize", y="Adj.Significance", 
                                    lab=DA_resMTX$Taxon,
                                    FCcutoff= effZ, pCutoff=sig, legendPosition = "bottom",
                                    caption = vlc.caption, captionLabSize = 10,
                                    ylab="-Log10(Adj.Significance)", titleLabSize = 12, 
                                    title=vlc.title, xlab="EffectSize",
                                    subtitle = vlc.subtitle, legendLabSize=10 ,
                                    axisLabSize = 10, labSize = volcLabSize); 
    # volcano.plot
  }else{
    volcano.plot <- NULL
  }
  
  ##' ------------ binomial plot (to be used only in 2-level comparisons)
  DA_resMTX.sig$log2fc <- log2(exp(DA_resMTX.sig$EffectSize))
  # print(DA_resMTX.sig)
  binom.xlab=paste("log(2)FC abundance")
  # binom.title = paste(DA_type, "top", ntax, "differential TAXa ( at",agg_at, "level) b/w", fcts[[1]], "groups")
  
  tax_fill="Class"
  sel_taxa <- c(unique(DA_resMTX.sig[,tax_fill]))
  sel_colors <- rev(colVector)[1:length(sel_taxa)]
  names(sel_colors) <- sel_taxa
  # print(sel_colors)
  
  
  DA_resMTX.sig.pos <- DA_resMTX.sig %>% top_n(ntax, abs(log2fc)) #top SIG +&-
  DA_resMTX.sig.neg <- DA_resMTX.sig %>% filter(log2fc > 2.5) %>% top_n(5, log2fc)## grab some + values
  DA_resMTX.sig.bin <- rbind(DA_resMTX.sig.pos, DA_resMTX.sig.neg) %>% arrange(log2fc) %>%
    distinct(TAX, .keep_all = TRUE) # removes rbind()  redundancy
  # print(DA_resMTX.sig.bin)
  # print(DA_resMTX.sig.neg)
  # print(DA_resMTX.sig.bin)
  
  binom.plot <- ggplot(DA_resMTX.sig.bin, aes(x = log2fc, y = Taxon, fill=Class)) +
    theme_bw() + scale_fill_manual(values = sel_colors, name=tax_fill) + # Set custom colors
    geom_bar(stat = "identity") + # Remove the legend
    xlab(binom.xlab) + ylab("Feature") # + labs(title=binom.title)
  
  
  
  ##' --------  normalize abundance of raw counts
  # counts <- amp_obj$abund
  source(paste0(myRpath, "dataCorrections.R"))
  nrmCounts <- normMyData(countsMTX = counts, norm_method = norm, pseudocount = norm_add)
  if(check==T){   print("---- checkpoint0.5 dataNorm")   ; 
    print(nrmCounts[[1]][1:10, 1:10] %>% t() %>% as.data.frame()); print(DA_resMTX.sig$TAX)  }
  txfCounts <- nrmCounts[[1]] %>% t() %>% as.data.frame() %>% select(DA_resMTX.sig$TAX)
  ylabel    <- nrmCounts[[2]]; 
  
  
  ##' --------  boxplot of EffectSize (DF4plot)
  if(check==T){ print("---- checkpoint1 DF4plot"); print(DA_resMTX.sig %>% head) ; print(DA.sigTAX) }
  
  topSigTAX <- DA.sigTAX[1:ntax] %>% na.omit() %>% c()
  df4plot0<-NULL
  
  # try({  
  for(i in topSigTAX){
    tmp<-NULL
    tmp <- data.frame(txfCounts[,i], meta, rep(i),
                      rep(paste0(DA_resMTX.sig$Taxon[DA_resMTX.sig$TAX == i], " (", i,")" ) ),
                      rep( DA_resMTX.sig$EffectSize[DA_resMTX.sig$TAX== i]) , 
                      rep( DA_resMTX.sig$Adj.Significance[DA_resMTX.sig$TAX==i]) )
    # rep(paste0("EffectSize = ", sprintf("%.4g",
    # top_sig_tax$effect[top_sig_tax$TAX== i]), \n"; qval.BH = ", 
    # sprintf("%.1e", top_sig_tax$glm.eBH[top_sig_tax$TAX==i]) ) ),
    if(is.null(df4plot0)){df4plot0<-tmp} else { df4plot0<-rbind(df4plot0,tmp)}
  }
  # print(df4plot0)
  
  if(!is.null(df4plot0)){ ## avoids errors 
    colnames(df4plot0)<-c("Value", c(fcts, rand, corr), "TAX" , "Taxon", 
                          "EffectSize", "Adj.Significance"); 
    df4plot0$Taxon %<>% factor( levels=unique(df4plot0$Taxon) )
  }
  # })
  
  ##' --------  EffectSize boxplot for sigTAXa 
  if(is.null(df4plot0)){ df4plot <- NULL; foldChangeBox.plot <- NULL }else{
    df4plot <- df4plot0 %>% dplyr::mutate(SigVal=paste0("Eff.size=", formatC(EffectSize, digits=2),
                                                        "; pVal =", formatC(Adj.Significance,
                                                                            format="E", digits=2)) )
    df4plot$Value %>% range() %>% min() -> minflt ## grabbing the minimum post-norm value (representing abundance == 0)
    # df4plot %<>% filter(!Value == minflt) ## removing samples with abundance == 0
    df4plot %<>% unique()
    
    ##' try DF4plot
    if(check==T){
      print("---- checkpoint2 DF4plot w meta")
      print(df4plot %>% head) 
      print(paste("min abundance is:", minflt, "; ylabel is:", ylabel))
      
    }
    
    # fc.color=NULL; fc.ncol=4; fc.nrow=NULL;
    if(is.null(fc.color)){fc.color = fcts[[1]] }
    fc.title = paste(DA_type, "top differential TAXa (", agg_at, "level)")
    fc.subtl = paste("Contrasts  b/w", fcts[[1]], "groups")
    fc.captn = ifelse(dim(DA_resMTX.sig)[1] > ntax, paste("Plotted:", ntax, "/" , dim(DA_resMTX.sig)[1]) , "")
    # fc.captn = paste("Plotted:", ntax, "/" , dim(DA_resMTX.sig)[1])
    
    myTheme <- theme_bw() + theme(
      # axis.text.x = element_blank(), axis.ticks.x = element_blank(),   ## facetes x.axis settings
      axis.text.x = element_text(angle= x.angle, hjust=0.5, size=10),              ## facetes x.axis settings
      legend.position = "top",  legend.title= element_text(size=12, face = "bold"), legend.text = element_text(size = 10), 
      axis.title = element_text(size = 12),  ## global axises settings
      strip.text.x = element_text(size = 10, colour = "black"),  ## facet text strip settings
      plot.title = element_text(size=14, face = "bold"), plot.subtitle = element_text(size = 10), 
    )
    myLabsBOX <- labs(y=ylabel, x=fcts[[1]], color=fc.color , title = fc.title, 
                      subtitle = fc.subtl, caption = fc.captn)
    if (length(fcts) > 1){
    myLabsLIN <- labs(y=ylabel, x=fcts[[2]], color=fc.color , title = fc.title, 
                      subtitle = fc.subtl, caption = fc.captn) }
    
    ## info about the aes_string() replacement: https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
    foldChangeBox.plot <- ggplot(unique(df4plot),  aes_string(fcts[[1]], "Value", color= fc.color )) + 
      facet_wrap( ~ Taxon + SigVal, scales="free_x", ncol=fc.ncol, nrow=fc.nrow  ) +  ###labeller = label_wrap_gen(width=45),
      geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.3) +
      geom_boxplot(outlier.shape = NA, fill = NA) + 
      myLabsBOX + scale_color_manual(values = colVector) + myTheme 
  }
  
  
  ##' --------  linear regression plot (for continuous (or longitudinal) fct2)
  if ( length(fcts) >1 && mkLongiPlot == T) {
    df4plot[[fcts[[2]]]] %<>% as.character() %>% as.numeric()
    longi.plot <- ggplot(df4plot, aes_string(fcts[[2]], "Value")) +
      facet_wrap(~Taxon+SigVal, scales="free", ncol=fc.ncol, nrow=fc.nrow ) +
      geom_point(aes_string(color=fcts[[1]]),  size = 1, alpha = 0.3, 
                 position=position_jitterdodge(jitter.width = 0.3))  + 
      geom_smooth(aes_string(color=fcts[[1]]), method = fit, alpha= 0.1) + ##gam for non-linear fit, glm for lin fit
      myLabsLIN + scale_color_manual(values = colVector) + myTheme 
    
  } else { longi.plot=NULL }
  
  
  res_plt<- list(volcano.plot,  foldChangeBox.plot, binom.plot,   longi.plot ,  df4plot)
  pltNames <- c("VolcanoPlot",    "fcBoxPlot",     "binomPlot",   "longiPlot", "df4plot")
  
  all <- do.call(base::c, list(res_plt, res_obj, res_inf))
  names(all) <-       c( pltNames, objNames, infNames)
  
  return(all)
}

# dfi4plot <- desiP$DF4plot %>% unique() ##  removes redundant zero 
# # dfi4plot %>% dim
# `as.numeric(as.character(Mom_PARITYcn))` should be as.character() %>% as.numeric() beforehand or wrong
# p1 <- ggplot(dfi4plot, aes(x=as.numeric(as.character(Mom_PARITYcn)), y=Value)) +
#   facet_wrap(~Taxon+SigVal, scales="free", nrow=4) +
#   geom_smooth(aes(color= Mom_prpmABXc), method = "lm", alpha= 0.1) +
#   geom_point(aes(color=Mom_prpmABXc), position=position_jitterdodge(jitter.width = 0.3))  +
#   theme_bw() +  labs(y="Log2(RelAbundance)", x="Mom_PARITYcn", color="Mom_prpmABX" ) +
#   theme(text=element_text(size=16), axis.text.x = element_text(angle=45, hjust=1, size=14),
#         legend.position = "bottom",  legend.title= element_text(size=12), #legend.position = "none"
#         legend.text=element_text(size=10)) +
#   scale_color_manual(values = myColors$discrete[c(30,22)]); p1





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
  meta_sub <-  meta %>% select(any_of(c(fcts, corr)))# %>% mutate_all(factor)
  # if(is.null(ref_lvl)){reff=NULL}else{ reff <- paste0( ref_lvl , ",", fcts[[1]] ) }
  # reff
  
  #'------------------------ data PREP ---------------------------------------------
  countdata = round(as(counts, "matrix"), digits = 0) +1
  DA_type="DESeq2"
  #'------------------------ DDS calculation ------------------------
  if(is.null(corr)){
    fml <- as.formula(paste(" ~", paste(rev(fcts), collapse="+")))
    redfml <- as.formula(paste("~1") )
  }else{
    fml <- as.formula(paste(" ~ ", paste(rev(corr), collapse = " + "), "  +  ", paste(rev(fcts), collapse="+")  ) ) # paste(rev(fcts), collapse = " + "))
    redfml <- as.formula(paste("~ 1+", paste(rev(corr), collapse  = "+") ) ) #~1 doesnt change, but safer with
  }
  library(DESeq2)
  library(BiocParallel)
  threads=detectCores()-4
  register(MulticoreParam(threads))
  
  dds <- DESeqDataSetFromMatrix(countdata, meta_sub, fml)
  set.seed(123)
  dsqRAW = DESeq(dds, test="LRT", reduced = redfml, parallel = T);
  # dispPlot <- DESeq2::plotDispEsts(dsqRAW)
  #' -------------------------- extracting DDS results ---------------------------------------------
  #' If results is run without specifying contrast or name (name is for continuous variables), \
  #' it will return the comparison of the last level of the last variable in the design formula \
  #' over the first level of this variable. (log2(lastGroup/firstGroup) => 
  #' **firstGroup (deonominator) listed in the levels of a factor, will be considered CONTROL** 
  #' To set contrast manually, set contrast=c("condition","Treated","Control/Negative/Normal")
  dsqRES = results(dsqRAW, cooksCutoff = F, contrast = reffLvL) ## refLVL = c("condition","Treated","Control")
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
  
         all<- list(dsqRES_tax, agg_at,  DA_type,   amp_obj,   fcts, corr,      dds,   dsqRAW, dsqRES)
  names(all) <- c("DA_resMTX", "agg_at", "DA_type", "amp_obj", "fcts", "corr", "ddsObj", "dsqRAW", "dsqRES")
  
  return(all)
}






#' ------------------------------------------------------------------------------------------------------------------
# --------------- my Maaslin funciton --------------------------------------------

myMaAsLin2 <- function(amp_obj, fcts, rand = NULL, maasDIR= "maasOUT", minPrev=0.1, minAbu=0.0, 
                       reffLvL=NULL,  agg_at = "OTU", model = "LM", transf="NONE", norm="CLR"){
  ## transf="LOG", or norm="CLR", not both
  ## for reffLvL=paste0("HPVstatus,HPVneg")
  # reffLvl=levels(amp_obj$metadata[[fcts[[1]]]])
  
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
  set.seed(123)
  # suppressMessages( suppressWarnings(
  maasOUTa = Maaslin2::Maaslin2(input_data = counts, input_metadata = meta, output = maasDIR,
                                fixed_effects = fcts, ## order of vars does NOT matter, all vars are adj for ea/o equal
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






#' ----------------------------------------------------------------------------------------------------------------
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
  
  if(model %in% c("LM", "CPLM","NEGBIN","ZINB") ){norm=NULL; transf=NULL } # wrapper will auto do individual combinations 
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
  
  AICs <- data.frame(AIC=unlist(lapply(maasOUT[5]$fits, function(x) tryCatch(AIC_forcplm(x), error=function(e) NA))))
  maasOUT <- merge(maasOUT$results, round(AICs, 2), by.x="feature", by.y=0)
  
  # print(maasOUT %>% head() )
  
  ##' --- output formatting
  set.seed(123)
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


#' ------------------------------------------------------------------------------------------------
# ----------------------------------my ALDEX function (fancy) --------------------------------------
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
  
  meta   <- amp_obj$metadata
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
        library(BiocParallel); library(parallel)
        if(is.null(threads)) {threads=detectCores()-4 }
        register(MulticoreParam(threads))
        set.seed(123)
  ###' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###' Aldex pre-calculations
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



#' ------------------------------------------------------------------------------------------------
# ------------------ my LMER DA function ------------------
##' linear mixed effects model, set as differential abundance fuction 
##' for longitudinal data
##' ensure reffLevel is the first of the fixed effect variable (its also the string_to_remove)
##' e.g. formula is ABX_AnyDuringNICU + PMA + (1|SubjectID) : must include fixed effect & longitudinal effect & random effect
##' Grab_term is to choose which term from the stat you want the effect size & p-value for
##' 
##' 
myLmerDA <- function(ampObj, formula, minPREV= 10, method = "anova", mainVar, #timeVar, subjVar,
                     agg_at = "OTU", norm = "CLR", pseudocount = 0.01 ){
  
  if(! agg_at %in% c("OTU")){
    myScripts("dataClean.R")
    ampObj <- amp_agg_abund(ampObj, agg_at = agg_at)
  }
  
  myScripts("dataCorrections.R")
  normCounts <- normMyData(ampObj$abund, norm_method = norm , pseudocount = pseudocount)[[1]]
  # normCounts[1:10, 1:10]
  normMatrix <- normCounts %>% t() %>% as.data.frame()
  # normMatrix[1:10, 1:10]
  
  ## for all taxa/aro##s
  library(MuMIn)
  
  ##' method 1 (ano): using anova R2 as effect size & p-value stats from the limMDL:
  ##' - conditional effect size (for the whole model), but lost directionality info
  ##' method 2 (sum): using summary estimate & p-value stats from limMDL:
  ##' - this will have effect size directional (+/- effect size relative to reference level)
  ##' - this will have multiple effect sizes and multiple p-values- one for each pair of refLvL-vs-fctLvl for each taxa
  ##' - here I am selecting to represent the stats for the pair with biggest effsiz (& smaller p-value) for each taxon
  
  lmr_all <- as.data.frame(matrix(ncol = 3, nrow = length(colnames(normMatrix)) ), row.names = colnames(normMatrix) )
  colnames(lmr_all) <- c( "EffectPair", "EffSiZ", "pval")
  
  
  ## fixed_effect name to remove from effPair
  # string_to_remove<- ifelse(is.null(string_to_remove), strsplit(formula, split = "\\+| ")[[1]][[1]], string_to_remove) 
  
  model_formula <- as.formula(paste("normMatrix[[i]]", formula))  
  
  if (method == "anova"){
    ##' method 1 (anova): using anova R2 as effect size & p-value stats from the limMDL:
    ##' - conditional effect size (for the whole model), but lost directionality info
    print("calculating linear mixed effects & grabbing anova(model) stats for each feature")
    for (i in colnames(normMatrix) ){
      limMDL <- lmerTest::lmer(model_formula , data = ampObj$metadata)
      stat <- anova(limMDL, by = "margin")
      d<- dim(stat)
      avoP <- stat[d[[1]], d[[2]] ]  ## the last-row & last-col is p-value of interest Def is the fixed_effect 
      # eta_sqr <-  stat[["Sum Sq"]][d[[1]]]  / sum(stat[["Sum Sq"]]) ## variance proportion captured by that 1 variable. We dont want this for DA
      avoZ <- r.squaredGLMM(limMDL)[[2]] ## variance explained by full model (all effects & conditions) == R2conditional == effectsize (non-directional)
      lmr_all[i, ] <- list("global", avoZ, avoP) %>% as.data.frame()
    }
  }else if (method == "summary") {
    ##' method 2 (sum)mary: using the estimate & pvalue from the summary of limMDL:
    ##' - this will have effect size be directional (+/- effect size relative to reference level)
    ##' - this will have multiple effect sizes and multiple p-values- one for each pair of refLvL-vs-fctLvl for each taxa
    ##' - here I am selecting to represent the stats for the pair with biggest effsiz (& smaller p-value) for each taxon
    print("calculating linear mixed effects & grabbing the summary(model) stats for each feature")    
    for (i in colnames(normMatrix) ){
      limMDL <- lmerTest::lmer(model_formula , data = ampObj$metadata)
      stat <- summary(limMDL)$coefficients
      d <- dim(stat); l <- d[[1]]- (length(levels(ampObj$metadata[[mainVar]])) -2 ) ## start-row to grab from
      sum1 <- stat[ c(l:d[[1]]) , c(1,5)] %>% as.data.frame() %>% rownames_to_column("EffPair") %>% 
        mutate(EffPair = str_remove(EffPair, ".*\\:"))  %>% mutate(EffPair= str_remove(EffPair , mainVar))
      ## grabbing the last N rows, representing the interactive or additive values for the mainVAR
      sum1 <- sum1 %>% slice_max(order_by = desc((Estimate) ), n = 1) ### grabbing the pair with biggest effect size
      lmr_all[i, ] <- sum1 
    }
  }else{ print("method has to be one of 'anova' or 'summary'"); stop()   }
  
  print("example stat output:")
  print(stat) ## prints the last stat (anova or summary) select proper grab_term if needed
  
  # lmr_all %>% heil()
  lmr_all %<>% mutate(padj = p.adjust(pval, method = "BH") ) ## adjsting the pvals using Benjamini-Hochberg correction
  lmr_all %<>% mutate(prev = rowSums(ampObj$abund >0))          ## adding a PREV column, for downstream filtering
  
  ## adding taxRanks
  tu <- ampObj$tax %>% select(any_of(c(agg_at, "Class", "Phylum"))) %>% 
    rename_with(~c("TAX", "Class" , "Phylm")) %>% mutate(Taxon = TAX, .before = TAX) ; #colnames(tu)
  
  lmr_alltu <-  bind_cols(lmr_all, tu)                    ## adding in the taxonomy info
  lmr_alltu %<>% filter(prev > minPREV)                   ## the downstream filtering
  
  
  ## Gathering steps
  lmerOUT <- list()
  lmerOUT$DA_resMTX <- lmr_alltu
  lmerOUT$amp_obj   <- ampObj
  lmerOUT$formula   <- formula
  lmerOUT$agg_at    <- agg_at
  lmerOUT$DA_type   <- "lmerDA"
  
  return(lmerOUT)
}

# ============ a little tutorial on LMER() =================================
# #### A little tutorial on lmer()
# - anova(lmer()) vs summary(lmer())
# ANOVA is for testing lmer() with continuous or categorical variables (like timepoint or estrogen levels)
# SUMMARY is for testing lmer() models with continuous only variables (), but gives more info
# 
# - the random effect term
# Usually + (1|SubjID). Required by lmer()
# TimePoint is NOT a random term. its longitudinal term
# 
#  - How to find longitudinal effect
#       do:
#       anova(lmer(PCoN + TimePoint + (1|PID)) 
#      
#  - how to know which model better fits the data:
#        lets say we have:
#        m1<- lmer(PCoN ~ mainVar+TP + (1|PID) ) [additive model]
#        m2<- lmer(PCoN ~ mainVar*TP + (1|PID) ) [interactive model]
#          
#     **method 1**: check interaction
#        do **anova(m2)**: 
#     - if interaction is SIG => m2 the more complex model is significantly better. 
#     - if interaction NS, m1 (simpler, additive model) is preferred due to parsimony (simplicity)
#        
#     **method 2**: check model fit 
#      do **AIC(m1, m2)** and **BIC(m1,m2)**
#     - the lower AIC or BIC value represents the better model fit
#          
#     Note: it is normal to have sig anova and higher AIC/BIC model fit. means that the interaction does\
#     explain sig more variation but at the cost of added complexity to the model fit.
#     If the interaction is theoretically important (e.g. biological reasoning), keep it despite the higher AIC.
#     Example: the different between "correction for the time effect" (additive model) vs "as time progresses" effect (interactive model)
#      
#     
#  - how to know **if variable is a confounder** and needs to be corrected for 
#          setup:
#            m1 <- lmer(PCoN ~ mainVAR + (1|PID) ) [no confounder]
#          m2 <- lmer(PCoN ~ mainVAR + pConf + (1|PID)) [added potential confounder]
#          do **anova(m1, m2)**
#          - if anova(m1,m2) is SIG => adding the pConf is necessary (pConf is a confounder)
#          - if anova(m1,m2) is NS, adding the pConf is not needed
#          or do **compare the F-values** of m1 and m2 
#          - if `m2$F` >10% off from `m1$F`, then pConf is confounder (explains more variance).
  
#' ===================================================================================


#' ===================================================================================
# ---------- my Byasian additive mixed effects model ------------------------------------------
##' original by Katie, saved in NICU project, untitled.R
##' works with DA plotter
myBamDA<- function(ampObj, FctVar, TimeVar, SubjVar = "SubjectID", selTaxRanks = taxRanks, minPREV = 10, agg_at = "OTU" ) {
  source(paste0("~/Documents/MyApps/RScripts/myRscripts/dataClean.R"))
  source(paste0("~/Documents/MyApps/RScripts/myRscripts/otherFunctions.R"))
  
  if(agg_at!= "OTU"){
    print(paste("aggregating TAXa to", agg_at, "level"))
    ampObj <- amp_agg_abund(ampObj, agg_at = agg_at)
    # print(amp_obj$tax %>% head(12))
  }
  
  phyObj <- amp2phy(ampObj)
  
  taxon <- psmelt(phyObj) %>% mutate(
    SubjectID = factor(get(SubjVar)), 
    TimeVar = get(TimeVar),
    Var = factor(get(FctVar))) %>% 
    group_by(OTU)  ## groups set here is how do() knows what to apply to
  taxon_gam <- taxon %>% do(broom::tidy(bam(log(Abundance+1) ~ Var + s(TimeVar, by=Var) + s(SubjectID, bs="re"), data=.), parametric=T)) %>% 
    filter(grepl("Var", term)) %>% ungroup() %>% mutate(p.adj = p.adjust(as.numeric(p.value), "fdr")) %>% select(!term)  
  # taxon_gam
  tax_info <- ampObj$tax %>% select(selTaxRanks) %>% data.frame() %>% rownames_to_column()
  tab_gam <- taxon_gam %>%  arrange(OTU, p.value) %>% distinct(OTU, .keep_all = T) %>%   left_join(tax_info)
  tab_gam %<>% dplyr::rename(TAX = Species, Taxon = Genus, Phylm = Phylum) %>% column_to_rownames("rowname")
  tab_gam %<>% mutate(prev =  rowSums(ampObj$abund > 0), .after = "p.adj") %>% arrange(prev) %>% filter(prev > minPREV)
  # tab_gam  %>% heil()
  
  gamdaOUT <- list()
  gamdaOUT$DA_resMTX <- tab_gam 
  gamdaOUT$DA_type <- "bamDA"
  gamdaOUT$amp_obj <- ampObj %>% amp_subset_taxa(tab_gam$OTU)
  gamdaOUT$fcts <- FctVar
  gamdaOUT$rand <- SubjVar
  gamdaOUT$corr <- TimeVar
  gamdaOUT$agg_at <- agg_at
  
  return(gamdaOUT)
}


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

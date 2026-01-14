#' ------------------ betaDiv plot (ordination) ---------------------------------------
# manditory intput: amp object, and 1 factor of interest
# optional inputs: 
#   a 2nd factor of interest,
#   a type_plot other than PCA, and
#   a distance matrix other than Bray c("bray", "euclidean", "canbrera")
beta_plot <- function(amp , shape_var, color_var= shape_var, type= "PCoA", dist = "bray",
                      sampLab=NULL, abund_cutOff = 0.01, pt_size=4, labSize=3,  
                      colVector= myColors$discrete, pchVector=myPchR$all, alpha = 0.8, 
                      efit_fct=NULL, efit_num=NULL, efit_Pval = 0.05, x_axis = 1, y_axis = 2,
                      efit_scale = 1, efit_color = "darkgray", reverse_x = F, reverse_y = F,
                      title=paste(type,"(",dist ,")", "for", color_var)){
  
  pca_ucli <- amp_ordinate(amp, filter_species = abund_cutOff, type = type, detailed_output = TRUE,
                           distmeasure = dist, transform = "none", x_axis =x_axis, y_axis = y_axis, 
                           envfit_factor= efit_fct, envfit_numeric = efit_num, envfit_signif_level = efit_Pval, 
                           envfit_numeric_arrows_scale = efit_scale, envfit_arrowcolor = efit_color, envfit_textcolor = efit_color,
                           sample_point_size = pt_size, sample_label_size = labSize, opacity = alpha, 
                           sample_shape_by = shape_var,  sample_color_by = color_var, sample_label_by=sampLab) 
  
    if(reverse_x == T){revX <- scale_x_reverse(labels = \(x) x * -1)}else{ revX <- NULL }
    if(reverse_y == T){revY <- scale_y_reverse(labels = \(y) y * -1)}else{ revY <- NULL }
  
  if (is.factor(amp$metadata[[color_var]]) || is.character( amp$metadata[[color_var]] )  ) {
    print("color_var is recognized as categorical")      
    pca_ucli$plot <- pca_ucli$plot + scale_shape_manual(values=pchVector) + labs(title=title) +
      scale_color_manual(values=colVector) + 
      revX + revY + coord_fixed(ratio = NULL) # + theme(aspect.ratio = 1) to make it square
  }else{
    print("Color_factor is recognized as numeric")
    pca_ucli$plot <- pca_ucli$plot + scale_shape_manual(values=pchVector) + labs(title=title) +
      scale_color_gradient(low = colVector[[1]], high = colVector[[2]]) + 
      revX + revY + coord_fixed(ratio = NULL) #theme(aspect.ratio = 1)
  }#else if(is.finite.POSIXlt(amp$metadata[[color_var]][[1]]) && all(is.finite(as.numeric(amp$metadata[[color_var]]))) ){
  #       print("color_var is recognized as Date")
  #       pca_ucli$plot= pca_ucli$plot + scale_shape_manual(values=pchVector) + labs(title=title) +
  #         scale_color_datetime(low = colVector[[1]], high = colVector[[2]]) 
  
  # return(pca_ucli$plot) ## $plot does not contain $plot$data anymore :(
  ## use + scale_x_reverse() or + sscale_y_reverse(labels = \(y) y*-1) to reverse coordinates
  pca_ucli$data <- pca_ucli$dsites %>% select(1: PCo10)
  return(pca_ucli)
}

# amp_ordinate(btchOUTspp, type = "PCoA", distmeasure = "bray", 
# sample_color_by = "Nreads", sample_shape_by = "Batch") + 
#   labs(title = "PCoA with Bray dist for number of reads in CORR data") +
#   scale_color_gradient(low = "violet", high = "violetred4")

#' ----------------- betaPlot with TaxLoadings labels -----------------
#' You will need to input the ampObj & pcoObj outputted from beta_plot
#' Generally, leave arrows as F, as I dont like them
#' the taxLabs_scale adjusts the distance of the taxloading labels
#' leave ntax as low value, as you will always get ~3* more taxa
#' method can be spearman, pearson or kendall

beta_arrowLoads <- function(ampObj, pcoaObj, ntax=3, method = "spearman", fillVector = myColors$pastel, 
                            sampCol="SampleID", colorRank ="Order", taxLabs_scale = 0.5, sel_by = "top",
                            fill_title = "Tax contributions",
                            arrows = F, taxLabs_alpha = 0.85, taxLabs_size = 3, selRank = "Species"){
 
      pco <-pcoaObj$data %>% select(all_of(c(sampCol, "PCo1", "PCo2"))) %>% column_to_rownames(sampCol)
      otutab <- ampObj$abund 
      tax <- ampObj$tax 
      arrow_width <- ifelse(arrows == T, 0.8, 0) 
      
      merged <- merge(pco, t(otutab), by=0)
      otus <- tax$OTU
      sav <- list()
      
      for(i in 1:length(otus)) {
        taxon <- otus[i]
        pco1 <- cor(merged$PCo1, merged[,otus[i]], method=method)*taxLabs_scale ## scaling factor for distancing the PCos
        pco2 <- cor(merged$PCo2, merged[,otus[i]], method=method)*taxLabs_scale ## scaling factor for distancing the PCos
        sav[[i]] <- cbind(taxon, pco1, pco2)
      }
      
      cor.results <- data.frame(do.call(rbind, sav)) %>%
        mutate(across(2:3, as.numeric)) %>% arrange(pco1, pco2)
      final <- merge(cor.results, tax, by.x="taxon", by.y="OTU") %>% 
        arrange(desc(pco1), desc(pco2));          # dim(final)
      
      if (sel_by == "slice"){
            abs1 <- final %>% slice_max(order_by = abs(pco1), n = ntax) ## top_n selects only + values => using slice_max
            abs2 <- final %>% slice_max(order_by = abs(pco2), n = ntax) ## top_n(ntax, pco2)
            abs <- rbind(abs1, abs2) %>% distinct(selRank, .keep_all = T)
            # abs %>% arrange(desc(abs(pco1)))
      }else if (sel_by == "top"){
            ## it is more equal representation of negative & positive correlation values, but there is more taxa
            absp1 <- final %>% top_n (ntax, pco1);                              ##top_frac(0.005, pco2);  subset(abs(pco2)>0.85);
            absn1 <- final %>% top_n (round(ntax/2, 0), desc(pco1))             ## using half of the negative correlation taxa
            absp2 <- final %>% top_n (ntax, pco2); #                            ##top_frac(0.005, pco2);     subset(abs(pco2)>0.85);
            absn2 <- final %>% top_n (round(ntax/2, 0), desc(pco2))             ## using half of the negative correlation taxa
            abs <- rbind(absp1, absn1, absp2, absn2) ;  # dim(abs)
            abs %<>% distinct(get(selRank), .keep_all = T);  # dim(abs)
            # abs %>% arrange(desc(abs(pco1)))
      }else{ message("choose taxa selection method: top or slice"); return()    }    
      
      
      betaTAX <- pcoaObj$plot + 
        geom_segment(data = abs, mapping = aes(x = 0, y = 0, xend = pco1, yend = pco2), inherit.aes = F, ## geom_segment is needed for proper legend
                     arrow = arrow(length = unit(0.2, "cm")), linewidth = arrow_width, alpha = 0.35, color = "gray")  +
        geom_label_repel(data = abs, aes(x = pco1, y = pco2, fill = get(colorRank), label = get(selRank)),  
                         inherit.aes = F, color = "black", size = taxLabs_size, box.padding = 0.3, alpha = taxLabs_alpha) + 
        scale_fill_manual(values = fillVector) + labs(fill = fill_title) + #: \n", colorRank, "/", selRank)) + 
        guides(color = guide_legend(order = 3), shape = guide_legend(order = 2),
               fill = guide_legend(order = 1, override.aes = list(label = "  ", color = NA) ) )
      
      outObj <- list(betaTAX, abs); names(outObj) <- c("plot", "tax_loadings")  
  
  return(outObj)
}

#' ------------------ lmer on PCos with random effect -------------------------------
#' This function is to test PCo variance, longitudinally for a factor, adjusting for random SubjID. So
#' instead of :
#' anova(lmerTest::lmer(betENV$data$PCo1 ~ SeasonDay * HeadImpact_Catg+ (1|SubjectID), data= betENV$data)) ## interaction is NS
#' anova(lmerTest::lmer(betENV$data$PCo1 ~ SeasonDay + HeadImpact_Catg+ (1|SubjectID), data= betENV$data)) ## additions also NS
#' anova(lmerTest::lmer(betENV$data$PCo2 ~ SeasonDay * HeadImpact_Catg+ (1|SubjectID), data= betENV$data)) ## interaction is NS
#' anova(lmerTest::lmer(betENV$data$PCo2 ~ SeasonDay + HeadImpact_Catg+ (1|SubjectID), data= betENV$data)) ## additions also NS
#' do the lmerPCo_model()
#' 
#' Ensure your formula starts with **"~"** , then the longi & continuous & all correction variables,
#' with interactive or additive terms **"+"** or __"*"__: the model can be either interactive or additive, 
#' The **1-before-LAST VAR in the model should be the mainVAR of interest** (for collecting pVals & explained proportion variance)
#' Finally, the last argument is the **random effect (1|SubjectID)**
#' template lm_formula = "~ Time + corrVAR + mainVAR + (1|SubjectID)"
#' example  lm_formaula = "~ Day + GameSeason + (1|SubjectID)"
#' essentially **put var of interest just before random effect**, so proper pVals & explained variance proportion is captured

lmerPCo_model <- function(lm_formula, dat, npcs = 10){
        pc <- paste0("PCo", 1:npcs)
        res_tab <- setNames( data.frame(matrix(ncol = 2, nrow = npcs)) , c("prop", "pval") )
        row.names(res_tab) <- pc 
        lmr_aov_list <- list()
        
        for (i in pc) {
          model <- paste0(i, lm_formula)
          lmr_aov <- anova(lmerTest::lmer(model, data=dat))# this anova is type III already, doesnt need by = "margins"
          # grep values from interaction: always last line of output
          d <- dim(lmr_aov)
          res_tab[[i, "pval"]] <- lmr_aov[d[[1]], d[[2]] ] ## last row, last column
                ## eta_squred == standard effect size (0-1) similar to R2
                ss_eft  <- lmr_aov[["Sum Sq"]][d[[1]]]
                ss_ttl  <- sum(lmr_aov[["Sum Sq"]])
          eta_sqr <-  ss_eft / ss_ttl
          res_tab[[i, "prop"]] <-  round(eta_sqr, 4)
          lmr_aov_list[[i]] <- lmr_aov
        }
  res_tab %<>%  add_significance("pval", "sig")
  print(res_tab)
  outObj <- setNames(list(res_tab, lmr_aov_list), c("res_tab", "lmer_aov_all"))
  return(outObj)
}



#' ------------------ betaDiv lmer & lmer Plot for PCos (longitudinal)  ------------------------
#' This function presents tax loadings on PCo1 and PCo2 (despite significance), along with
#' uses longitudinal trend changes in tax loadings across groups on any of 1-10 PCos. 
#' It plots the tax trend changes PCo1, PCo2 and PCo-significance based on the lm_formula. 
#'
#' the lmerPCo_model(lm) calcs the longitudinal tax trend changes for PCo1-10 [def], based on model 
#' for the anova(lmer(model)) result, the last column & last line is the pvals & propVar of interest
#' this plt2 could still use some work

beta_lmerPCo <- function(ampObj, time_var="Time", rand_var = "SubjectID",  #catg_var = "SampleType",
                         lm_formula = "~ Time * SampleType + (1|SubjectID)",
                         colrVector = myColors$discrete, legend_rows = 2,  
                         dist="bray", pco_sig = 0.05, pco_prop = 0.5, npcs = 10,
                         lfit = "glm", facet_nrows = NULL, facet_ncols = NULL, 
                         se_fill = "gray", se_alpha = 0.2, se_level = 0.60  ) { 
  
  dist <- vegdist(t(ampObj$abund), method = dist)
  ordi <- cmdscale(dist, k = npcs)
  colnames(ordi) <- paste0("PCo", 1:npcs)
  df <- merge(ampObj$metadata, ordi, by=0)
  # ampObj$tax$SppASV <- paste0(ampObj$tax$Species, " (", ampObj$tax$OTU, ")")
  
  results_table <- lmerPCo_model(lm_formula = lm_formula, dat = df, npcs = npcs)$res_tab
  
  # here loop will produce loadings for PCo1, PCo2, and all the significant ones with coef >0.5
  sig_pcos = unique(c(rownames(results_table)[c(1,2) ],                            # grabbing PCo1 and PCo2 even if NS
                      rownames(results_table)[
                      results_table[,"pval"] <= pco_sig & results_table[, "prop"] >= pco_prop] ) ) # grabbing significant other PCos
  
  lin <- all_cors <- list()
  for(ii in sig_pcos) {
    lin[[ii]] <- ggplot(df, aes(x=!!sym(time_var), y=!!sym(ii), color =!!sym(rand_var) ) )  + 
      geom_point( alpha = 0.5) + 
      stat_smooth(method=lfit, se = T, level = se_level, fill = se_fill, alpha = se_alpha) + 
      theme_bw() + ylab(ii) + guides(color = guide_legend(nrow = legend_rows )) +
      scale_color_manual(values = colrVector) + 
      ## hjust > 1 pushes text **left from the right edge**; vjust > 1 pushes text **down from the top edge*
      annotate("text", x=Inf, y=Inf, hjust=1.025, vjust=1.15, 
               label=paste0("p.Val = ", sprintf("%.2e", results_table[ii,2])))
    
    all_cors[[ii]] <- t(cor(df[,ii], t(ampObj$abund), method="s")) %>% na.omit() %>% data.frame() %>% 
      merge(ampObj$tax, by=0) %>% dplyr::rename(corr = '.') %>%
      column_to_rownames("Row.names") %>% 
      slice_max(order_by = abs(corr), n = 10) 
    ## draft below
    #                     sel_cors <- slice_max(all_cors, abs(corr), n=10) ## 
    # box[[ii]] <- ggplot(sel_cors, aes(x=corr, y=reorder(SppASV, corr), fill=Class)) + ## not sure I like SppASV
    #   geom_col() +  theme_bw() + scale_fill_manual(values=colrVector) + ylab("Taxa") + xlab(ii)
  }
  linPlot <- ggarrange(plotlist = lin, legend = "bottom", ncol = facet_ncols, nrow = facet_nrows, common.legend = T)
  # boxPlot <- ggarrange(plotlist = box,  legend = "bottom", ncol = facet_ncols, common.legend = T)
  
  outObj <- list(all_cors, linPlot); names(outObj) <- c("taxa_loadings", "plot")
  return(outObj)
}


# beta_div_pcs(amp_fcl, lm_formula="~Week*Genotype + (1|SubjectID)") 

#' ---------------- boxplot of taxonomic loadings from PCoA ----------
##' it's basically a large-ish loop, where you correlate each taxon with PC1
##' and then plot the largest/most significant correlations."
##' plots a boxplot of taxonomic loadings on each PCo

beta_boxLoads <- function(ampObj, pcoaObj, pcoN="PCo1", ntax=20, method = "spearman", 
                          sampCol="SampName", fillRank="Class" , 
                          y.angle = -90, text.sz = 12, y.position="left"){
    pco <-pcoaObj$data %>% select(all_of(c(sampCol, pcoN))) %>% column_to_rownames(sampCol)
    
    otutab <- ampObj$abund 
    tax <- ampObj$tax 
    
    merged <- merge(pco, t(otutab), by=0)
    otus <- tax$OTU
    sav <- list()
    
    for(i in 1:length(otus)) {
      sav1 <- otus[i]
      sav2 <- cor(merged[[pcoN]], merged[,otus[i]], method=method)
      sav[[i]] <- cbind(sav1, sav2)
    }
    
    cor.results <- data.frame(do.call(rbind, sav)) %>%
      mutate(sav2= as.numeric(as.character(sav2))) %>% 
      arrange(sav2)
    final <- merge(cor.results, tax, by.x="sav1", by.y="OTU") %>% 
      arrange(desc(sav2));  
    absp <- final %>% top_n (ntax, sav2); # ##top_frac(0.005, sav2);  subset(abs(sav2)>0.85); 
    absn <- final %>% top_n (round(ntax/2, 0), desc(sav2)) ## using half of the negative correlation taxa
    abs <- rbind(absp, absn)
    # final[c(-3,-8)] %>% tail(n=8)
    # cor.results %>% View()
    
    ## plotting
    library(rlang) ## allows me to call the fillRank content as string
    lp <- ggplot(abs, aes(x=sav2, y = fct_reorder(Species, sav2), fill=!!sym(fillRank) )) + 
      geom_bar(stat = "identity") + theme_bw()+ 
      scale_fill_manual(values=rev(myColors$discrete), name=fillRank ) + 
      xlab(paste(pcoN, "loadings")) + ylab("Feature") + 
      theme(axis.text.y = element_text(angle=0, size=text.sz))

    if(pcoN == "PCo2"){
    lp <- lp + scale_y_discrete(position = y.position) + coord_flip() +
        theme(axis.text.x = element_text(angle=y.angle, size =text.sz, hjust = 0)) 
    }
    
    outObj <- list(lp, abs); names(outObj) <- c("plot", "matrix")
    return(outObj)
    # return(lp)
}

##' usage
# beta_loadings(amp, p)


#' ------------------ heatmap plot ---------------------------------------

heat_plot<- function(amp_obj, facet_fct, group_fct = facet_fct, ntax=35, normalize = T, 
                     tax_agg = "Species", tax_add = "Class",  leg.position = "bottom", measure = "mean", 
                     title =paste("\n Top", ntax,"taxa (at", tax_agg, "level)") , showRemn = T, value_size = 3,
                     x.angle = 0, min_abund = 0.01, colscale = "sqrt", #def: "log10", but sqrt better
                     subtitle = paste0("Grouped by ", group_fct, ", faceted by ", facet_fct))
  {
  amp_obj$abund %<>% as.data.frame()
  amp_htmp <- amp_heatmap(amp_obj, group_by = group_fct, facet_by = facet_fct, tax_add = tax_add, 
                          tax_aggregate = tax_agg, tax_show=ntax, plot_values = T, measure=measure,
                          min_abundance = min_abund, plot_colorscale = colscale,  #c("log10", "sqrt"),
                          plot_values_size = value_size, normalise=normalize, showRemainingTaxa = showRemn)  + 
    labs(title= title, subtitle = subtitle) + #scale_y_discrete(labels \(y) gsub(y.break, "\n", y)) +
    theme(axis.text.x = element_text(angle = x.angle, size=12, vjust = 0.5, hjust=0.5), #element_blank() , # 
          axis.text.y = element_text(size=10), legend.position=leg.position,
          plot.title = element_text(size=14), plot.subtitle = element_text(size=12) );
  return((amp_htmp))
}
# heat_plot + 
# scale_fill_gradient2(low = "steelblue", mid="lemonchiffon", high = "darkorange", 
# midpoint = 10, transform = "sqrt") ## def:"log10" ,but sometimes sqrt is better


#' ------------------------ profile bar plot -----------------------------------------------------------------
#' 
#' bar_plot <- function(amp_obj, facet_fct, group_fct, ...){
# transform abund to TSS with normMyData
# ggplot(indata, aes(x=sample, y=`%unambiguousReads` + `%ambiguousReads`, fill=`#name`)) + 
# geom_bar(stat="identity", color="black") + 
#   scale_fill_brewer(palette="Set3") +
#   theme(axis.text.x=element_text(angle=-90, vjust=0, hjust=0)) +
#   facet_wrap(~LysisTube+ExtractionKit+BeadBeater, scales="free_x", nrow = 1) +
#   scale_y_continuous(breaks=seq(0,100,10)) +
#   geom_text(data=subset(indata, !Description %in% ""), aes(x=sample, y=101, label="*"), size=8, check_overlap = T)
#' 
#' }
#' 
#' 





#' ------------------------ alpha Indices calculations -------------------------------

alpha_indices <- function(abund){
  library(vegan)
  tmp=NA
  #Observed species richness and maximum sp. richness:
  tmp$minindiv <- min(colSums(abund))  #alpha$minindiv
  tmp$obs_sprch_max=max(vegan::specnumber(abund, MARGIN = 2))  
  tmp$maxsamp= which(specnumber(abund, MARGIN = 2)== 
                       max(specnumber(abund, MARGIN = 2)))
  #expected species richness and maximum expected:
  tmp$exp_sprch=rarefy(abund, tmp$minindiv, se=F, MARGIN = 2);  #tmp$exp_sprch
  tmp$exp_sprch_max=max(rarefy(abund, tmp$minindiv, se=F, MARGIN = 2));  #tmp$exp_sprch_max
  # cat("---- Most Species rich sample is ", names(tmp$maxsamp), " with ", tmp$obs_sprch_max, " observed species. 
  #     Expected maximum number of species in a sample:", round(tmp$exp_sprch_max,0))
  # 
  #Shannon, Simpson and Inverse Simpson,  Pielou's Evenness, dominance:
  alpha=list()
  alpha$reads=colSums(abund);  # alpha$reads #Number reads per sample
  alpha$SppRichness=vegan::specnumber(abund, MARGIN = 2);  # alpha$SppRichness #number of species
  alpha$Shannon=vegan::diversity(as.matrix(abund), "shannon", MARGIN = 2);  # head(alpha$Shannon)
  alpha$InvSimpson=vegan::diversity(abund, "inv", MARGIN = 2);   #head(alpha$InvSimpson)
  alpha$Evenness=alpha$Shannon/log(specnumber(abund, MARGIN = 2)); #head(alpha$Evenness) #this is Pielou's Evenness index
  alpha$chao1=apply(abund, 2, fossil::chao1); #head(alpha$chao1)
  
  
  #Amending the metadata table with the produced indices (based on countsR data)
  alpha$table=cbind(alpha$SppRichness, alpha$Evenness, alpha$reads, 
                    alpha$Shannon, alpha$InvSimpson, alpha$chao1)
  colnames(alpha$table)=c("SppRichness", "Evenness", "NumbReads", 
                          "Shannon", "InvSimpson", "chao1")
  
  # check normalcy of the indices
  als <- alpha$table %>% as.data.frame
  try({ for (a in colnames(als)){
    cat("Shapiro-Wilk normality distribution test for", a, "\t")
    shp <- shapiro.test(als[[a]])
    print( shp$p.value   )
  } })
  cat("SIG p-Values == abnormal data distribution => non-parametric tests")
  
  return(round(alpha$table, 2))
}

#' ------------------------ AlphaDiv Boxplot -------------------------------
#' Ensure the factors of interest are in ordered levels of interest (manually)
#'  wilcox.test == Mann-Whitney U test (non-paramentric)
alpha_plot<- function(metadata, fct1, fct2 = fct1, sigplot=T, legend.pos = "right", x.angle=45,
                      points=F , boxes =T, yScale=NULL, ylim=1, wrap_rows=1, medLine=F, 
                      sig_test = "wilcox.test", ## or sig_test = "t.test", lines = F, #linesVAR = "SubjectID", #<- problem still
                      which_alphas= "main", y.lab="Alpha Index Value", titlelab=NULL,  alphaColr=0.5,
                      colVector = myColors$discrete, pchVector=rep(myPchR$all[1], 10) ){
  sst=list()
  sst$meta=metadata
  sst$meta$group1= factor(sst$meta[[fct1]] )
  sst$meta$group2= factor(sst$meta[[fct2]] )
  # sst$meta$group3= factor(sst$meta[[linesVAR]]) ## if SubjectID doesnt exist, this will be a problem
  alphas0=c("SppRichness", "Evenness", "chao1" , "Shannon", "InvSimpson","NumbReads")
  
  alpha_options = c("all", "main", "old", "div")
  if     ( which_alphas[1] == "main"){alphas=alphas0[c(1,3,4,5)] } # alphas w Chao1
  else if( which_alphas[1] == "old"){ alphas=alphas0[c(1,2,3,4)] } # alphas w Evennness
  else if( which_alphas[1] == "all"){ alphas=alphas0             } # alphas all
  else if( which_alphas[1] == "div"){ alphas=alphas0[c(4,5)]     } # alphas Shannon & InvSimp ONLY
  else { alphas = which_alphas }  
  
   
  sst$meta=sst$meta[, c("group1", "group2",alphas)] # "group3", 
  sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2") )# "group3", 
  print(sst$melt %>% head())
  if(boxes==T){gm_box=geom_boxplot(aes(x=group1, y=as.numeric(as.character(value), stat="identity", position="dodge"), 
                                            fill=group2), outlier.shape=NA, alpha=alphaColr)  }else{gm_box=NULL}

  if(points==T){gm_point=geom_point(aes(group=group2, shape=group2), show.legend=T, ##, color=NULL, 
                                    position=position_jitterdodge(jitter.width = 0.5) ) }else{gm_point=NULL} 
  # if(lines==T) {gm_line=geom_line(aes(group = group3), color = "gray80",  alpha = 0.7,
  #                                   position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75) ) }else(gm_line=NULL)
   
  if(medLine==T){gm_mean=stat_summary(aes(x=group1, y=value, group=group2), fun = "median",
                                      geom="crossbar",  width = 0.5, fatten = 2, color="red",
                                      show.legend = F, position = position_dodge(width = 0.8))}else{gm_mean=NULL}
  
  if(!is.null(yScale)){yScale=scale_y_sqrt(limits=c(0,ylim))}

  p1=ggplot(sst$melt, aes(group1, value)) + facet_wrap(~variable, scales="free", nrow=wrap_rows) + 
      gm_box + gm_point + gm_mean +  yScale  + #gm_line + 
      scale_color_manual(values = colVector ) + scale_fill_manual(values=colVector )  +  
      scale_shape_manual(values=pchVector) +  
      labs(y=y.lab, x=fct1, fill=fct2 , title=titlelab, shape=fct2) + theme_bw() + 
      theme(text=element_text(size=18),  
              axis.text.x = element_text(angle=x.angle, size=14, hjust=0.5, vjust=0.5), #hjust=1, 
              legend.position = legend.pos,  legend.title= element_text(size=12), #legend.position = "none"
              legend.text=element_text(size=10) ) ;
     
   if(sigplot == T){
        a<-length(unique(sst$melt$group1) )
        if (a > 1){
        p1= p1 + geom_signif(comparisons=combn(a, 2, simplify = F), 
                             test=sig_test,# t.test is paramentric, wilcox.test nonparametric
                             map_signif_level = c("***"= 0.001, "**"=0.01, "*"= 0.05, "." = 0.065, " "=1),
                             step_increase = 0.05, color="gray60", tip_length = 0.01, vjust = 0.55);  
                              # map_signif_level = F for actual p-values or
                  }
                      }
  return(p1)
}
##' for paired points do example:
##' Where BMTpair has uniqueIDs for unpaired (e.g. row names), and redundant IDs for paired (e.g. PID)
# geom_point(aes(shape = SampleTypeA, fill = SampleTypeA, color = SampleTypeA, group = BMTpair), size = 2, 
#            position = position_dodge(width = 0.5)) +
#   geom_line(aes(group = BMTpair), color = cols[[2]], alpha = 0.6,   position = position_dodge(width = 0.5))


##' instructions on how to modify plot to look as above but add sig.Labels in every pair of sub-groups (group2) in X
##' Essentially, each sub-group is in its own facets, and I have formatted the theme to look like above
##' Look into Audrey's Salmonella project for more alterations

# a <- alpha_plot(submeta, fct1 = "DPI",fct2 = i, sigplot = F, points = T, y.lab = "Percent host reads", 
#                 which_alphas = "hostContent", legend.pos = "top"); plot(a)
# 
# p<- ggplot(a$data, aes(group1, value)) + facet_wrap(~group1, scales = "free_x", nrow = 1, strip.position = "bottom") +
#   geom_point(aes(x=group2, color=group2, y= value), show.legend=F, position=position_jitterdodge(jitter.width = 0.5)) +
#   geom_boxplot(aes(x=group2, fill=group2, y= value), outlier.shape=NA, alpha=0.5, show.legend = T, width=1) + theme_bw() +
#   labs(title="Host reads content", x=bquote(bold("DPI")), y="Fraction of sample reads", fill="SampleType") + scale_color_discrete(guide="none") +
#   theme(panel.border = element_blank(), strip.background = element_blank(), panel.spacing.x = unit(0.25, "lines"), 
#         axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
#         strip.text.x = element_text(size=12, angle = 45), text=element_text(size=12),  
#         plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"), legend.position = "top") ; plot(p)
# 
# b <- length(unique(a$data$group2))
# p + geom_signif(aes(x=group2, y=value), comparisons=combn(b,2,simplify = F), test = "wilcox.test")


#' ------------------------ AlphaDiv Linear Model Stats -------------------------------
#' Ensure the factors of interest are in ordered levels of interest (manually)
#' Not a boxplot. This presents a continuous X-axis (e.g. longitudinal TP), 
#' plots samples as points colored by categorical factor
#' Calculates the change (slope; regression LM) of Alpha Indices over time,
#' between the categorical factor variables. 
#' Samples need to be repeated measurements of same IndivID
#' 
#' 
alpha_lmplot <- function(metadata, longi_fct , other_fct = longi_fct, subjID="SubjectID", adjst_fct = NULL,
                         x.angle = 45, 
                         yLab = "value", calc_stats =T, yVARs = c("Shannon","InvSimpson"), wrap_rows=1){
      sst=list()
      sst$meta=metadata
      sst$meta$longi= sst$meta[[longi_fct]]            # use the longitudinal factor
      sst$meta$group= sst$meta[[other_fct]]
      sst$meta$sbjID= sst$meta[[subjID]]
      # if(is.null(adjst_fct)){sst$meta$adjst <- "none" } else { sst$meta$adjst <- sst$meta[[adjst_fct]] } 
      # adjfcts <- ifelse(is.null(adjst_fct), "none", paste(adjst_fct, collapse = "+")  ) 

      sst$meta  %<>% select(any_of(c( "group", adjst_fct ,"longi", "sbjID", yVARs)) )
      sst$melt <- melt(as.data.table(sst$meta), id.var=c("longi", adjst_fct,  "group","sbjID") ) %>%
        mutate(across(c(value, longi), ~as.numeric(as.character(.)))) %>% 
        mutate(across(c(group, adjst_fct, sbjID), ~ as.factor(.)))
      
      # isadj <- ifelse(is.null(adjst_fct), "", "+ adjst")
      isadj <- ifelse(is.null(adjst_fct), "", paste0("+", paste(adjst_fct, collapse = "+")) ) 
      fml <- as.formula(paste("value ~longi * group" , isadj, "+ (1|sbjID)" ) )
      if (calc_stats == T){
        # stat1 <- sst$melt %>%  filter(variable %in% yVARs[[1]]) %>% 
        #   lmerTest::lmer(fml, data=. ) %>% summary()
        #   stat1 <- stat1$coefficients
        stat1 <- sst$melt %>%  filter(variable %in% yVARs[[1]]) %>% 
          lmerTest::lmer(fml, data=. ) %>% anova()
                           ## data= . → returns the entire inputted dataframe
        print(stat1)
        # stat1 <- stat1[[dim(stat1)[1], dim(stat1)[2]]]
        
        # stats <- sst$melt %>%  group_by(variable) %>%
        #   summarize(pVal = summary(lmerTest::lmer(fml))$coef[4,5] )
        stats <- sst$melt %>%  group_by(variable) %>%    ## pick(everything()) → returns all columns in the current group slice
          summarize(pVal = anova(lmerTest::lmer(fml, data = pick(everything()) ))[ "longi:group", "Pr(>F)" ] )
      }else{
        stat1=paste0("not calculated")
        stats <- data.frame(x="", pVal = Inf )
      }
      
      
      p1=ggplot(sst$melt, aes(x=as.numeric(as.character(longi)), y=value)) + 
        facet_wrap(~variable, scales="free", nrow=wrap_rows) + 
        geom_point (aes(color= group), position=position_jitterdodge(jitter.width = 0.3), alpha = 0.3, size = 0.7)  +
        geom_smooth(aes(color= group), method = "lm", alpha= 0.1) +
        theme_bw() +  labs(y=yLab, x=longi_fct, color=other_fct ) +
        theme(text=element_text(size=16), axis.text.x = element_text(angle=x.angle, hjust=1, size=14),
              legend.position = "right",  legend.title= element_text(size=12), #legend.position = "none"
              legend.text=element_text(size=10)) +
        geom_text(aes(x=Inf, y=Inf, label=paste0("p = ", formatC(pVal, 2, format = "e"))),
                  data=stats, vjust=2, hjust=1.1, inherit.aes = F) 
      #' hjust = 1.1,  # >1 moves text left; <1 moves right
      #' vjust = 1.5   # >1 moves text down; <1 moves up
      out <- list(p1, stat1)
      names(out) <- c("plot_lm", "stats_lm")
      return(out)
}
# example use: alpha_lmplot(meta_fcl, "Week", "Genotype")



#' ------------------------ AlphaDiv Statistical tests -------------------------------
#' ---------------------------- alpha LM intersect -----------------------------------
#' tests for significant interaction btween alphaDiv variables 
alpha_lmX <- function(metadata, fct1, fct2 = fct1, subjID="SubjectID"){
  sst=list()
  sst$meta=metadata
  sst$meta$group1= as.numeric(as.character( factor(sst$meta[[fct1]] ) )) # use the longitudinal factor
  sst$meta$group2= factor(sst$meta[[fct2]] ) 
  sst$meta$subjID= sst$meta[[subjID]]
  sst$meta=sst$meta[, c("group1", "group2", "SppRichness","chao1" , 
                        "Shannon", "InvSimpson", "subjID")]
  sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2","subjID") )
  
  stats <- sst$melt %>%  group_by(variable) %>%
    summarize(pVal = summary(lmerTest::lmer(as.numeric(value) ~ group1*group2 + 
                                              (1|subjID)))$coef[4,5] )
  return(stats)
}
#' alpha_lmplot(meta_fcl, "Week", "Genotype")
#' output of stats looks like:
#' Fixed effects:
#                                            Estimate   Std.Err   df       t value     Pr(>|t|)    
# (Intercept)                                67.6979    11.5437 152.0000   5.864     0.0000271 ***
# Week                                       4.2880     3.2870 152.0000   1.305        0.194    
# GenotypeKO                                 -2.7372    14.7154 152.0000  -0.186        0.853    
# Week:GenotypeKO                            -0.6579     4.1902 152.0000  -0.157        0.875   
#' so we'd always want the 5th column (the p-value) & the last row (the interaction) $coef[x,5]

#' ----------------------- alpha LM independent -----------------------------
#' alternative to KW and ANOVA, but when lognitudinal change within SubjID is involved

alpha_lm <- function(metadata, fct1, fct2 = fct1, subjID="SubjectID"){
  sst=list()
  sst$meta=metadata
  sst$meta$group1= factor(as.numeric(as.character(sst$meta[[fct1]] ) )) # use the longitudinal factor, as numeric
  sst$meta$group2= factor(sst$meta[[fct2]] )
  sst$meta$subjID= sst$meta[[subjID]]
  sst$meta=sst$meta[, c("group1", "group2", "SppRichness","chao1" , "Shannon", "InvSimpson", "subjID")]
  sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2","subjID") )

  stats <- sst$melt %>%  group_by(variable) %>%
  summary(lm(fct1 ~ fct2, data=sst$meta))
  stats <- sst$melt %>%  group_by(variable) %>%
    summarize(pVal = summary(lm(fct1 ~ fct2, data=sst$meta)) )
  return(stats)
}


#' ------------------------ PARAMETRIC ANOVA tests -------------------------------
#'  PARAMETRIC ANOVA tests
#'  alpha_anova function can perform multiple types of ANOVA tests, based on a list of factors and a model:
#'    - one-way anova:: evaluates 1 factor with multiple levels
#'    - two-way anova:: evaluates additive un-related fcts) with +, see: https://tinyl.io/87cX
#'    - two-way anova:: evaluates interactive fcts with * , see https://tinyl.io/87yF)
#'    - nested-factor anova:: evaluates nested-by-study-design fcts, Info at: https://tinyl.io/87cZ
#'      Note: nests cannot be time points. Longitudinal analysis is ranked => considered 'repeated measure' not 'nest'
#'  I could also add pw-multiple comparisons.Info at: https://tinyl.io/87cM


#' I am switching to library(car)::Anova(lm(), type = I, II or III) function, where
#' Type I: tests terms sequentially & accumulatively (marginal) ("Sequential Sums of Squares"), but expects no interaction 
#'        => it tests the effects of variables in a fixed order
#'        base::aov() Applicable only for balanced designs, where !is.list(replications(model, meta))  == T 
#'        car::Anova(lm(fct), type = I) can also handle unbalanced desings.
#'        Note: base::aov(XYZ + B) == car::Anova(lm(B+XYZ), type = I) are the same for factor B only! 
#'        Values for XYZ is unreliable with aov() as they are marginal of their priors
#' Type II is Anova(lm(), type= "II") - tests terms sequentially & accumulatively like in Type I, but 
#'        can handle interactions among preceeding terms. Type II cannot handle ranked terms (only categorical, where
#'        all categories are equal ranks)
#' Type III is Anova(lm(), type="III") - tests the effects of predictor variables in a way that is independent of
#         the order in which they were entered into the model. The test assumes that the effects of each variable 
#         can only be understood in the context of the other variables in the model ("Partial sum of squares"), 
#         and interactions among variables are expected. TypeIII is useful for complex interactions
#'        Note: an error about "aliased coefficients in the model",  means that some coefficients in model are 
#'        perfectly correlated
#' info at: https://tinyl.io/8Fn5


#' ------------------------  KW and AovI AlphaDiv test -------------------------------
#' This performs one-way anova (Parametric) with KW test (non-par)
#' This function is perfect for inside loops testing 1 factor at a time
anova1 <- function(meta, fct){
  sst=list()
  sst$meta=meta
  variables = c("SppRichness","chao1" , "Shannon", "InvSimpson")
  sst$meta=sst$meta[, c(fct, variables)]
  sst$melt=melt(as.data.table(sst$meta), id.var=c(fct) )
  
  # cat("==== Kruskal-Wallis test on AlphaDiv contrasts for factor", i, "====== \n")
  fml=as.formula(paste0("value ~", fct))
  
  #' KW test 
  suppressWarnings({
  sst$kw_models=sst$melt %>% group_by(variable) %>% do(Model=kruskal.test(fml, data=.) )
  names(sst$kw_models$Model)=paste(sst$kw_models$variable)#, "~", fct) #"KW test::", 
  })
  
  krusk_raw <- sst$kw_models$Model
  
  krusk_pvals <- Reduce(rbind, lapply(sst$kw_models$Model, '[', "p.value" ))
  row.names(krusk_pvals)<- paste(sst$kw_models$variable)#, "~", fct) , "~", fct)
  krusk_pvals <- krusk_pvals %>% data.frame() %>% rename_at(1, ~"KW.pVal") %>%
    mutate(KW.pVal = round(as.numeric(KW.pVal), 8) ) %>% rownames_to_column("Comparison")
  
  #' ANOVA type I test
  suppressWarnings({
  sst$nma=sst$melt %>% group_by(variable) %>% do(Model=aov(fml , data=.) )
  sst$ano_models= lapply(sst$nma$Model, summary) #summary?
  names(sst$ano_models)=paste(sst$nma$variable)#, "~", fct) , "~", fct) #"ANOVA test::", 
  })
  ano1_raw = sst$ano_models
  
  ano1_pvals = as.data.frame(unlist(sst$ano_models)) %>% rename_at(1, ~"ANOVA1.pVal") %>%  
    rownames_to_column("Comparison") %>% filter(stringr::str_detect(Comparison, '.Pr..F..')) %>% 
    mutate_all(~sub('.Pr..F..', "", .)) %>%   mutate(ANOVA1.pVal = round(as.numeric(ANOVA1.pVal), 8) ) %>% 
    na.omit(ANOVA1.pVal)
  
  merged_pvals <- full_join(ano1_pvals, krusk_pvals, by = "Comparison")
  
  anv0=list()
  anv0=list(krusk_raw, ano1_raw, merged_pvals)
  names(anv0) <- c("kw_raw", "ano1_raw", "merged_pVals")
  
  return(anv0)
  
}

#' ------------- type III ANOVA where all fcts are corrected for e/a ------------------
#' BEWARE: easy over-correction can occur if a lot of factors are string together
#' This can be used on Shannon vs Continuous variable
#' glm & lm take only fixed effect, no random effect fixes

# library(car)

anova3 <- function(meta , model = "~ fct1 + fct2", type = "III"){
 
  fcts <- all.vars(as.formula(model))
  sst=list()
  sst$meta=meta
  sst$meta=sst$meta[, c(fcts ,"SppRichness","Evenness" , "chao1" , "Shannon", "InvSimpson")]
  sst$melt=melt(as.data.table(sst$meta), id.var=c(fcts) )
  
  fml1=as.formula(paste("value ", model))
  suppressWarnings({
  sst$nma=sst$melt %>% group_by(variable) %>% do(Model=car::Anova(lm(fml1 , data=.), type = type) )
  })
  sst$models= sst$nma$Model 
  names(sst$models)=paste(sst$nma$variable)#, "~", model) #"ANOVA test::", 
  anv3 = list()
  anv3$raw_output = sst$models
  
  
  anv3$p.values = as.data.frame(sst$models) %>%  t() %>% as.data.frame() %>% rownames_to_column("ANOVA3.pVal") %>%
    select(ANOVA3.pVal, fcts) %>% filter(stringr::str_detect(ANOVA3.pVal, '.Pr..F.')) %>% 
    mutate(across("ANOVA3.pVal", str_replace, '.Pr..F.', " ~ adj.Factor")) %>% column_to_rownames("ANOVA3.pVal") %>% 
    setNames(paste0("adj.", names(.), " ")) %>% rownames_to_column("ANOVA3.pVal")
    # 
    # as.data.frame(sst$models) %>%  rownames_to_column("ANOVA3.Comparisons") %>%
    # select(c(ANOVA3.Comparisons, contains('.Pr..F.')) ) %>% 
    # rename_with(~gsub('.Pr..F.', "~ ", .)  ) %>%
    # filter(!stringr::str_detect(ANOVA3.Comparisons, c("(Intercept)" ))) %>% na.omit()
  
  
  names(anv3) <- c("raw_output", "p.values")
  return(anv3)
}

 
####' rank transformation / non-parametric / ANOVA-like test
####' for when a non-parametric version of ANOVA3 is needed, where
####' the main variable needs to be adjusted for a confounder
anolm3 <- function(meta , model = "~ fct1 + fct2",
                   alphas = c("SppRichness","Evenness", "NumbReads" , "chao1" , "Shannon", "InvSimpson")){
  sst=list()
  sst$meta=meta
  fcts <- all.vars(as.formula(model)) ## extracts the variables from a formula
  sst$meta=sst$meta[, c(fcts , alphas)]
  sst$melt=melt(as.data.table(sst$meta), id.var=c(fcts) )
  
  fml1=as.formula(paste("rank(value)", model))
  suppressWarnings({
    sst$nma=sst$melt %>% group_by(variable) %>% do(Model=anova(lm(fml1, data=.)) )
  })
    sst$models= sst$nma$Model
    names(sst$models) <- paste(sst$nma$variable)
    lmv3 = list()
    lmv3$raw_output = sst$models
    lmv3$p.values   = as.data.frame(sst$models) %>%  t() %>% as.data.frame() %>% rownames_to_column("anolm.pVal") %>%
      select(anolm.pVal, fcts) %>% filter(stringr::str_detect(anolm.pVal, '.Pr..F.')) %>% 
      mutate(across("anolm.pVal", str_replace, '.Pr..F.', " ~ ")) %>% column_to_rownames("anolm.pVal") %>% 
      setNames(paste0("adj.", names(.), " ")) %>% rownames_to_column("anolm.pVal")
  
  names(lmv3) <- c("raw_output", "p.values")
  return(lmv3)
}






#' -------------- Effect of Confounding interaction evaluation -----------------------
#' Function performs logistic regression with lm(X ~ v + j), to find potential confounding factors (j) 
#' to variables of interest (v) as it interacts w numeric variable AlphaDiv (X). The way to assess if j has
#' effect on v's contribution to X, is to assess the change in "Estimate" coefficient for any of the levels of V.
#' Model allows for mainVars & putvConf vas can be numerical, but AlphaIDXs to have 
#' normal distribution (NS shapiro test), since lm() is technically a parametric test
#' if abnormal distribution on AlphaDiv, a Spearman correlation test cor.test(method = "spearman")
#' is needed (for numeric V and J)
#' 
#' It is a bunch of nested loops, where for each alpha index, function calculates the change between 
#' Estimate coefficients of the mainVar of   
#' [lm(AlphaIdx(X) ~ mainVar(V))$Estimate] and [lm(AlphaIdx(X)~mainVar(V) + putvConf(J))$Estimate],
#' Since multi-level mainVar would have Estimate coefficients for each group **relative to refLVL**, 
#' The loops then calculates if ANY of the changes within all the groups are **affected** by the putvConf. 
#' Function gives out a kable showing which putative confounders need to be corrected with in ANOVA3,
#' when calculating anova(Shannon~mainVar)
#' **Important: factor(levels=c("refLevel", "OtherLvL1", "OtherLvL2", ...))**
#' 
#'  Side notes: on what not to do and why
#'  1) Model like anova(glm(V~J, family="binomial"), test="Chisq") is for only for testing binomial (categorical) 
#'    variables for interactions with each other. If any of the variables are numeric, one needs to use gaussian 
#'    distributions with model: glm(V~j, family = "gaussian"), which is the same as lm(V~J). This gives
#'    a best model fit assessment (lowest AIC value). However, if you want to also add covariables like:
#'    lm(V~j + f), the AIC value needs to be replaced by BIC value, which adds penalty to the added variables 
#'    for better assessment of model fit.
#'  2) For assessing confounding interactions on AlphaDiv effect the anova() wrapper around glm/lm is not to be used.
#'  3) Also anova (, test="Chisq") is only useful for glm(, family = "binomial") [measuring numeric var inter]. 
#'    The default anova test is the F-test:  anova(..., test="F") 

alphaConf <- function(amp , vars , conf, alphas=c("Shannon", "InvSimpson")){
  
  tab <- finVal <- jab <- zab <- list()
  tab.title = "Confounding Effects Evals [lm(alpha~ X + y )]"
  bch=alphas ## those are non-categorical, numeric/continuous response variables 
  amp$metadata <- amp$metadata %>% mutate_at(bch, as.numeric)
  
  for (x in bch){
    for (v in vars){
      fmlv <- as.formula(paste(x, "~", v)) ## calc main effect
      lmv <- summary(lm(fmlv, data = amp$metadata))$coefficients %>% as.data.frame();
      nv=grep(v, rownames(lmv), value =T);
      lmv %<>% filter(rownames(.) %in% nv) %>% select(Estimate) ## grab coeff for MainVar EffSz bf interaction
      
      for (j in conf) {
        fmlc=as.formula(paste(x , "~", v, "+", j)) ## calc influence of confVar
        lmc <- summary(lm(fmlc, data = amp$metadata))$coefficients %>% as.data.frame(); 
        nc=grep(v, rownames(lmc), value =T) 
        lmc %<>% filter(rownames(.) %in% nc) %>% select(Estimate) ## grab coeff of MainVar after interaction
        jab <- cbind(lmv, lmc) 
        colnames(jab) <- c("var", "conf")
        change <- abs(jab$conf-jab$var)/abs(jab$var)*100 #calc change of effSz
        jab$change <- ifelse(change >=10, "yes", "no") #add it as column
        finVal <- ifelse(any(jab$change == "yes"), "affected", "unaffected") ## assess if any changes
        zab <- cbind(alphaVar = x, mainVar = v, confVar = j, change=finVal)
        tab <- rbind(tab, zab);
        rm(finVal, zab, change, nc, lmc, fmlc)
      }
    }
    rm(fmlv, lmv, nv)
  }
  tab %<>% as.data.frame() %>% kable(caption = tab.title) %>% kable_styling(c("striped", "hover")) %>%
    scroll_box(width="80%" , height = "100%")
  return(tab)
}


#' anova_anova (old) testing 1 term at a time
#' ------------- type III ANOVA where all fcts are corrected for e/a ------------------
#' BEWARE: easy overcorrection can occur if a lot of factors are string together
# function is same as above but the output format is a bit different
alpha_anova <- function(meta , fcts, model = paste(fcts, collapse = "+" ), type = "III"){
  sst=list()
  sst$meta=meta
  sst$meta=sst$meta[, c(fcts, "Shannon", "InvSimpson")] #"SppRichness","chao1" ,
  sst$melt=melt(as.data.table(sst$meta), id.var=c(fcts) )
  
    fml1=as.formula(paste("value ~", model))
    sst$nma=sst$melt %>% group_by(variable) %>% do(Model=car::Anova(lm(fml1 , data=.), type = type) )
    sst$models= sst$nma$Model #lapply(sst$nma$Model, anova)
    names(sst$models)=paste(sst$nma$variable, "~", model) #"ANOVA test::", 
    anv1 = list()
    anv1$p.values = as.data.frame(unlist(sst$models)) %>% rename_at(1, ~"ANOVA.pVal") %>%  
      rownames_to_column("Comparison") %>% filter(stringr::str_detect(Comparison, '.Pr..F.2')) %>%  #Pr..F.1 = (Intercept)
      mutate_all(~sub('.Pr..F.2', "", .)) %>% mutate(ANOVA.pVal = as.numeric(ANOVA.pVal))
    anv1$raw_output = sst$models
    return(anv1)
}



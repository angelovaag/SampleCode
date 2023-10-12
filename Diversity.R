#' ------------------ betaDiv plot (ordination) ---------------------------------------
# manditory intput: amp object, and 1 factor of interest
# optional inputs: 
#   a 2nd factor of interest,
#   a type_plot other than PCA, and
#   a distance matrix other than Bray c("bray", "euclidean", "canbrera")
beta_plot <- function(amp , fct1, fct2=fct1, type= paste("PCoA"), dist = paste("bray"),
                      sampLab="SampleID", abund_cutOff = 0.01, colVector= myColors$discrete, pchVector=myPchR$all,
                      title=paste(type, "with",dist, "distances for", fct1, "and", fct2)){
  pca_ucli <- amp_ordinate(amp, filter_species = abund_cutOff, type = type, detailed_output = TRUE,
                           distmeasure = dist, transform = "none", #sample_colorframe = fct2,
                           sample_point_size = 4, sample_label_size = 3, #na.rm=T, 
                           sample_shape_by = fct1,  sample_color_by = fct2, sample_label_by=sampLab) 
  pca_ucli$plot= pca_ucli$plot +  scale_color_manual(values=colVector) + # stat_ellipse(level = 0.2) +
    scale_shape_manual(values=pchVector)  +  labs(title=title ) 
  return(pca_ucli$plot)
}

#' ------------------ betaDiv lm and lm plot  ------------------------
#' This function presents tax loadings on PCo1 and PCo2 (despite significance), along with
#' uses longitudinal trend changes in tax loadings across groups on any of 1-10 PCos. 
#' It plots the tax trend changes PCo1, PCo2 and PCo-significant. 
#' 
#' 
#' the rep_model(lm) calcs the longitudinal tax trend changes for PCo1-10, based on model 
#' for the LM model result, the last line is one of interest
#' This formula is not yet optimized at the aes() level

beta_div_pcs <- function(ampObj, mydist="bray", ntaxa=15,
                         fct="SampleType", time_var="Week", tax_fill="Class", #check aes()
                         lm_formula = "~ Week * SampleType + (1|SubjectID)" ) {
  
  dist <- vegdist(t(ampObj$abund), method = mydist)
  ordi <- cmdscale(dist, k = 10)
  colnames(ordi) <- paste0("Axis.", 1:10)
  df <- merge(ampObj$metadata, ordi, by=0)
  ampObj$tax$SppASV=paste0(ampObj$tax$Species, " (", ampObj$tax$OTU, ")")
  
  #sub-function for LM of PC interactions
  rep_model <- function(pc, dat , lm_fml= lm_formula) { 
    model <- paste0(pc, lm_fml)
    lm_res <- summary(lmerTest::lmer(model, data=dat) )
    # greps values from the factor/interaction of interest: the last line of the output
    lm_fct <- length(lm_res$coeff[, 5] )
    res_pval <- round(lm_res$coef[lm_fct, 5], 3)
    res_coef <- round(lm_res$coef[lm_fct, 1], 3)
    return(list(coef=res_coef, pval=res_pval)) }
  
  pc_res <- lapply(paste0("Axis.", 1:10), rep_model, dat=df)
  names(pc_res) <- paste0("Axis.", 1:10)
  results_table <- do.call(rbind.data.frame, pc_res) %>% data.frame()
  print(results_table)
  axis_plots <- list()
  
  sel_taxa <- c(unique(ampObj$tax[,tax_fill]))
  sel_colors <- rev(col)[1:length(sel_taxa)]
  names(sel_colors) <- sel_taxa
  
  # here loop will produce loadings for axis1,2, and all the significant ones
  which_pcos = unique(c(rownames(results_table)[c(1,2) ],                       #grabbing PCo1 and PCo2 for the loadings
                      rownames(results_table)[results_table[,"pval"] < 0.05]) ) #grabbing significant other PCos
  for(ii in which_pcos) {
    all_cors <- t(cor(df[,ii], t(ampObj$abund), method="s")) %>% 
      na.omit() %>% data.frame() %>% merge(ampObj$tax, by=0) %>% dplyr::rename(corr = '.')
    
    plt1 <- ggplot(df, aes(x=Week, y=df[,ii], color=Genotype)) +
      geom_point() +  geom_smooth(method="lm") + theme_bw() + ylab(ii) +
      annotate("text", x=Inf, y=Inf, hjust=1, vjust=1,
               label=paste0("p.Val = ", round(results_table[ii,2],3)))
    
    plt2 <- ggplot(slice_max(all_cors, abs(corr), n=ntaxa),
                   aes(x=corr, y=reorder(SppASV, corr), fill=Class)) + #reorder(Speices , corr) not sure I like SppASV
      geom_col() +  theme_bw() +
      scale_fill_manual(values=sel_colors) + ylab("Taxa") + xlab(ii)
    # plot(plt1)
    print(ggarrange(plt1, plt2, nrow=1, widths=c(0.5,1)))
  }
}

# beta_div_pcs(amp_fcl, lm_formula="~Week*Genotype + (1|SubjectID)") 



#' ------------------ heatmap plot ---------------------------------------

heat_plot<- function(amp_obj, fct1, fct2 = fct1, ntax=35, normalize = T, tax_agg = "Species", tax_add = "Class",  leg.position = "bottom",
                     title =paste("\n TAX profiles of top",ntax,"taxa (at", tax_agg, "level)") , x.angle = 0, min_abund = 0.01, 
                     colscale = "log10", subtitle = paste("Presenting groups in", fct1, " subgrouped by", fct2, "levels"))
  {
  amp_obj$abund %<>% as.data.frame()
  amp_htmp <- amp_heatmap(amp_obj, group_by = fct2, facet_by = fct1, tax_add =tax_add, 
                          tax_aggregate = tax_agg, tax_show=ntax, plot_values = T, 
                          min_abundance = min_abund, plot_colorscale = colscale,  #c("log10", "sqrt"),
                          plot_values_size = 3, normalise=normalize, showRemainingTaxa = T)  + 
    labs(title= title, subtitle = subtitle) + 
    theme(axis.text.x = element_text(angle = x.angle, size=12, vjust = 0.5, hjust=0.5), #element_blank() , # 
          axis.text.y = element_text(size=10), legend.position=leg.position,
          plot.title = element_text(size=14), plot.subtitle = element_text(size=12) );
  return((amp_htmp))
}
# heat_plot <- htmp_plot 

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
  return(round(alpha$table, 2))
}

#' ------------------------ AlphaDiv Boxplot -------------------------------
#' Ensure the factors of interest are in ordered levels of interest (manually)
#'  wilcox.test == Mann-Whitney U test (non-paramentric)
alpha_plot<- function(metadata, fct1, fct2 = fct1, sigplot=T, legend.pos = "right", x.angle=45, points=F ,yScale=NULL, ylim=1,
                      which_alphas= "main", y.lab="Alpha Diversity", titlelab=NULL,  alphaColr=1, colVector = myColors$discrete){
  sst=list()
  sst$meta=metadata
  sst$meta$group1= factor(sst$meta[[fct1]] )
  sst$meta$group2= factor(sst$meta[[fct2]] )
  alphas0=c("SppRichness", "Evenness", "chao1" , "Shannon", "InvSimpson","NumbReads")
  
  alpha_options = c("all", "main", "old", "div")
  if( which_alphas == "main"){alphas=alphas0[c(1,3,4,5)] } # alphas w Chao1
  if( which_alphas == "old"){ alphas=alphas0[c(1,2,3,4)] } # alphas w Evennness
  if( which_alphas == "all"){ alphas=alphas0             } # alphas all
  if( which_alphas == "div"){ alphas=alphas0[c(4,5)]     } # alphas Shannon & InvSimp ONLY
  if(!which_alphas %in% alpha_options){ alphas = which_alphas }  
 
  sst$meta=sst$meta[, c("group1", "group2", alphas)]
  sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2") )
  
  # if(points==T){gm_point=geom_point(aes(color=group2), position=position_dodge2(width = 0.3), alpha=alphaColr) }else{gm_point=NULL} 
  if(points==T){gm_point=geom_point(position=position_dodge2(width = 0.2)) }else{gm_point=NULL} 
  if(!is.null(yScale)){yScale=scale_y_sqrt(limits=c(0,ylim))}

  p1=ggplot(sst$melt, aes(group1, value)) + facet_wrap(~variable, scales="free") + 
    geom_boxplot(aes(y=as.numeric(as.character(value), stat="identity", position="dodge"), 
                     x=group1, fill=group2), outlier.shape=NA, alpha=alphaColr)  +      gm_point +  yScale + 
    labs(y=y.lab, x=fct1, fill=fct2 , title=titlelab) +  scale_fill_manual(values=colVector )  + 
    theme_bw() + theme(text=element_text(size=18),  axis.text.x = element_text(angle=x.angle, size=14, hjust=0.5, vjust=0.5), #hjust=1, 
          legend.position = legend.pos,  legend.title= element_text(size=12), #legend.position = "none"
          legend.text=element_text(size=10) ) ;
     
   if(sigplot == T){
        a<-length(unique(sst$melt$group1) )
        if (a > 1){
        p1= p1 + geom_signif(comparisons=combn(a, 2, simplify = F), 
                             test="wilcox.test",# map_signif_level = F for actual p-values or
                             map_signif_level = c("***"= 0.001, "**"=0.01, "*"= 0.05, "." = 0.065, " "=1),
                             step_increase = 0.05, color="gray60", tip_length = 0.01, vjust = 0.55);  
                  }
                      }
  return(p1)
}

#' ------------------------ AlphaDiv Linear Model Stats -------------------------------
#' Ensure the factors of interest are in ordered levels of interest (manually)
#' Not a boxplot. This presents a continuous X-axis (e.g. longitudinal TP), 
#' plots samples as points colored by categorical factor
#' Calculates the change (slope; regression LM) of Alpha Indices over time,
#' between the categorical factor variables. 
#' Samples need to be repeated measurements of same IndivID

alpha_lmplot <- function(metadata, fct1, fct2 = fct1, subjID="SubjectID"){
  sst=list()
  sst$meta=metadata
  sst$meta$group1= factor(sst$meta[[fct1]] ) # use the longitudinal factor
  sst$meta$group2= factor(sst$meta[[fct2]] ) 
  sst$meta$subjID= sst$meta[[subjID]]
  sst$meta=sst$meta[, c("group1", "group2", "SppRichness","chao1" , "Shannon", "InvSimpson", "subjID")]
  sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2","subjID") )
  
  stats <- sst$melt %>%  group_by(variable) %>%
    summarize(pVal = summary(lmerTest::lmer(as.numeric(value) ~ 
                                              as.numeric(as.character(group1))*group2 + (1|subjID)))$coef[4,5] )
  
  p1=ggplot(sst$melt, aes(x=as.numeric(as.character(group1)), y=value)) + 
    facet_wrap(~variable, scales="free") + 
    geom_smooth(aes(color= group2, fill = group2), method = "lm", alpha= 0.1) +
    geom_point(aes(color=group2), position=position_jitterdodge(jitter.width = 0.3))  +
    theme_bw() +  labs(y="Alpha Diversity", x=fct1, fill=fct2 ) +
    theme(text=element_text(size=16), axis.text.x = element_text(angle=45, hjust=1, size=14),
          legend.position = "right",  legend.title= element_text(size=12), #legend.position = "none"
          legend.text=element_text(size=10)) +
    geom_text(aes(x=Inf, y=Inf, label=paste0("pVal=", round(pVal,3))),
              data=stats, vjust=1, hjust=1, inherit.aes = F)
  return(plot(p1))
}
# example use: alpha_lmplot(meta_fcl, "Week", "Genotype")



#' ------------------------ AlphaDiv Statistical tests -------------------------------
#' ---------------------------- alpha LM --------------------------------------------
alpha_lm <- function(metadata, fct1, fct2 = fct1, subjID="SubjectID"){
  sst=list()
  sst$meta=metadata
  sst$meta$group1= as.numeric(as.character( factor(sst$meta[[fct1]] ) )) # use the longitudinal factor
  sst$meta$group2= factor(sst$meta[[fct2]] ) 
  sst$meta$subjID= sst$meta[[subjID]]
  sst$meta=sst$meta[, c("group1", "group2", "SppRichness","chao1" , "Shannon", "InvSimpson", "subjID")]
  sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2","subjID") )
  
  stats <- sst$melt %>%  group_by(variable) %>%
    summarize(pVal = summary(lmerTest::lmer(as.numeric(value) ~ group1*group2 + (1|subjID)))$coef[4,5] )
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
#' so we'd always want the 5th column (p-value) & the last row (the interaction) $coef[x,5]




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
  sst$kw_models=sst$melt %>% group_by(variable) %>% do(Model=kruskal.test(fml, data=.) )
  names(sst$kw_models$Model)=paste(sst$kw_models$variable, "~", fct) #"KW test::", 
  
  krusk_raw <- sst$kw_models$Model
  
  krusk_pvals <- Reduce(rbind, lapply(sst$kw_models$Model, '[', "p.value" ))
  row.names(krusk_pvals)<- paste(sst$kw_models$variable, "~", fct)
  krusk_pvals <- krusk_pvals %>% data.frame() %>% rename_at(1, ~"KW.pVal") %>%
    mutate(KW.pVal = round(as.numeric(KW.pVal), 8) ) %>% rownames_to_column("Comparison")
  
  #' ANOVA type I test
  sst$nma=sst$melt %>% group_by(variable) %>% do(Model=aov(fml , data=.) )
  sst$ano_models= lapply(sst$nma$Model, summary) #summary?
  names(sst$ano_models)=paste(sst$nma$variable, "~", fct) #"ANOVA test::", 
  
  ano1_raw = sst$ano_models
  
  ano1_pvals = as.data.frame(unlist(sst$ano_models)) %>% rename_at(1, ~"ANOVA.pVal") %>%  
    rownames_to_column("Comparison") %>% filter(stringr::str_detect(Comparison, '.Pr..F..')) %>% 
    mutate_all(~sub('.Pr..F..', "", .)) %>%   mutate(ANOVA.pVal = round(as.numeric(ANOVA.pVal), 8) ) %>% 
    na.omit(ANOVA.pVal)
  
  merged_pvals <- full_join(ano1_pvals, krusk_pvals, by = "Comparison")
  
  anv0=list()
  anv0=list(krusk_raw, ano1_raw, merged_pvals)
  names(anv0) <- c("kw_raw", "ano1_raw", "merged_pVals")
  
  return(anv0)
  
}

#' ------------- type III ANOVA where all fcts are corrected for e/a ------------------
#' BEWARE: easy over-correction can occur if a lot of factors are string together
# library(car)
anova3 <- function(meta , fcts, model = paste(fcts, collapse = "+" ), type = "III"){
  sst=list()
  sst$meta=meta
  sst$meta=sst$meta[, c(fcts, "SppRichness","chao1" , "Shannon", "InvSimpson")]
  sst$melt=melt(as.data.table(sst$meta), id.var=c(fcts) )
  
  fml1=as.formula(paste("value ~", model))
  sst$nma=sst$melt %>% group_by(variable) %>% do(Model=car::Anova(lm(fml1 , data=.), type = type) )
  sst$models= sst$nma$Model #lapply(sst$nma$Model, anova)
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


#' anova_anova (old) testing 1 term at a time
#' ------------- type III ANOVA where all fcts are corrected for e/a ------------------
#' BEWARE: easy overcorrection can occur if a lot of factors are string together
# function is same as above but the output format is a bit different
alpha_anova <- function(meta , fcts, model = paste(fcts, collapse = "+" ), type = "III"){
  sst=list()
  sst$meta=meta
  sst$meta=sst$meta[, c(fcts, "SppRichness","chao1" , "Shannon", "InvSimpson")]
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



#' -------------- Effect Size and Significance Evaluation ------------------------------------
#' Function performs anosim and adonis (PERMANOVA) on bch factors, using amp object
#' Both COUNTS and METADATA is needed from amp object
#' 
#' PERMANOVA tests whether distance differ between groups.
#' ANOSIM tests whether distances between groups are greater than within groups. =>
#' If var w/i group <   var b/w group, ANOSIM:: SIG (bw/wi >   1, SIG) \
#' If var w/i group >/= var b/w group, ANOSIM:: NS  (bw/wi </= 1, NS)
#' 
#' 
#' to use function: 
#' source("~/Documents/MyApps/RScripts/myRscripts/effectEval.R")
#' amp <- data$DZs$amp
#' bch <- c("SampleType", "Disease"), etc
#' effectEval(bch, amp)
#' a table will appear in the Viewer of R studio
# method=c("bray", "euclidean", "canbrera", "jaccard")
effectEval <- function(bch, amp, dist = paste("bray"), strata_var=NULL,
                       tab.title = "Independent Effect Size Evals (% & sig) [ PERMANOVA & ANOSIM ]"){
  
  tbl_func <- function(x, amp_obj, strata=strata_var) {
    CMTmtxA <- vegdist(t(amp_obj$abund), method=dist)
    if(!is.null(strata_var)){strata=amp_obj$metadata[[strata_var]]}else{strata = strata_var}
    
    out  <- anosim(CMTmtxA,   as.factor(amp_obj$metadata[[x]]), strata= strata, parallel = 12) ; #out
    out2 <- adonis2(CMTmtxA ~ as.factor(amp_obj$metadata[[x]]), strata= strata, parallel = 12) ; #out2
    return(list(bw.R2pct= percent(out2$R2[1], accuracy=0.1, suffix = NULL),
                bw.Fval = round(out2$F[1],  3),
                bw.Pval = round(out2$`Pr(>F)`[1], 3),
                wi.Rpct = percent(out$statistic[1], accuracy=0.1, suffix=NULL),
                wi.Pval = round(out$signif[1],    3)        ) )  
  }
  szs=as.data.frame(t(sapply(bch, tbl_func, amp, strata_var)) )
  return(
    as.data.frame(szs) %>% dplyr::mutate_at(colnames(.), as.numeric) %>%
      dplyr::arrange(desc(bw.R2pct) ) %>%
      add_significance("bw.Pval", "bw.sig")  %>%
      add_significance("wi.Pval", "wi.sig" ) %>% relocate("bw.sig", .after = "bw.Pval") %>%
      kable(caption=tab.title) %>% #, label = "Test Label"
      add_header_above(c(" " =1,  "PERMANOVA" = 4, "ANOSIM" =3)) %>%
      column_spec(c(5,8), bold =T) %>%
      kable_styling(c("striped", "hover"))%>%
      scroll_box(width = "80%", height = "100%")
  )
}

#' -------------------- Modifiers vs Confounders -------------------------------------------
#' More info: 
#' https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704-ep713_confounding-em/bs704-ep713_confounding-em_print.html
#'
#' **Modifiers**:
#' Some interactions are biological phenomena (modifying interactions) producing a statistically significant \
#' modification of the community response to the factor of interest. Those interactions are biologically interesting \
#' and add detail into the understanding of the community response (cmtMTX) to the main factor (FCT). When there \
#' is effect modification, analysis of the pooled data (in global or stratum-specific measures) can be misleading. Therefore, a common way of dealing with \
#' such statistical interactions, is to explore cmtMTX ~ FCT response for **each level of MODIFIER variable (MODI)**.
#'  
#' **Confounders**
#' Other interactions are more of a byproduct from biases within the study groups (with respect to external variables) \
#' that can produce a distortion in the explored interactions between cmtMTX ~ FCT response (*confounding var, (Conf)*). 
#' Confounders need to be adjusted for in study designs and/or statistical models. 
#' 
#' A note on **Random Effects**
#' Random effects such as technical variations (PID or RunID or batch), are not to be explored as confounders nor modifiers \
#' As they are not part of study design, nor biological phenomena
#' 
#' --------------- Effect of modifier/ statistical interaction evaluation -----------------------------------------
#' Function performs regular adonis, but evaluates the response of the cmtMTX to FCT ~as modified by~ MODI
##'    Try to use only modifiers with 2+ levels

effectModi <- function(bch, amp, modi, dist = paste("bray"), strata_var=NULL){
  tab.title="Modifier Effect Evals [adonis2(CMTmtx ~ X*y) ]"
  tab=data.frame()
  CMTmtx <- vegdist(t(amp$abund), method = dist)
  
  if(is.null(strata_var)){ strata = NULL}else{strata=amp$metadata[[strata_var]] }
  
  for (x in bch){
    for (j in modi ) {
      try({
        amp$metadata <- amp$metadata %>% mutate_at(c(x, j), factor)
        pmn <- adonis2(CMTmtx ~ amp$metadata[[x]] * amp$metadata[[j]], strata=strata, parallel = 12); # print(pmn)
        z = cbind(effect = x, modifier = j, Pval=round(pmn$`Pr(>F)`[3], 4) )
        tab=rbind(tab, as.data.frame(z) )
      })  
    }
  }
  return(
    tab %>% mutate_at("Pval", as.numeric) %>% 
      rstatix::add_significance("Pval", "sig") %>%   #dplyr::filter(Pval <0.05) %>%
      kable(caption=tab.title) %>% 
      kable_styling(c("striped", "hover")) %>% column_spec(3, bold=T) %>%
      scroll_box(width = "80%", height = "100%") 
  )
}
#' If the factor on its own is significant, and its modification effect is also significant,
#' its also probably a confounder. But it would be a modifier


#' -------------- Effect of Confounding interaction evaluation -----------------------
#' Function performs logistic regression with anova(glm(, family = "binomial")), to find confounding factors to factors of interest. 
#'  Significance of the interaction is assessed with ChiSequared test. Is **important for X to be a binomial factor**, not continuous var
#'  The microbiome counts data is NOT actually used.

effectConf <- function(bch, amp, conf){
  tab=data.frame()
  tab.title = "Confounding Effects Evals  [anova(glm(X~y, binomial)), chisq]"
  for (x in bch){
    for (j in conf) { try({
      amp$metadata <- amp$metadata %>% mutate_at(c(x, j), factor)
      fml=as.formula(paste(x, "~", j))
      gmn = anova(glm(fml, data = amp$metadata, family="binomial"), test = "Chisq");  # print(gmn)
      z = cbind(main.factor = x, confounder = j, Pval=round(gmn$`Pr(>Chi)`[2], 4) )
      tab=rbind(tab,as.data.frame(z) )
    })      }
  }
  tab.kable <-  tab %>% mutate_at("Pval", as.numeric) %>% 
    add_significance("Pval", "sig", symbols = c("****", "***", "**", "*", "ns")) %>% 
    kable(caption=tab.title) %>% 
    kable_styling(c("striped", "hover")) %>%
    scroll_box(width = "80%", height = "100%")  
  
  return(  tab.kable  )
}


#' ------------------- effect adjustments for Confounders ------------------------------------------
effectEval_adjConf <- function(bch, amp, adj, dist = "bray", strata_var=NULL,
                               tab.title = "Adjusted Effect Size Evals (% & sig) [adonis2(CMTmtx ~ X+conf+... ,by='margin') ]" ){
  
  CMTmtx <- vegdist(t(amp$abund), method = dist)
  if(!is.null(strata_var)){ strata = amp$metadata[[strata_var]] }else{strata=strata_var}
  amp$metadata %<>% mutate_at(c(bch,adj), factor)
  
  tab=data.frame()
  for (x in bch){  try({
    fml <- as.formula(paste("CMTmtx ~", x, "+", paste(adj, collapse = " + ") ) )
    pmn <- adonis2(fml, by="margin", data = amp$metadata, parallel = 12, strata=strata);
    ans <-  anosim(CMTmtx, amp$metadata[[ x ]],           parallel = 12, strata= strata)
    z = cbind(effect = x, adjustment = paste(adj, collapse = " + "),
              bw.R2pct=percent(pmn$R2[1], accuracy=0.1, suffix = NULL),
              bw.Fval=round(pmn$F[1], 3),
              bw.Pval=round(pmn$`Pr(>F)`[1], 4),
              wi.Rpct = percent(ans$statistic[1], accuracy=0.1, suffix=NULL),
              wi.Pval = round(ans$signif[1],    3)            )
    tab=rbind(tab, as.data.frame(z) )
  }) }

    selcols <- colnames(tab)[-c(1,2)]
  tab.kable <- as.data.frame(tab) %>%  dplyr::mutate_at(selcols, as.numeric) %>% dplyr::arrange(desc(bw.R2pct) ) %>%
    rstatix::add_significance("bw.Pval", "bw.sig") %>% relocate("bw.sig", .after = "bw.Pval") %>% 
    rstatix::add_significance("wi.Pval", "wi.sig" ) %>% relocate("wi.sig", .after = "wi.Pval") %>%
    kable(caption=tab.title) %>%  add_header_above(., c(" " =2, "adjPERMANOVA" = 4, "ANOSIM" =3)) %>% 
    kable_styling(c("striped", "hover")) %>% column_spec(c(5,8), bold=T) %>%
    scroll_box(width = "80%", height = "100%")  
  
  return( tab.kable )
}



#' --------------- Pairwise Effect Size Evaluation -------------------------
#' need to find a way to give it a title

# effectPair <- function(bch, amp_obj, dist = paste("bray")){
#   CMTmtxA <- vegdist(t(amp_obj$abund), method=dist)
#   out <- pairwise.adonis(CMTmtxA, as.factor(amp_obj$meta[[bch]] ))
#   return(list(out )  )
# }

effectPair2 <- function(bch, amp_obj, dist = paste("bray"), strata=NULL){
  CMTmtxA <- vegdist(t(amp_obj$abund), method=dist)
  fml <- as.formula(paste0("CMTmtxA ~", bch))
  out <- pairwise.adonis2(fml, data=amp_obj$metadata, strata=strata ) 
  
  pvals <- as.data.frame(out)[1, ] %>% select(contains(c("Pr..F." )) ) %>% 
    rename_with(~str_remove(., ".Pr..F.") ) %>% t() %>%
    as.data.frame() %>% rstatix::add_significance(bch, "sig")  
  rvals <- as.data.frame(out)[1, ] %>% select(contains(c("R2" )) ) %>% 
    rename_with(~str_remove(., ".R2") ) %>%  t() %>% as.data.frame() %>% 
    rename_with( .cols=1, ~ "R2") %>% mutate(R2=percent(R2, accuracy=0.01))
  call <- as.data.frame(out)[1, ] %>% select(contains(c("call" )) ) 
  
  out1 <- cbind(rvals, pvals)
  out2<- list(out1, call)  
  names(out2) <- c("results", "call")
  return(out2)
}

#' -------------- Effect Dispersion Evaluation ------------------------------------
#' Function performs PERMDISP test on bch factors, using amp object
#' Both COUNTS and METADATA is needed from amp object
#' Since PERMANOVA and ANOSIM might take a while I make this separate function
#' 
#' to use function: 
#' source("~/Documents/MyApps/RScripts/myRscripts/effectEval.R")
#' amp <- data$DZs$amp
#' bch <- c("SampleType", "Disease"), etc
#' effectDisp(bch, amp)
#' a table will appear in the Viewer of R studio

effectDisp <- function(bch, amp, dist = paste0("bray")){
  tbl_func <- function(x, amp_obj) {
    CMTmtxA <- vegdist(t(amp_obj$abund), method=dist)
    tab.title <- "Evaluating Effect Homogeneity [ PERMDISP ]"
    bdis <- betadisper(CMTmtxA, as.factor(amp_obj$metadata[[x]]))
    bdisp <- vegan::permutest(bdis, pairwise = T, parallel = 12)
    return(list(disp.Size = round(bdisp$tab[1,4], 3),
                disp.Pval = round(bdisp$tab[1,6], 3)    ) )  
  }
  szs=as.data.frame(t(sapply(bch, tbl_func, amp)) )
  return(
    as.data.frame(szs) %>% dplyr::mutate_at(colnames(.), as.numeric) %>%
      dplyr::arrange(desc(disp.Size) ) %>%
      add_significance("disp.Pval", "sig") %>% relocate("sig", .after = "disp.Pval") %>% 
      kable(caption=tab.title) %>%
      add_header_above(c(" " =1,  "DISPERSION" = 3)) %>%
      column_spec(c(4), bold =T) %>%
      kable_styling(c("striped", "hover"))%>%
      scroll_box(width = "40%", height = "100%")
  )
}
#' 	- The PERMITATIONAL dispersion test, is technically only applicable if PERMANOVA
#' 	 is significant and separate clusters are not observed on the PCoA ordination.
#'  - (because if clusters are observably separated on the PCoA/nMDS, => the differences 
#'  between groups that PERMANOVA detects, are location-based, duh!)
#'  - These differences could ALSO be based on dispersion of the clusters 
#'  (how tight/dispersed the clusters are from the centroid of their cluster)
#'  
#'  PERMDISP only tests homogeneity of dispersion of the clusters, but location-based differences can be inferred from its stat output

#'  Condition 1: PERMANOVA needs to be significant (differences between groups need to be detected!)

#'  PERMDISP: 
#'  - if p<0.05 and F>>2 => very significant  dispersion differences, and it is possible location-based differences are not observed
#'  - if p<0.05 and F<<2 =>  dispersion is significant, but it is also possible that location-based differences are observed
#'  
#'  - if p>0.05 and F<<2 => the clusters are homogenous (have similar dispersion) and it is likely that any 
#'  differences between groups (b/c PERMANOVA is sig) are due to location separation of clusters
#'  - if p>0.05 and F>>2 => what are we doing here? your groups of samples are not only 
#'  homogeneously distributed, but also likely overlap. It is very likely your PERMANOVA 
#'  was not significant 



#' --------------- Nested Effect evaluation -----------------------------------------
#' Function performs regular adonis, but evaluates 
#'    the response of the microbiome to the effect of any factors of interest 
#'    ~across each level of~  any nest factor
#'    Function uses both COUNTS.mtx and factors of interest
#'    Try to use only modifiers wtih >2 levels
#'    
#'    About the by term: 	
#'  by = "terms" will assess significance for each term (sequentially from first to last), 
#'  setting by = "margin" will assess the marginal effects of the terms, where each marginal
#'  term is analysed in a model with all other variables),
#'  and by = NULL will assess the overall significance of all terms together.
#'    
#' to use function: 
# source("~/Documents/MyApps/RScripts/myRscripts/effectConf.R")
# amp=data$DZs$amp               # amp: uses both abundance and meta data
# bch=fct$all[c(1,2,4)]          # factors of interest
# nest <- fct$all[7]             # the nesting factors (modifiers),
#' where we already know nest/modifier has a significant effect
# effectNest(bch, amp, modi)
effectNest <- function(bch, amp, nest, dist = paste0("bray")){
  tab=data.frame()
  CMTmtx <- vegdist(t(amp$abund), method = dist)
  for (x in bch){
    for (j in nest ) {
      try({
        amp$metadata <- amp$metadata %>% mutate_at(c(x, j), factor)
        pmn <- adonis2(CMTmtx ~ amp$metadata[[j]] / amp$metadata[[x]]); # print(pmn)
        # pmn <- adonis2(CMTmtx ~ amp$metadata[[x]] + amp$metadata[[j]],
                       # by="margin", parallel = 12); # print(pmn)
        z = cbind(effect = x, nest = j, 
                  eff.Fval=round(pmn$F[1], 3), 
                  eff.Pval=round(pmn$`Pr(>F)`[1], 4),
                  nst.Fval=round(pmn$F[2], 3)      ) #,
                  # nst.Pval=round(pmn$`Pr(>F)`[2], 4)   )
        tab=rbind(tab, as.data.frame(z) )
      })  
    }
  }
  return(
    tab %>% mutate_at("eff.Pval", as.numeric) %>% 
      rstatix::add_significance("eff.Pval", "sig") %>% relocate("sig", .after = "eff.Pval") %>% 
      kable(caption="Eval. Effect Size% & Sig ~across~ Nests [adonis2(CMTmtx ~ y/X ]") %>% 
      kable_styling(c("striped", "hover")) %>% column_spec(5, bold=T) %>%
      scroll_box(width = "80%", height = "100%") 
  )
}


#' -------------- Effect Size and Significance Evaluation (with adjustments) ------------------------------------
#' Function performs ~~anosim and~~ adonis (PERMANOVA) on metadata factors, using amp object
#' Both COUNTS and METADATA is needed from amp object
#' 
#' PERMANOVA tests whether distance differ between groups.
#' ANOSIM tests whether distances between groups are greater than within groups. =>
#' If var w/i group <   var b/w group, ANOSIM:: SIG (bw/wi >   1, SIG) \
#' If var w/i group >/= var b/w group, ANOSIM:: NS  (bw/wi </= 1, NS)
#' methods of choice = c("bray", "euclidean", "canbrera", "jaccard")

#' ensure that the continuous/categorical variables are set as as.numeric()/as.factor() before use of function. 
#' if the degrees-of-freedom of PERMANOVA for a numeric factor are NOT == 1, the numeric factor is calculated as categorical (wrong). 
#' ensure degrees-of-freedom==1 for numeric factor PERMANOVAs

#' If a variable has SIG influence on the microbiome (effectEval) &AND& significant influence on the factor of interest (effectConf),
#'  then it is a **confounder** to the factor of interest and needs to be corrected for. 

effectEval <- function(bch, amp, adj= NULL, dist = "bray", strata_var = NULL, tab.title = NULL ){
  if(is.null(tab.title)){
    title="Microbiome responce to variable (+adjustments)"
  }else{ title= tab.title }
  labl = "[PERMANOVA CMTmtx ~ X + adjmt , by = 'margin']"  
  
  
  if(!is.null(strata_var)){ stra = amp$metadata[[strata_var]] }else{stra=strata_var}
  # amp$metadata %<>% mutate_at(c(bch,adj), factor)
  
  tab=data.frame()
  adjmt <- paste(adj, collapse = " + ")
  
  # pmn <- matrix(nrow = 3, ncol = 6) %>% as.data.frame %>% 
  #   dplyr::rename_with( ~ c("Df", "SumSq", "R2", "F", "Pr(>F)", "sig")) %>% 
  #   { .[1, ] <- 1; . } ## creating empy pmn with fake 1 values, for when adonis fails due to no groups in X 
  # 
  
  for (x in bch){  try({
    if(is.numeric(amp$metadata[[x]])){  ## sets categorical vs continuous variable
      ftype="continuous"; amp$metadata[[x]] %<>% as.numeric()
      }else{
      ftype="categorical"; amp$metadata[[x]] %<>% as.factor()  }
    print(paste("----- evaluating ", x, " as ", ftype))
    
    
    amp <- amp_subset_samples(amp, !is.na(amp$metadata[[x]]) )                     ## removes samples with NA for metadata[[x]]
    CMTmtx <- vegdist(t(amp$abund), method = dist, na.rm = T)
    

    if (is.null(adj)){ 
      fml <- as.formula(paste("CMTmtx ~", x ) ); adjmt <- "none"
    }else{  
      fml <- as.formula(paste("CMTmtx ~", x, "+", adjmt ) ) }
    
    pmn <- adonis2(fml, data = amp$metadata, parallel = 8, strata=stra, by="margin")
    # print(pmn)
    # ans <-  anosim(CMTmtx, amp$metadata[[ x ]],            parallel = 8, strata= stra)
    
    
    z = cbind(`variable` = x, `+ adjmts` = adjmt, bw.Df=pmn$Df[1],
              `bw.R2%`=percent(pmn$R2[1], accuracy=0.1, suffix = NULL),
              bw.Fval=round(pmn$F[1], 3), bw.Pval=round(pmn$`Pr(>F)`[1], 4) #,
              # wi.Rpct = percent(ans$statistic[1], accuracy=0.1, suffix=NULL),
              # wi.Pval = round(ans$signif[1],    3)
    )
    tab=rbind(tab, as.data.frame(z) )
  }) }
  
  selcols <- colnames(tab)[-c(1,2)]
  tab.kable <- as.data.frame(tab) %>%  dplyr::mutate_at(selcols, as.numeric) %>% dplyr::arrange(desc(`bw.R2%`) ) %>%
    rstatix::add_significance("bw.Pval", "bw.sig") %>% relocate("bw.sig", .after = "bw.Pval") %>% 
    # rstatix::add_significance("wi.Pval", "wi.sig" ) %>% relocate("wi.sig", .after = "wi.Pval") %>%
    kable(caption=title, label = labl) %>%  add_header_above(., c(" " =2, "PERMANOVA" = 5)) %>% ##, "ANOSIM" =3)) %>% 
    kable_styling(c("striped", "hover")) %>% column_spec(c(6), bold=T) %>% ##c(5,8)
    scroll_box(width = "70%", height = "90%")  
  
  return( tab.kable )
}



#' --------------- Evaluate Effect Changes (modifiers vs confounders ) ---------------------------
#' This function tests for changes in effect size (F) away from "self" (shouldnt it be R2)?
#' If F changes more than 10% when varB is added, varB is a confounding factor.
#' For modifier we are looking at the p-value in the "modi" test. If the p-value is smaller, varB is a modifier
 

effectTests <- function(mainVar, amp, testVar, dist = "bray", strata_var=NULL, seed=123, test_type="both"){
  tab=data.frame()
  CMTmtx <- vegdist(t(amp$abund), method = dist)
  
  if(is.null(strata_var)){ stra = NULL}else{stra=amp$metadata[[strata_var]] }
  
  for (x in mainVar){
    for (j in testVar ) {
      try({
        # amp$metadata <- amp$metadata %>% mutate_at(c(x, j), factor)
        set.seed(seed)
        
        ## checks mainVar effect
        pmn <- adonis2(CMTmtx ~ amp$metadata[[x]],                     strata=stra, parallel = 8, by = "margin"); 
        ## checks confounder effects of testVar on mainVar
        pcf <- adonis2(CMTmtx ~ amp$metadata[[x]] + amp$metadata[[j]], strata=stra, parallel = 8, by = "margin"); 
        ## checks modifier effects of testVar on mainVar
        pmd <- adonis2(CMTmtx ~ amp$metadata[[x]] * amp$metadata[[j]], strata=stra, parallel = 8, by = "margin"); 
        ## grabs interaction row.number (contains ":"), if exists, or returns NA
        intR <- grep(":", row.names(pmd), value = F);  intR <- ifelse(is.integer(intR), intR, NA) ;  ## for grabbing the interaction row, if it exists
        
        ## getting effect sizes:
        # mVAR <- grep(x, row.names(pmn), vaue = F) ## but mVAR row number will always be 1 in this model setup
        pmn_R2p <- percent(pmn$R2[1], accuracy=0.01, suffix = NULL); pmn_Fval = round(pmn$F[1],2); pmn_pval = round(pmn$`Pr(>F)`[1], 4)
        pcf_R2p <- percent(pcf$R2[1], accuracy=0.01, suffix = NULL); pcf_Fval = round(pcf$F[1],2); pcf_pval = round(pcf$`Pr(>F)`[1], 4) 
        pmd_R2p <- percent(pmd$R2[1], accuracy=0.01, suffix = NULL); pmd_Fval = round(pmd$F[intR],2); pmd_pval = round(pmd$`Pr(>F)`[intR], 4) # grabs intR
        
        ## binding results
        m = cbind(mainVAR = x, testVAR = "", test.for="self", R2p = pmn_R2p, Pval = pmn_pval) %>% as.data.frame 
        c = cbind(mainVAR = x, testVAR = j , test.for="conf", R2p = pcf_R2p, Pval = pcf_pval) %>% as.data.frame
        f = cbind(mainVAR = x, testVAR = j , test.for="modi", R2p = pmd_R2p, Pval = pmd_pval) %>% as.data.frame
        b = cbind(mainVAR = "----", testVAR = "----" , test.for="----", R2p= "----", Pval=NA ) %>% as.data.frame
        # m = cbind(mainVAR = x, testVAR = "", test.for="self", Fval=pmn_Fval , Pval = pmn_pval ) %>% as.data.frame 
        # c = cbind(mainVAR = x, testVAR = j , test.for="conf", Fval=pcf_Fval , Pval = pcf_pval ) %>% as.data.frame
        # f = cbind(mainVAR = x, testVAR = j , test.for="modi", Fval=pmd_Fval , Pval = pmd_pval) %>% as.data.frame
        # b = cbind(mainVAR = "----", testVAR = "----" , test.for="----", Fval= "----", Pval=NA )   %>% as.data.frame
        tab=rbind(tab,  m,c,f, b)
      })  
    }
  }
  # print(tab)
  if(test_type != "both"){ #subsets table to test type & removes redundant rows
    tab %<>% filter(test.for %in% c("self", test_type )) %>% distinct(mainVAR, testVAR, test.for, .keep_all=T)
  }

  bd <- ifelse(test_type=="modi", 5, ifelse(test_type=="conf", 4, c(4,5) ) )
  tab.title="Evals Effect Interactions [adonis2(CMTmtx ~ mainVar {+,*} testVar) ]"
  
  return(
    tab %>% mutate_at("Pval", as.numeric) %>% 
      rstatix::add_significance("Pval", "sig") %>% 
      kable(caption=tab.title) %>% 
      kable_styling(c("striped", "hover")) %>% column_spec(bd, bold=T) %>%
      scroll_box(width = "100%", height = "100%") 
    # add_row(mainVAR = "-------", testVAR = "-----", effect = "----", Pval = NA , sig = "") %>%
    
  )
}

#' -------------------- Modifiers vs Confounders -------------------------------------------
#' More info: 
#' https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704-ep713_confounding-em/bs704-ep713_confounding-em_print.html
#' 
#' **Covariables /biologically non-interesting interactions (e.g. demographics) ** 
#' **can be adjusted for as their influence is not of importance**
#' Usually things like the demographic variables (age, race, bmi, gender), they a formed by biases within the 
#' study groups.  They usually significantly covary with the main factors of interest within the dataset. 
#' Correction for them in PERMANOVA and other models, usually increases the accuracy of the study. 
#' They can produce distortion in the explored interactions between cmtMTX ~ FCT response, but generally 
#' they are "less serious confounders" (PERMANOVA MCV ~ FCT + confounders + covariables) 
#' 
#' **Confounders (serious, independent variable - e.g. antibiotics)**
#' These are also a byproduct from biases within the study groups, but can produce serious distortion in the 
#' explored interactions between microbial community variation (MCV) ~ FCT response. Confounders need to be
#' adjusted for in study designs & statistical models. (PERMANOVA MCV ~ FCT + confounders)
#' How to determine: in PERMANOVA model MCV ~ varA + varB, if outcome of varA sig. changes when varB is added,
#' then varB is confounder of varA.
#'  
#' **Modifiers (serious, dependent or biological variable - e.g. opportunistic infection)**:
#' Some interactions are biological phenomena (modifying interactions) producing a statistically significant \
#' modification of the community response to the factor of interest. Those interactions are biologically interesting \
#' and add detail into the understanding of the community response (MCV) to the main factor (FCT). Modifiers are 
#' identified through model **(PERMANOVA MCV~FCT*modifier)**. When there is significant effect modification, analysis 
#' of the pooled data (in global or non-stratum-specific measures) can be misleading. Therefore, a common way of 
#' dealing with such statistical interactions, is to explore MCV ~ FCT response for **each level of MODIFIER variable**
#' ** individually (MODI)**. Another way to identify is - exploration of MCV~FCT at each strata (if at each strata 
#' MCV~FCT was NS, then PERMANOVA MCV~FCT*modifier was NS)
#' How to determine: in PERMANOVA model MCV ~ varA * varB, if outcome of varA sig changes when varB is added,
#' then varB is modifier to varA
#'  
#' 
#' 
#' A note on **Random Effects**
#' Random effects such as technical variations (PID or RunID or batch), are not to be explored as confounders nor modifiers \
#' As they are not part of study design, nor biological phenomena
#' 
#' 
#' 
#' 
#' 
#' --------------- Effect of modifier/ statistical interaction evaluation -----------------------------------------
#' Function performs regular adonis, but evaluates the response of the cmtMTX to FCT ~as modified by~ MODI
##'    Try to use only modifiers with 2+ levels
#' If the factor on its own is significant, and its potential modifier is also significant,
#' the potential modfier is probably also a confounder.

effectModi <- function(bch, amp, modi, corr = NULL , dist = paste("bray"), strata_var=NULL){
  tab.title="Modifier Effect Evals [adonis2(CMTmtx ~ X*y+corr) ]"
  tab=data.frame()
  CMTmtx <- vegdist(t(amp$abund), method = dist)
  
  if(is.null(strata_var)){ strata = NULL}else{strata=amp$metadata[[strata_var]] }
  
  for (x in bch){
    for (j in modi ) {
      try({
        # amp$metadata <- amp$metadata %>% mutate_at(c(x, j), factor)
        rr <- if(is.null(corr)){NULL}else{paste("+", corr) }
        
        fml0 <- as.formula(paste("CMTmtx", "~" ,x, rr) )
        pmn0 <- adonis2(fml0, data = amp$metadata, strata=strata, parallel = 12, by = "margin")
        l0  <- grep(x, row.names(pmn0)); # greps for the main factor line (corrected or not)
        
        fml  <- as.formula(paste("CMTmtx", "~",x ,"*", j , rr) )
        pmn <- adonis2(fml, data = amp$metadata, strata=strata, parallel = 12, by = "margin")
        l <- grep(":", row.names(pmn)); l<- ifelse(is.integer(l), l, NA)
        adjr <- ifelse(is.null(corr), "none", corr)
        z0= cbind(effect = x, modifier = "none", adjustor = adjr,  Pval=round(pmn0$`Pr(>F)`[l0], 4) )
        z = cbind(effect = x, modifier = j,      adjustor = adjr,  Pval=round( pmn$`Pr(>F)`[l ], 4) )
        tab=rbind(tab, as.data.frame(z0), as.data.frame(z) )
      })  
    }
  }
  # tab <- tab[!duplicated(tab$modifier), ]
  tab %<>% distinct(modifier, .keep_all=T)
  return(
    tab %>% mutate_at("Pval", as.numeric) %>% 
      rstatix::add_significance("Pval", "sig") %>%   #dplyr::filter(Pval <0.05) %>%
      kable(caption=tab.title) %>% 
      kable_styling(c("striped", "hover")) %>% column_spec(4, bold=T) %>%
      scroll_box(width = "80%", height = "100%") 
  )
}



#' -------------- Effect of Confounding interaction evaluation -----------------------
#' Function performs logistic regression with anova(glm( X ~ j, family = "binomial"), test="Chisq"), 
#' to find potential confounding factors (J) to factors of interest (X). 
#' This function **only applies where X is a binomial (2-level categorical factor)**, not continuous (e.g. Shannon)
#' the potential confounder (J) is free to be either numeric or categorical with multiple levels. 
#' Significance of the interaction is assessed with ChiSequared test, which applies only to glm()
#' Potential confounding factors are to be evaluated only if they are also SIG contributors to CMvar. 
#' *The microbiome counts data is NOT actually used*.

effectConf_adonis <- function(bch, amp, conf, sigFLT=0.065, type = "binomial"){
  tab=data.frame()
  
  for (x in bch){
    for (j in conf) { try({
      amp$metadata <- amp$metadata %>% mutate_at(c(x, j), factor)
      fml=as.formula(paste(x, "~", j))
      if (type == "binomial"){
        tab.title = "Confounding Effects Evals  [anova(glm(X~y, binomial)), chisq]"
        gmn = anova(glm(fml, data = amp$metadata, family="binomial"), test = "Chisq"); ## for binomial X values 0<=X<=1
        n=grep(j, rownames(gmn), value =F)
        pvl <- gmn$`Pr(>Chi)`[n];# print(pvl)
      } else if (type == "multinomial"){
        tab.title = "Confounding Effects Evals  [anova(multinom(X~y, multinomial)), chisq]"
        gmn = car::Anova(nnet::multinom(fml, data = amp$metadata), type = 3)  ## for multinomaial X
        n=grep(j, rownames(gmn), value =F)
        pvl <- gmn$`Pr(>Chisq)`[n]; #print(pvl)
      } else { print("---Chose type 'binomial' or 'multinomial' ") }
      
      z = cbind(main.factor = x, confounder = j, Pval= formatC(pvl, format = "e", digits = 2) )
      tab=rbind(tab,as.data.frame(z) )
      # print(tab)
    })      }
  }
  tab <- tab %>% mutate(Pval = as.numeric(Pval)) %>% filter(Pval <= sigFLT) %>%
    add_significance("Pval", "sig", symbols = c("****", "***", "**", "*", "ns"))
  tab.kable <- tab %>% mutate(Pval = as.character(Pval)) %>%
    kable(caption=tab.title) %>%
    kable_styling(c("striped", "hover")) %>%
    scroll_box(width = "80%", height = "100%")

  return(  tab.kable  )
}

#' ------------------- testing candidate confounders with anova(lmer()) ---------------------------
#' meant for longitudinal dataset (having a timeVar) to test **responceVAR ~ mainVAR*timeVAR + potConf + (1|SubjID)**
#' where reponceVAR can be a PCo1.. or Alpha index (like Shannon or SppRichness)
#' The function explores the confounding effects of a list of potential confounders (extra_term), 
#' on the **mainVAR*timeVAR** interaction. In the results table, "baseline" is the p-value for mainVAR*timeVAR interaction alone, 
#' and all other p-values are affected by the added candidate confounder (extra_term) to the interaction formula.
#' If p-value changes a lot, or significance dissappears, the extra_term added is a confounder
#' The input dataset is a beta_ordination_plot object from beta_plot() function
#' the find_str is defaulted at ":" for function to report the interactive main_term
#' if non-interaction model is used (dont), set find_str to the mainVAR string
#' 
#' If the ordi$data has alpha diversities, change the responceVAR to one of the alpha indices and you get alpha div confounding

effectConf_lmer <- function(ordi_data = ordiDataObj, candidate_confounders_list= c("Sex", "Ethnicity"), responceVAR = "PCo1",
                            find_str = ":", main_term = "mainVAR*timeVAR", rand_var = "(1|SubjectID)" ) {
  
  # Function: fit model and return interaction results
  fit_with <- function(extra_term = NULL, label = "baseline", pco = responceVAR, data = ordi_data) {
    base_fml <- as.formula(paste(pco , "~",                  main_term, "+", rand_var )) 
    conf_fml <- as.formula(paste(pco , "~", extra_term, "+", main_term, "+", rand_var ))
    fml <- if (is.null(extra_term)) { base_fml } else { conf_fml }
    m <- lmer(fml, data = ordi_data)
    an <- anova(m) %>% as.data.frame() %>%  ## stats::anova does sequential but main_term is last so OK
      tibble::rownames_to_column("term") %>%
      filter(str_detect(term, find_str ) ) %>% mutate(model = label)
  }
  
  # Run baseline + confounder models
  results <- bind_rows(
    fit_with(NULL, "baseline"),
    map_dfr(candidate_confounders_list, ~ fit_with(.x, .x))# %>% 
    # add_significance(p.col = "Pr(>F)", output.col = "sig")
  )
  
  return(results)
}

###' examples
###' effectConf_lmer(data = ordi$data, responceVAR = "PCo3", candidate_confounders_list = confounders, main_term = "PMA * Cohort2")
###' effectConf_lmer(data = ordi$data, responceVAR = "PCo3", candidate_confounders_list = confounders, main_term = "PMA + Cohort2", find_str = "Cohort2")
###' effectConf_lmer(data = ordi$data, responceVAR = "Shannon", candidate_confounders_list = confounders, main_term = "PMA * Cohort2")



#' --------------- Pairwise Effect Size Evaluation -------------------------
#' calculates pairwise contrasts between groups of a variable
#' I have made it lapply-style to accommodate adjustments into the model 
#' (could also be done with for loop, but i want to learn to be fancy)
#' (since pairwise.adonis2() could not take in adjustments into the model)
#' 


effectPair2 <- function(bch, amp_obj, adj = NULL,  dist = "bray", strata_var = NULL){
  
  metadat <- amp_obj$metadata
  abundat <- amp_obj$abund
  CMTmtxA <- as.matrix(vegdist(t(abundat), method=dist))
  
  adj_term <- ifelse(is.null(adj), "", paste("+", paste(adj, collapse = "+")) )
  fml <- as.formula(paste0("sub_dist ~", bch, adj_term))
  
  groups <- unique(metadat[[bch]])
  pairs <- combn(groups, 2, simplify = F)
  pairsR <- sapply(pairs, \(x) paste(x, collapse = " ~ "))
  
  ## calculates adonis and collates the results
  results <- lapply(seq_along(pairs), function(i) {
    p <- pairs[[i]]
    n <- pairsR[i]
    
    sub_meta <- metadat[metadat[[bch]] %in% p, ]
    sub_dist <- CMTmtxA[rownames(sub_meta), rownames(sub_meta)]
    if(is.null(strata_var)){sub_stra <- NULL}else{sub_stra <- sub_meta[[strata_var]]}
    
    ads <- adonis2(fml, data = sub_meta, strata = sub_stra, parallel = 12, by = "margin")
    pval <- ads[get("bch"), 5]; rval <- percent(ads[get("bch"),3], accuracy = 0.11)
    data.frame(contrastedPairs = n, R2prc = rval , pVal = pval) 
  })
  out2 <- bind_rows(results) %>% add_significance("pVal", "sig", 
                                                  cutpoints = c(0,      1e-04, 0.001, 0.01, 0.05, 0.065, 1),
                                                  symbols =   c("****", "***",  "**",  "*",  ".", "ns"))
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
    bdis <- betadisper(CMTmtxA, as.factor(amp_obj$metadata[[x]]))
    bdisp <- vegan::permutest(bdis, pairwise = T, parallel = 12)
    return(list(disp.Size = round(bdisp$tab[1,4], 3),
                disp.Pval = round(bdisp$tab[1,6], 3)    ) )  
  }
  tab.title <- "Evaluating Effect Homogeneity [ PERMDISP ]"
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

### working on LM with mixed effects (tutorial)

#### A little tutorial on lmerTest::lmer()

# - anova(lmer()) vs summary(lmer())
# ANOVA is for testing lmer() with continuous or categorical variables (like timepoint or estrogen levels)
# SUMMARY is for testing lmer() models with continuous only variables (), and gives more info
# 
# - the random effect term
# Usually + (1|SubjID). Required
# TimePoint is NOT a random term. its longitudinal term
# 
# - How to find longitudinal effect
# do **anova(lmer(PCoN + TimePoint + (1|PID))**
#          
#          - how to know which model better fits the data:
#            lets say we have:
#             m1<- lmer(PCoN ~ mainVar+TP + (1|PID) ) [additive model]
#             m2<- lmer(PCoN ~ mainVar*TP + (1|PID) ) [interactive model]
#          
#          **method 1**: check interaction
#             do **anova(m2)**: 
#             - if interaction is SIG => m2 the more complex model is significantly better. 
#             - if interaction NS, m1 (simpler, additive model) is preferred due to parsimony (simplicity)
#          
#          **method 2**: check model fit 
#           do **AIC(m1, m2)** and **BIC(m1,m2)**
#            - the lower AIC or BIC value represents the better model fit
#          
#          - how to know if variable is a confounder and needs to be corrected for 
#             setup:
#               m1 <- lmer(PCoN ~ mainVAR + (1|PID) ) [no confounder]
#               m2 <- lmer(PCoN ~ mainVAR + pConf + (1|PID)) [added confounder]
#             do **anova(m1, m2)**
#               - if anova(m1,m2) is SIG => adding the pConf is necessary (pConf is a confounder)
#               - if anova(m1,m2) is NS, adding the pConf is not needed
#             or do **compare the F-values** of m1 and m2 
#               - if m2$F >10% off from m1$F, then pConf is confounder.


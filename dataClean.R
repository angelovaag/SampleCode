
#' -------------------- explore Depths --------------------------------------
#' Requires SampleID column in metadata, even if not called in fct1 or fct2
#' both fct1 and fct2 are customizable
exploreDepths <- function(amp, fct1="SampleID", fct2="SampleType", title="",
                          stepsize = 1e3, xmin = NULL, id.col="SampleID"){
  # rarecurve
  xmax=max(sort(colSums(amp$abund)))
  if(is.null(xmin)){xmin=min(unlist(lapply( colSums(amp$abund), \(x) x[x!=0]  ) ))  }
  suppressWarnings({
    rarecurve <- amp_rarecurve(amp, color_by = fct1, stepsize = stepsize , facet_by = fct2) + 
      xlim(0,xmin) + guides(color=guide_legend(ncol=12)) +
      labs(title=paste("Rarecurve.", title), x="Reads (seqDepth)", y="ASVs (sppRichness)")   }) #+ 
  # geom_text(aes(label = SampleID ), cex = 3, show.legend = F, check_overlap=T,
  # position = position_dodge(width=1)) # this works but is ugly 
  
  
  # octave plot 1
  octave1=amp_octave(amp, group_by = fct1, tax_aggregate = "Species", num_threads = 4) + 
    labs(title = paste("Octave plot ~", fct1, ".", title ))  
  # octave plot 2
  octave2=amp_octave(amp, group_by = fct2, tax_aggregate = "Species", num_threads = 4) +
    labs(title = paste("Octave plot ~", fct2, ".", title ))    
  
  # seqDepth info (based on fct2 categories)
  seqDepths=as.data.frame(sort(colSums(amp$abund))) %>% 
    rename_at(1, ~"seqDepth") %>% rownames_to_column("SampleID") %>% 
    mutate("{fct2}" := amp$metadata[[fct2]][match(SampleID, amp$metadata[[id.col]] )])
  
  # standard seqPlot no point in plot if its not based on SampleID
  # aes(x= .data[[fct1]], y=seqDepth)) & # aes_string(x= fct2, y="seqDepth"))  # work too
  seqPlot <- ggplot(seqDepths, aes(x=SampleID, y=seqDepth)) +
    geom_point(color="black") + theme_bw() + 
    theme(text=element_text(size=12), axis.text.x = element_text(angle=90));
  
  # return output
  all = list()
  all <- list(rarecurve, octave1, octave2, seqDepths, seqPlot)
  names(all) = c("rare", "octave1", "octave2" , "seqDepths", "seqPlot")
  
  return(all)
}


#' ------------------ my dataset CLEANING/ FLT function ------------------------------------------
#' for cleaning datasets, flt asvs and samples
#' provide AMP object, seqDepth threshold, asvPRESence in %, asvPREValence in % of samples, aggregation rank

myFilter <- function(amp, seqDepth= 10e3, asvPRES = 0.01, asvPREV = 3, agg_at = "Species"){
  
  # seqDepth at 10K reads, PRESENCe at 0.01%, PREVALENCE at 3%, aggregate and count Species
  spp_raw <- aggregate_abund(amp$abund, amp$tax, tax_aggregate = agg_at, calcSums = T, format= "abund")
  vbf.spp <- nrow(spp_raw)
  # try({
  suppressMessages({    
    amp_flt=amp_subset_samples(amp, minreads=seqDepth, removeAbsentOTUs = F, normalise =F) 
    })  
  # })
  # if (!exists("amp_flt")){amp_flt <- amp} 

  sbf=ncol(amp$abund) # number ssamples before flt
  saf=ncol(amp_flt$abund) # samples after flt
  srem<- sbf - saf
  
  # ========= ASVs presence-based filtering (relative abundance)
  vbf <- nrow(amp$tax) #asvs before filtering
  suppressMessages({   amp_flt <- filter_species(amp_flt, filter_otus = asvPRES )   })
  vaf.PRES <- nrow(amp_flt$tax)
  vrem <- vbf - vaf.PRES
  
  # ===== ASVs prevalence-based filtering (percent prevalence)
  asvPRV = round(ncol(amp_flt$abund)*(asvPREV/100), 0 ) # converting % in # of samples
  prvdf =  apply(amp_flt$abund, MARGIN=1, FUN = function(x){sum(x>0)}) # calc prev# for each ASV
  amp_flt$abund = amp_flt$abund[ prvdf >= asvPRV , ]#;   #(1)          # filter based on prev#
  amp_flt$tax <- amp_flt$tax[ row.names(amp_flt$abund), ] # cleaning TAX table
  vaf.PREV <- nrow(amp_flt$abund)
  vrem2 <- vaf.PRES - vaf.PREV 
  
  # ========= aggregating abundances at agg_at level (abund table only)
  amp_flt$agg=aggregate_abund(amp_flt$abund, amp_flt$tax, tax_aggregate = agg_at, calcSums = T, format= "abund")
  vaf.spp <- nrow(amp_flt$agg)
  
  
  # ===== stats
  flt_info <- paste("RAW dataset started with",  sbf, "samples,  represented by",  vbf, "ASVs (",  vbf.spp, "spp). \n",  
                    "- seqDepth filtering at",  seqDepth, "reads per sample,  removed",  srem, "samples,  retaining",  saf, "samples. \n", 
                    "- ASV presence filtering at",  asvPRES, "% relative abundance,  removed",  vrem, "ASVs,  retaining",  vaf.PRES, "ASVs. \n", 
                    "- ASV prevalence filtering at",  asvPREV, "% across all samples (<",  asvPRV, "samples),  removed additional",  vrem2, "ASVs. \n",  
                    "Finally,  \n",  
                    "FLT dataset is represented by",  saf, "samples with",  vaf.PREV, "well represented ASVs and",  vaf.spp,  "Species.\n")
  
  out <- list(amp_flt, flt_info)
  names(out) <- c("FLTamp", "FLTinfo" )
  return(out)
  
}


#' -------------------- amp aggregate Species -----------------------------------

amp_agg_abund <- function(amp, agg_at = "Species"){
      taxRanks=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      cutRanks=taxRanks[1:which(taxRanks==agg_at)]
      library(ampvis2)
      agg_spp <- aggregate_abund(amp$abund, tax= amp$tax, calcSums = T,
                                 tax_aggregate = agg_at,
                                 tax_add = cutRanks, format = "abund")
      # agg_spp %>% head
      agg_spp$Taxonomy <- row.names(agg_spp)
      row.names(agg_spp) <- paste0("tax", 1:nrow(agg_spp))
      # agg_spp %>% head
      
      agg_abu <- agg_spp %>% select(-Taxonomy)
      agg_tax <- agg_spp %>% select(Taxonomy) %>% 
        separate( col = "Taxonomy", into=cutRanks, sep="; ")
      # agg_tax
      
      amp_agg <- suppressWarnings( amp_load(agg_abu, amp$metadata, agg_tax) ) 
      amp_agg$biom <- agg_spp
      return(amp_agg)
}
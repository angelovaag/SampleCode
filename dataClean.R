
##' --------- explore Depths --------------------------------------
#' Requires SampleID column in metadata, even if not called in fct1 or fct2
#' both fct1 and fct2 are customizable
exploreDepths <- function(amp, fct1="SampleID", fct2="SampleType", title="", prop=1e-7,
                          stepsize = 1e3, xmax_quantile= 0.75, id.col="SampleID"){
  ## xmax_quantile is the xmax at the % of max reads (when some samples are much higher 
  ## seq depth than others)
  # rarecurve
  # if(is.null(xmax)){xmax=min(unlist(lapply( colSums(amp$abund), \(x) x[x!=0]  ) ))  }
  xmax=quantile(unlist(lapply( colSums(amp$abund), \(x) x[x!=0] )), xmax_quantile) 
  suppressWarnings({
    rarecurve <- amp_rarecurve(amp, color_by = fct1, stepsize = stepsize , facet_by = fct2) + 
      xlim(0,xmax) + guides(color=guide_legend(ncol=12)) + 
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
  
  # Katie's Read filtering plots (tax distribution plot)
  # y-lab is how many TAXa/ASVs fall in the category of being represented by X-axis # of reads
  ## the base you want to reduce, the [stepenta] you want to increase to move xint to left.
  taxSums <- rowSums(amp$abund)
    # prop <- 0.0000001
    xint <- prop*sum(taxSums)
    taxPlot <- ggplot(data.frame(taxSums), aes(x=taxSums)) + 
              geom_histogram(bins=50, fill="cornflowerblue", color="black") + theme_bw() +
              xlab("Taxa Sums (log10 scale)") + scale_x_log10() +
              geom_vline(xintercept = xint, color="red")
    
  # return output
  all = list()
  all <- list(rarecurve, octave1, octave2, seqDepths, seqPlot, taxPlot)
  names(all) = c("rare", "octave1", "octave2" , "seqDepths", "seqPlot", "taxPlot")
  
  return(all)
}


##' -------- seqData FLT & cleaning ---------------------
#' for cleaning datasets, flt asvs and samples
#' provide AMP object, seqDepth threshold,  aggregation rank
#' asvPRESence input is in % relative abundance of ASV in each sample, 
#' asvPREValence input is % of samples in dataset,

myFilter <- function(amp, seqDepth= 10e3, asvPRES = 0.01, asvPREV = 1, agg_at = "Species"){
  
  # seqDepth at 10K reads, PRESENCe at 0.01%, PREVALENCE at 3%, aggregate and count Species
  spp_raw <- aggregate_abund(amp$abund, amp$tax, tax_aggregate = agg_at, 
                             calcSums = T, format= "abund")
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
  prvdf =  apply(amp_flt$abund, MARGIN=1, function(x){sum(x>0)}) # calc prev# for each ASV
  amp_flt$abund = amp_flt$abund[ prvdf >= asvPRV , ]#;   #(1)          # filter based on prev#
  amp_flt$tax <- amp_flt$tax[ row.names(amp_flt$abund), ] # cleaning TAX table
  vaf.PREV <- nrow(amp_flt$abund)
  vrem2 <- vaf.PRES - vaf.PREV 
  
  # ========= aggregating abundances at agg_at level (abund table only)
  amp_flt$agg=aggregate_abund(amp_flt$abund, amp_flt$tax, tax_aggregate = agg_at, calcSums = T, format= "abund")
  vaf.spp <- nrow(amp_flt$agg)
  
  
  # ===== stats
  flt_info <- paste(" RAW dataset started with",  sbf, "samples,  represented by",  
                    vbf, "ASVs (",  vbf.spp, "spp). \n",  
                    "- seqDepth filtering at",  seqDepth, "reads per sample,  removed",  srem,
                    "samples, retaining",  saf, "samples. \n", 
                    "- ASV presence filtering at",  asvPRES, "% relative abundance,  removed",  vrem,
                    "ASVs,  retaining",  vaf.PRES, "ASVs. \n", 
                    "- ASV prevalence filtering at",  asvPREV, "% across all samples (<", 
                    asvPRV, "samples), removed additional",  vrem2, "ASVs. \n",  
                    "Finally,",  
                    "FLT dataset is represented by",  saf, "samples with",  vaf.PREV, "ASVs and",  
                    vaf.spp,  "Species.\n")
  
  out <- list(amp_flt, flt_info)
  names(out) <- c("FLTamp", "FLTinfo" )
  return(out)
  
}


##' --------- amp TAX aggregate -----------------------------------
##' corrected Mar 2025
amp_agg_abund <- function(ampObj, agg_at = "Species"){
        taxRanks=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
        cutRanks=taxRanks[1:which(taxRanks==agg_at)]
        
        library(ampvis2)
        agg_abu <- aggregate_abund(ampObj$abund, tax= ampObj$tax, calcSums = T,
                                   tax_aggregate = agg_at, format = "abund")
        
        agg_tax <- ampObj$tax %>% distinct(get(agg_at), .keep_all = T) %>% 
          remove_rownames() %>% mutate(rn = get(agg_at)) %>% column_to_rownames("rn") %>%
          select(all_of(cutRanks))
        
        agg_mrg <- merge(agg_abu, agg_tax, by = 0)  %>% column_to_rownames("Row.names")
        agg_met <- ampObj$metadata
        
        amp_agg <- suppressWarnings( amp_load(agg_mrg, agg_met) ) 
        amp_agg$biom <- agg_mrg
        
        print(paste("input: ", dim(ampObj$abund)[[1]], "taxa"))
        print(paste("output:", dim(agg_abu)[[1]],     "taxa (", agg_at, ") level") )
        
        return(amp_agg)
}


##' --------- FIX species ------------
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
  
  #' fill out the rest of the empty cells
  taxdf %<>% mutate(Genus = ifelse(is.na(Genus), Species, Genus))
  taxdf %<>% mutate(Family = ifelse(is.na(Family), Species, Family) )
  taxdf %<>% mutate(Order = ifelse(is.na(Order), Species, Order) )
  taxdf %<>% mutate(Class = ifelse(is.na(Class), Species, Class) )
  taxdf %<>% mutate(Phylum = ifelse(is.na(Phylum), Species, Phylum) )
  
  
  return(taxdf)
}



##' -------- agglomerate Ganon tsv2tax table -----------------
##' Function takes the tsv2tax table from Ganon and reconciles the lineages 
##' from the different databases, agglomerates count values per sample for 
##' each unique species and chooses a lineage to represent each species.
##' The strain codes from GTdb are removed to accomplis agglomeration

fixGanonTAB <- function(tsv2taxFile) {
  
  ## import and fix up original ganon table
  gantab <- read.table(tsv2taxFile, sep = "\t", header = T, quote = "") %>% as.data.frame() %>% 
    rename_with(~gsub(".tre", "", .x)) %>%                                 ## removes .tre suffix on the sample names
    mutate(SPP = paste0("spp", 1:nrow(.)), .before = lineage) %>%          ## gives temp spp numbers (will be replaced after sums)
    mutate(lineage = gsub("\\'|\\[|\\]|Candidatus " , "", lineage)) ## removes problematic strings producing errors, mismatch or bulky values
  
  print(paste("=== raw ganon tsv2tax table has", dim(gantab)[[1]], "lineages and", dim(gantab[-c(1:2)])[[2]] , "samples"))
  
  ## extract and fix up taxonomy (Species and lineage info)
  gantax <- gantab %>% select(c(lineage, SPP)) %>%                             ## grabs lineage & temp spp numbers
    mutate(lineage = gsub(".__", "", lineage)) %>%                       ## removes the RefSeq leading k__, p__, ... s__ prefix
    separate(lineage, into = taxRanks, sep = "; ", fill = "right") %>%   ## breaks the lineage, filing to the right from lefts
    mutate(Species  = gsub(" sp[0-9]+.*", "", Species)) %>%              ## removes unnecessary GTdb strain codes (suffixes)
    mutate(Kingdom = ifelse(Kingdom == "", NA, Kingdom)) %>%                              ## for some reason, Kingdom needs fixing to NA
    mutate(Species = ifelse(if_all(all_of(taxRanks), is.na), "Unknown", Species)) %>%     ## designates full NA lineages to unknown
    mutate(Species = ifelse(is.na(Species), coalesce(Genus,Family,Order,Class,Phylum,Kingdom), paste(Species) )) %>% ## coalesce SPP to first existing spp
    mutate(Species = ifelse(str_detect(Species, " "), Species, paste(Species, "sp.")))    ## adds " sp." to non-species level names
  # gantax %>% View()
  
  ## agglomerate counts based on fixed taxonomy/uniq Species names
  sumtab <- merge(gantab[-2], gantax[, c("Species", "SPP")], by = c("SPP")) %>%  ## merges abundances and Species names based on tmp spp number
    select(-SPP) %>% summarize(across(everything(), sum), .by = Species)  ## removes the tmp spp numbers & sums all the values based on Species names
  # sumtab %>% View()
  
  ## grab distinct lineages for unique Species 
  out <- list()
  out$tax <- gantax %>% distinct(Species, .keep_all = T) %>%                           ## captures the distinct lineages of each species
    mutate(SPP = paste0("spp", 1:nrow(.)), .after = Species) %>%                 ## gives new/final SPP numbers to each unique taxonomy
    mutate(rn = SPP) %>% column_to_rownames("rn")                                ## sets the SPP numbers to rownames for tax table
  # out$tax %>% View()
  
  ## grabs agglomerated counts
  out$abu <- merge(sumtab, out$tax[, c("Species", "SPP")], by = "Species") %>% 
    column_to_rownames("SPP") %>% select(-Species)
  
  print(paste("=== after SPP-based agglomeration, the abundance table has", dim(out$abu)[[1]], "taxa and", dim(out$abu)[[2]] , "samples"))
  return(out)
}


##' --------- export data ------------------------------------------------

myExport <- function(ampObj, outDIR, fileName = NULL, biom=F, seqs=F, fastaFile = NULL, idString = "TAB"){
  # create output DIR
  if (!dir.exists(outDIR) ){dir.create(outDIR) }
  
  prefix = ifelse(is.null(fileName), "", paste0(fileName, "_") )
  # export COUNTS+TAXranks (no need)
  otutab <- cbind(ampObj$abund, ampObj$tax)
  write.table(otutab,          file.path(paste0(outDIR, "abundtax_", prefix, idString,".txt")), sep="\t",
              col.names = NA, quote=F, na="")
  #export COUNTS (no need)
  # write.table(ampObj$abund,    file.path(paste0(outDIR,"abundata",prefix, idString,".txt")), sep="\t", 
  #             col.names = NA, quote=F, na="") ## we need row names here
  #export TAX (no need)
  write.table(ampObj$tax,       file.path(paste0(outDIR, "taxodata_", prefix, idString,".txt")), sep="\t",
              col.names = NA, quote=F, na="") ## we need row names here
  #export META
  write.table(ampObj$metadata, file.path(paste0(outDIR,"metadata_", prefix, idString,".txt")), sep="\t", 
              col.names = T, quote=F, na="", row.names = F) ## no row names are exported
  
  #export BIOM
  if (biom == T ){
    
    # export LINEAGE (used for biom)
    lineage <- ampObj$tax %>% unite("taxonomy", colnames(ampObj$tax),   sep="; ") # "Taxonomy"
    # write.table(lineage,      file.path(paste0(outDIR, prefix, "lingdata", idString,".txt")), sep="\t",
    #             col.names = NA, quote=F, na="")
    # # export COUNTS+LINEAGE (not used in biom)
    
    lintab <- cbind(ampObj$abund, lineage[1])
    write.table(lintab,       file.path(paste0(outDIR, "abuling_", prefix, idString,".txt")), sep="\t",
                col.names = NA, quote=F, na="")
    
    # library(biomformat) ##fixed based on PRNi's script bellow
    # biomObj<- make_biom(data=as.matrix(ampObj$abund), sample_metadata = ampObj$metadata, 
    #                     observation_metadata = lineage) ## use lineage, not ampObj$tax or lintab
    # write_biom(biomObj,  biom_file = paste0(outDIR, "abundlinJSON.biom"))
    # print("WARNING: please import biom with read_biom(), not with phyloseq, qiime or ampvis biom-import features")
    
    ## Poorani's script for biom   #made work with install.packages("import") && hash printdebug() lines
    myRpath <- paste0("~/Documents/MyApps/RScripts/myRscripts/")
    source(paste0(myRpath, "shared_wME/biomformat.R"))
    source(paste0(myRpath, "shared_wME/ioutils.R"))
    write_biom(dada2biom(otu = t(ampObj$abund), tax = ampObj$tax), biom_file = paste0(outDIR, prefix, "JSON.biom") )
    # rm(write_biom, dada2biom,dada2fasta,dada2taxonomy, dada2text)
    # I dont know why this works but my code doesnt 
  }
  
  
  # export FASTAs seqs
  if (seqs == T) {
    # Check if fastaFile is provided and exists
    if (!is.null(fastaFile) && file.exists(fastaFile)) {
      # Load required library
      library("Biostrings")
      # Imports seqs
      seqs <- readDNAStringSet(fastaFile)
      # Subsets seqs
      seqsFLT <- seqs[names(seqs) %in% ampObj$tax$OTU]
      # Exports filtered seqs
      writeXStringSet(seqsFLT, paste0(outDIR, "seqfasta", idString, ".fasta"))
    } else {
      stop("fastaFile is required if you wish to export SEQs (provide a valid path/to/file/Name.fasta).")
    }
  }
}


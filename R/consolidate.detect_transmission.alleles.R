#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Pancreatic Cancer Microbiome, Oral-Fecal Transmission
#
#Consolidate per-taxon jobs on oral-gut transmission detection from SNV data.
#
#
#=> iterate through taxa
#=> collect data in large frame for further analyses
#=> save and export individual datafiles in R format
#
#2018-10-07
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("plyr", warn.conflicts=F, quietly=T)
library("tidyverse", warn.conflicts=F, quietly=T)
library("effsize", warn.conflicts=F, quietly=T)

options(stringsAsFactors = FALSE)
################################################################################
################################################################################

################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.base <- gsub("R/", "", PARAM$folder.R)
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.metadata <- paste0(PARAM$folder.base, "metadata/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")
PARAM$folder.results_detection <- paste0(PARAM$folder.results, "detect.transmission/")

#Thresholds to filter detection results
PARAM$thresh.min_shared_obs <- 20
################################################################################
################################################################################

################################################################################
################################################################################
#Load sample data
load(paste0(PARAM$folder.metadata, "metaG.oral_fecal.data_sample.Rdata"))

#Load specI abundance data
load(paste0(PARAM$folder.data, "metaG.oral_fecal.data_mOTU_rel.Rdata"))

#Load coverage data
load(paste0(PARAM$folder.data, "metaG.oral_fecal.data_coverage.Rdata"))

#Load taxonomy data
load(paste0(PARAM$folder.data, "metaG.oral_fecal.data_taxonomy.Rdata"))

#Read list of taxa to process
all.taxa <- as.character(scan(paste0(PARAM$folder.parameters, "metaG.oral_fecal.tax_ids_abd.txt")))
################################################################################
################################################################################

################################################################################
################################################################################
#Preallocate data.frame to collect results
data.detect_transmission <- data.frame()
data.taxa.transmission <- data.frame()

#Define function to Z transform a distribution
z.transform <- function(x) {(x - mean(x, na.rm=T)) / sd(x, na.rm=T)}

#Get list of unique subjects
subjects <- as.character(unique(data.sample$subject))

#Iterate through taxa and process/collect results
for (tax.id in all.taxa){
  #Get current taxon name
  taxon <- data.taxonomy[data.taxonomy$repID == tax.id, "refOTUs"]
  print.taxon <- gsub(" ", "_", taxon)
  print.taxon <- gsub("/", "_", print.taxon)
  
  #Get current file names
  file.data_transmission <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".detect_transmission.RData")
  file.data_bg <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".background_scores.RData")
  
  #Skip if files for taxon don't exist
  if (! file.exists(file.data_bg)) {next()}
  if (file.info(file.data_bg)$size == 0) {next()}
  
  #Load data
  load(file.data_transmission)
  load(file.data_bg)
  
  #Remove all intra-individual pairs with too few shared observations
  data.s_t.oi <- data.s_t.oi[data.s_t.oi$n.obs.shared > PARAM$thresh.min_shared_obs & !is.na(data.s_t.oi$n.obs.shared), ]
  #Skip if no data remains
  if (nrow(data.s_t.oi) == 0) {next()}
  
  #
  #
  #Skip if data is still old
  if (ncol(data.s_t.oi) > 20) {next()}
  #
  #
  
  #Skip if too few points exist (in full bg) to properly estimate background
  if (nrow(data.bg) == 0) {next()}
  if (median(data.bg$n.comparisons) < 100) {next()}
    
  ###################################
  #Skip if too few subjects were score-able
  if (sum(!is.na(data.s_t.oi$transmission.score)) <= 5) {next()}
  ###################################
  
  ###################################
  #Get current data
  curr.frame <- cbind(
    data.frame(
      Tax_ID=rep.int(tax.id, nrow(data.s_t.oi)),
      Taxon=rep.int(taxon, nrow(data.s_t.oi)),
      stringsAsFactors=F
    ),
    data.s_t.oi
  )
  curr.frame <- curr.frame[!is.na(curr.frame$transmission.score), ]
  ###################################
  
  ###################################
  #Calculate per-subject p values
  #=> is a subject's transmission score significant relative to the subject-specific background?
  curr.frame[, "transmission.score.p.adj"] <- p.adjust(1 - pnorm(curr.frame$transmission.score), method="BH")
  ###################################
  
  ###################################
  #Annotate current taxon
  curr.data.taxa <- cbind(
    data.taxonomy[data.taxonomy$repID == tax.id, ],
    data.frame(
      n.indiv.scored=nrow(curr.frame),
      transmission.score.mean = mean(curr.frame$transmission.score, na.rm=T),
      transmission.score.sd = sd(curr.frame$transmission.score, na.rm=T),
      transmission.score.median = median(curr.frame$transmission.score, na.rm=T),
      transmission.score.cohen_d = cohen.d(curr.frame$transmission.score, unlist(S.bg))$estimate,
      transmission.score.cliff_delta = cliff.delta(curr.frame$transmission.score, unlist(S.bg))$estimate,
      transmission.score.wilcox_p.raw = ifelse(all(is.na(curr.frame$transmission.score)), NA, wilcox.test(curr.frame$transmission.score, unlist(S.bg), paired=F, alternative="greater")$p.value)
    )
  )
  data.taxa.transmission <- rbind(data.taxa.transmission, curr.data.taxa)
  ###################################
  
  ###################################
  #Append current data to filtered results collector frame
  data.detect_transmission <- rbind(data.detect_transmission, curr.frame)
  ###################################
  
  ###################################
  #Report
  writeLines(paste(date(), "=> done with", taxon))
  ###################################
}
################################################################################
################################################################################


################################################################################
################################################################################
#Post-process transmission detection data
################################################################################
################################################################################
#Remove all intra-individual comparisons with too few shared observations
data.detect_transmission <- data.detect_transmission[data.detect_transmission$n.obs.shared >= PARAM$thresh.min_shared_obs, ]

###################################
#Correct p.values on transmission Z scores and generate genus-level plots
data.taxa.transmission$transmission.score.wilcox.p <- p.adjust(data.taxa.transmission$transmission.score.wilcox_p.raw, method="BH")

#Set significance levels as factor
data.taxa.transmission$transmission.score.wilcox.sig <- "FALSE"
curr.no_na <- !is.na(data.taxa.transmission$transmission.score.wilcox.p)
data.taxa.transmission[data.taxa.transmission$transmission.score.wilcox.p < 0.05 & curr.no_na, "transmission.score.wilcox.sig"] <- "<0.05"
data.taxa.transmission[data.taxa.transmission$transmission.score.wilcox.p < 0.01 & curr.no_na, "transmission.score.wilcox.sig"] <- "<0.01"
data.taxa.transmission[data.taxa.transmission$transmission.score.wilcox.p < 0.001 & curr.no_na, "transmission.score.wilcox.sig"] <- "<0.001"

#Assign logical vector "transmitter"
data.taxa.transmission$transmitter <- F
data.taxa.transmission[data.taxa.transmission$transmission.score.wilcox.sig != "FALSE", "transmitter"] <- T
###################################

###################################
#Reorder taxa factors in results frames
#taxa.order.sci_name <- as.character(data.taxa$Scientific_Name[data.taxa$Scientific_Name %in% data.detect_transmission$Taxon])
#data.detect_transmission$Taxon <- factor(data.detect_transmission$Taxon, levels=taxa.order.sci_name)
###################################
################################################################################
################################################################################

#Save data
save(data.detect_transmission, data.taxa.transmission, file=paste0(PARAM$folder.results, "data.detect.transmission.RData"))
################################################################################
################################################################################



#q(save="no")






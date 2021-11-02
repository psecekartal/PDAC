#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Pancreatic Cancer Microbiome, Oral-Fecal Transmission
#
#Detect putative oral or fecal strain associations to disease and confounders
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load SNV tables per taxon
#=> calculate pairwise SNV distance between samples
#=> cluster by SNV distances and test for putative sub-species structure
#=> test whether SNV distances are associated to disease status and/or confounders
#=> save and export individual datafiles in R format
#
#2018-10-12
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
################################################################################
################################################################################


################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.base <- paste0("/g/bork1/tschmidt/pancreatic_cancer_microbiome/")
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.metadata <- paste0(PARAM$folder.base, "metadata/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Set parameters
#=> minimum allele prevalence in oral and/or fecal samples
PARAM$min.allele_prevalece <- 0.1
#=> minimum observed positions per sample
PARAM$min.sample_size <- 500
#=> minimum required shared alleles for comparisons
PARAM$min.shared_obs <- 100

#Read current taxon of interest from command line
args <- commandArgs(TRUE)
tax.id <- args[2]

#Report
writeLines(paste(date(), "=> Processing", tax.id))
################################################################################
################################################################################


################################################################################
################################################################################
#Load allele data
file.allele_data <- paste0(PARAM$folder.SNV_data, tax.id, ".allele_data.filtered.RData")
if (! file.exists(file.allele_data)) {writeLines(paste(tax.id, "file does not exist:", file.allele_data)); q("no")}
load(file.allele_data)

#Load sample data
load(paste0(PARAM$folder.metadata, "metaG.oral_fecal.data_sample.Rdata"))
samples <- intersect(colnames(mat.allele), rownames(data.sample))
#Remove biological replicates for 3 subjects, for the sake of simplicity
#samples <- samples[! samples %in% c("MMPC33706732OR", "MMPC11113549ST", "MMPC92345363OR", "MMPC48555115OR")]
n.samples <- length(samples)
data.sample <- data.sample[samples, ]

#Load taxonomy data
load(paste0(PARAM$folder.data, "metaG.oral_fecal.data_taxonomy.Rdata"))

#Get current taxon name
taxon <- data.taxonomy[data.taxonomy$repID == tax.id, "refOTUs"]
print.taxon <- gsub(" ", "_", taxon)
print.taxon <- gsub("/", "_", print.taxon)
################################################################################
################################################################################


################################################################################
################################################################################
#Filter allele matrix to include only SNVs satisfying a minimum prevalence criterion
mat.allele <- mat.allele[, colnames(mat.allele) %in% samples]

#Get logcial matrices of observed positions and observed alleles
obs.allele <- mat.allele > 0
inc.allele <- mat.allele > 10^-12

#Get current habitat vector
habitats <- data.sample[colnames(mat.allele), "body_site"]
subjects <- data.sample[colnames(mat.allele), "subject_id"]

#Get indices of positions satisfying minimum SNV prevalece criterion
idx.inc.oral <- rowSums(inc.allele[, habitats == "saliva", drop=FALSE]) > (PARAM$min.allele_prevalece * sum(habitats == "saliva"))
idx.inc.gut <- rowSums(inc.allele[, habitats == "feces", drop=FALSE]) > (PARAM$min.allele_prevalece * sum(habitats == "feces"))
idx.keep <- idx.inc.oral | idx.inc.gut

if (sum(idx.keep) < 1000) {writeLines(paste(date(), "=> Too few SNV positions at minimum prevalence:", tax.id, taxon)); q(save="no")}

#Apply filter
obs.allele <- obs.allele[idx.keep, ]
inc.allele <- inc.allele[idx.keep, ]
mat.allele <- mat.allele[idx.keep, ]

#Filter for samples with minimum observations
idx.samples.min_size <- colSums(obs.allele) > PARAM$min.sample_size

if (sum(idx.samples.min_size) < 20) {writeLines(paste(date(), "=> Too few samples at minimum required SNV size:", tax.id, taxon)); q(save="no")}

#Apply filter
obs.allele <- obs.allele[, idx.samples.min_size]
inc.allele <- inc.allele[, idx.samples.min_size]
mat.allele <- mat.allele[, idx.samples.min_size]
data.sample <- data.sample[idx.samples.min_size, ]
habitats <- habitats[idx.samples.min_size]
subjects <- subjects[idx.samples.min_size]

#Re-zero mat.allele
mat.allele[mat.allele == 10^-12] <- 0

#Get current total allele count
n.allele <- nrow(obs.allele)
################################################################################
################################################################################


################################################################################
################################################################################
#Get all sample pairs
curr.sample_pairs <- t(as.matrix(combn(colnames(obs.allele), 2)))
curr.pairs.st <- cbind(as.character(data.sample[curr.sample_pairs[,1], "subject_id"]), as.character(data.sample[curr.sample_pairs[,2], "subject_id"]))
curr.pairs.habitats <- cbind(as.character(data.sample[curr.sample_pairs[,1], "body_site"]), as.character(data.sample[curr.sample_pairs[,2], "body_site"]))
curr.pairs.intra <- curr.sample_pairs[curr.pairs.st[,1] == curr.pairs.st[,2], ]
###################################

###################################
#Define functions to compute pairwise SNV distances
get.dist <- function(pair) {
  #Get current indices of shared observations
  idx.obs <- rowSums(obs.allele[, pair]) == 2
  if (sum(idx.obs) < PARAM$min.shared_obs) {return(c(NA, NA))}
  
  #Scale mismatches by total comparisons
  tmp.sums <- rowSums(inc.allele[idx.obs, pair])
  dist.inc <- sum(tmp.sums == 1) / sum(idx.obs)
  
  curr.mat.allele <- mat.allele[idx.obs, pair]
  dist.freq <- sum(abs(curr.mat.allele[, 1] - curr.mat.allele[, 2])) / sum(idx.obs)
  
  c(dist.inc, dist.freq)
}
###################################

###################################
#Compute all pairwise SNV distances
#=> this will take time!
curr.dist <- t(apply(curr.sample_pairs, 1, get.dist))
colnames(curr.dist) <- c("dist.inc", "dist.freq")
###################################

save(curr.dist, curr.sample_pairs, file=paste0(PARAM$folder.results, "allele.distances/", tax.id, ".", print.taxon,  ".allele_distances.RData"))
################################################################################
################################################################################



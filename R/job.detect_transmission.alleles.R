#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Pancreatic Cancer Microbiome, Oral-Fecal Transmission
#
#Detect putative oral-gut transmission from SNV data.
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load SNV tables per taxon
#=> calculate SNV frequencies in oral and gut samples, per position
#=> for each individual, calculate deviation of observed SNV patterns from bg
#=> save and export individual datafiles in R format
#
#2018-10-07
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
#=> minimum required shared alleles for comparisons
PARAM$min.shared_obs <- 20
#=> maximum number of background comparisons (beyond which distributions are sampled)
PARAM$max.background_size <- 1000

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
################################################################################
################################################################################

################################################################################
################################################################################
#Filter down to include only positions which are a SNV in ≥1 saliva and ≥1 gut sample
mat.allele <- mat.allele[data.allele$covered.oral & data.allele$covered.gut, colnames(mat.allele) %in% samples]
data.allele <- data.allele[data.allele$covered.oral & data.allele$covered.gut, ]

#Get current habitat vector
habitats <- data.sample[colnames(mat.allele), "body_site"]
subjects <- data.sample[colnames(mat.allele), "subject_id"]

#Get current taxon name
taxon <- data.taxonomy[data.taxonomy$repID == tax.id, "refOTUs"]
print.taxon <- gsub(" ", "_", taxon)
print.taxon <- gsub("/", "_", print.taxon)

###################################
#Identify positions of interest
#=> minimum number of observations in both saliva and gut (≥10 samples each)
#=> allele in at least one oral and one gut sample (filtered already above)
obs.allele <- mat.allele > 0
obs.oral <- obs.allele[, habitats == "saliva"]
obs.gut <- obs.allele[, habitats == "feces"]

#Get an allele incidence matrix
#=> TRUE if allele was observed, FALSE if not (ref or "-1")
inc.allele <- mat.allele > 10^-12
inc.oral <- inc.allele[, habitats == "saliva"]
inc.gut <- inc.allele[, habitats == "feces"]

#Exit if positions are distributed across habitats too unevenly
if (nrow(obs.oral) < 20) {writeLines(paste(date(), "=> Too few oral positions:", tax.id)); q(save="no")}
if (nrow(obs.gut) < 20) {writeLines(paste(date(), "=> Too few gut positions:", tax.id)); q(save="no")}

#Get raw observation frequencies per habitat
obs.oral.sum <- rowSums(obs.oral)
obs.gut.sum <- rowSums(obs.gut)

#Get raw incidence frequencies per habitat
inc.oral.sum <- rowSums(inc.oral)
inc.gut.sum <- rowSums(inc.gut)

#Require that an allele...
#=> is observed in at least one oral and one gut sample (otherwise, would be trivial)
#=> is not observed in *all* oral or *all* gut samples (otherwise, would be trivial, as well)
#=> has position coverage in at least 5 oral and 5 gut samples
keep.alleles <- obs.oral.sum >= 5 & obs.gut.sum >= 5 & inc.oral.sum >= 1 & inc.gut.sum >= 1 & inc.oral.sum < ncol(inc.oral) & inc.gut.sum < ncol(inc.gut)

#Exit if too few alleles remain
if (sum(keep.alleles) < 100) {writeLines(paste(date(), "=> Too few observed alleles:", tax.id)); q(save="no")}

#Filter down to relevant positions
mat.allele <- mat.allele[keep.alleles, ]
###################################

###################################
#Re-turn mat.allele into an incidence matrix (after filtering)
#=> TRUE if position was observed (as current allele or other), FALSE if not ("-1")
obs.allele <- mat.allele > 0
obs.oral <- obs.allele[, habitats == "saliva"]
obs.gut <- obs.allele[, habitats == "feces"]

#Get an allele incidence matrix
#=> TRUE if allele was observed, FALSE if not (ref or "-1")
inc.allele <- mat.allele > 10^-12
inc.oral <- inc.allele[, habitats == "saliva"]
inc.gut <- inc.allele[, habitats == "feces"]

#Tidy up
rm(mat.allele); rm(data.allele); invisible(gc())

#Get current total allele count
n.allele <- nrow(obs.allele)

#Get current habitat vector
habitats <- data.sample[colnames(obs.allele), "body_site"]
subjects <- data.sample[colnames(obs.allele), "subject_id"]

#Get the number of remaining eligible subject.timepoints with paired saliva-gut observations
curr.sample_pairs <- as.matrix(expand.grid(colnames(obs.oral), colnames(obs.gut)))
curr.pairs.st <- cbind(as.character(data.sample[curr.sample_pairs[,1], "subject_id"]), as.character(data.sample[curr.sample_pairs[,2], "subject_id"]))
curr.pairs.intra <- curr.sample_pairs[curr.pairs.st[,1] == curr.pairs.st[,2], ]

#Skip current taxon if there are no paired subject.timepoints left
if (nrow(curr.pairs.intra) < 2) {writeLines(paste(date(), tax.id, taxon, "=> Too few intra-individual paired observations.")); q("no")}

#Get list of subjects with paired intra observations
curr.pairs.subject <- as.character(data.sample[curr.pairs.intra[,1], "subject_id"])
###################################

###################################
#Collapse timepoints in allele observation and incidence tables
#=> combine information from multiple timepoints of an individual
#=> retain any allele that was *ever* observed for a given subject
#=> the updated tables are alleles x subject, where subjects are collapsed longitudinally
#=> additionally, populate a data.subject frame for metadata
###################################
#Preallocate
subjects.with_oral <- unique(as.character(subjects[habitats == "saliva"]))
subjects.with_gut <- unique(as.character(subjects[habitats == "feces"]))
obs.allele.subject.oral <- obs.allele[, habitats == "saliva"]
inc.allele.subject.oral <- inc.allele[, habitats == "saliva"]
colnames(obs.allele.subject.oral) <- colnames(inc.allele.subject.oral) <- data.sample[colnames(obs.allele.subject.oral), "subject_id"]
obs.allele.subject.gut <- obs.allele[, habitats == "feces"]
inc.allele.subject.gut <- inc.allele[, habitats == "feces"]
colnames(obs.allele.subject.gut) <- colnames(inc.allele.subject.gut) <- data.sample[colnames(obs.allele.subject.gut), "subject_id"]

#Iterate through subjects
data.subject <- data.frame()
for (subj in union(subjects.with_oral, subjects.with_gut)) {
  #Collapse oral samples
  if (subj %in% subjects.with_oral) {available.oral <- TRUE} else {available.oral <- FALSE}
  
  #Collapse gut samples
  if (subj %in% subjects.with_gut) {available.gut <- TRUE} else {available.gut <- FALSE}
  
  #Collect info on current subject
  data.subject <- rbind(data.subject, data.frame(
    subject = subj,
    available.oral = available.oral,
    available.gut = available.gut
  ))
}
rownames(data.subject) <- as.character(data.subject$subject)

#Further filter alleles to exlucde trivial cases
#=> do not consider alleles which were observed in *all* oral or *all* gut samples
keep.alleles <- (rowSums(inc.allele.subject.oral) < rowSums(obs.allele.subject.oral)) & (rowSums(inc.allele.subject.gut) < rowSums(obs.allele.subject.gut))
if (sum(keep.alleles) < 100) {writeLines(paste(date(), tax.id, taxon, "=> Too few alleles remaining after filtering.")); q(save="no")}

#Apply filter
obs.allele.subject.oral <- obs.allele.subject.oral[keep.alleles,]
inc.allele.subject.oral <- inc.allele.subject.oral[keep.alleles,]
obs.allele.subject.gut <- obs.allele.subject.gut[keep.alleles,]
inc.allele.subject.gut <- inc.allele.subject.gut[keep.alleles,]

obs.allele <- obs.allele[keep.alleles,]
obs.oral <- obs.oral[keep.alleles,]
obs.gut <- obs.gut[keep.alleles,]

inc.allele <- inc.allele[keep.alleles,]
inc.oral <- inc.oral[keep.alleles,]
inc.gut <- inc.gut[keep.alleles,]

#Get current total allele count
n.allele <- nrow(obs.allele)
###################################

###################################
#Calculate subject-specific allele frequency backgrounds
bg.allele.subject.oral <- rowSums(inc.allele.subject.oral) / rowSums(obs.allele.subject.oral)
bg.allele.subject.gut <- rowSums(inc.allele.subject.gut) / rowSums(obs.allele.subject.gut)
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate background probabilities of observing alleles
#=> for each pair,
#=> sum log probabilities of observation for shared alleles (S_shared ≤ 0)
#=> sum -log probabilities of observation for non-shared alleles (S_nonshared ≥ 0)
#=> sum log probabilities for least probable shared obs per pos (S_max ≤ 0)
#=> calculate score by scaling S = (S_shared - S_obs) / S_max
###################################
#=> p(OR & ST)
bg.log_p.full.1_1 <- log10(bg.allele.subject.oral * bg.allele.subject.gut)
#=> p(OR & !ST)
bg.log_p.full.1_0 <- -log10(bg.allele.subject.oral * (1 - bg.allele.subject.gut))
#=> p(!OR & ST)
bg.log_p.full.0_1 <- -log10((1 - bg.allele.subject.oral) * bg.allele.subject.gut)
#=> p(!OR & !ST)
bg.log_p.full.0_0 <- log10((1 - bg.allele.subject.oral) * (1 - bg.allele.subject.gut))
###################################

###################################
#Calculate maximum theoretical scores per subject
#=> the maximum attainable score vs background is the min(log(p)) for the shared observations (1_1 and 0_0) per each allele
S.max_theo.bg <- pmin(bg.log_p.full.1_1, bg.log_p.full.0_0)
################################################################################
################################################################################

################################################################################
################################################################################
#Define current subjects & timepoints of interest
curr.s_t.oi <- curr.pairs.subject

#Reduce data.timepoint to only subjects of interest
data.s_t.oi <- data.subject[curr.s_t.oi, ]
#Define fields to add
data.s_t.fields <- c(
  "status",
  "sample.oral", "sample.gut",
  #Total number of alleles w/ observation in individual
  "n.obs", "n.allele",
  #Observed allele positions per body site
  "n.obs.oral", "n.allele.oral", "n.obs.gut", "n.allele.gut",
  #Number of shared observations and shared allele observations
  "n.obs.shared", "n.allele.shared",
  #Raw Jaccard index of shared allele observations
  "jac.shared",
  #Raw scores of current observation against background
  "S.agree", "S.disagree", "S.max_theo",
  #Probability scores against full background, across all pairs 
  "S",
  #Z scores against full and cohort-specific background
  "transmission.score"
)
#Preallocate additional columns
data.s_t.oi <- cbind(data.s_t.oi, as.data.frame(matrix(nrow=nrow(data.s_t.oi), ncol=length(data.s_t.fields), dimnames=list(rownames(data.s_t.oi), data.s_t.fields))))
###################################

###################################
#Get habitat, subject.timepoint, subject, family & background cohort for each pair of samples
#=> "cp" for "current pairs"
cp.habitat <- cbind(as.character(data.sample[curr.sample_pairs[,1], "body_site"]), as.character(data.sample[curr.sample_pairs[,2], "body_site"]))
cp.subject <- cbind(as.character(data.sample[curr.sample_pairs[,1], "subject_id"]), as.character(data.sample[curr.sample_pairs[,2], "subject_id"]))
cp.status <- cbind(as.character(data.sample[curr.sample_pairs[,1], "status"]), as.character(data.sample[curr.sample_pairs[,2], "status"]))

#Get a list of "eligible" pairs
#=> pairs which satisfy the "minimum shared observations" criterion
#=> this can take a long time
#cp.min_obs <- apply(curr.sample_pairs[1:20, ], 1, function(pair) {sum(rowSums(obs.allele[, pair]) == 2) > PARAM$min.shared_obs})
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate pairwise scores and theoretical maximum scores
##=> per background
#=> for all oral-gut pairs w/in background
#=> using subject-specific allele frequency profiles
################################################################################
################################################################################
#Score pairwise observations between all pairs of samples, for different backgrounds
S.bg <- S.bg.Z <- list()
data.transmission <- data.bg <- data.frame()

#Define function to score paired observations (OR-ST)
get.S <- function(pair, s_t, keep=rep.int(T, length(S.max_theo.bg))) {
  #Get indices of alleles to consider
  idx.obs <- keep & rowSums(obs.allele[, pair]) == 2
  
  #Return NA if there are too few shared observations
  if (sum(idx.obs) < PARAM$min.shared_obs) {return(as.numeric(rep.int(NA, 4)))}
  
  #Subset data to current pair
  curr.inc.allele <- inc.allele[, pair]
  
  #Get indices for direct comparisons
  idx.1_1 <- idx.obs & (rowSums(curr.inc.allele) == 2)
  idx.1_0 <- idx.obs & curr.inc.allele[, 1] & (!curr.inc.allele[, 2])
  idx.0_1 <- idx.obs & curr.inc.allele[, 2] & (!curr.inc.allele[, 1])
  idx.0_0 <- idx.obs & (!curr.inc.allele[, 1]) & (!curr.inc.allele[, 2])
  
  #Compute raw scores
  S.agree <- (sum(bg.log_p.full.1_1[idx.1_1], na.rm=T) + sum(bg.log_p.full.0_0[idx.0_0], na.rm=T))
  S.disagree <- (sum(bg.log_p.full.1_0[idx.1_0], na.rm=T) + sum(bg.log_p.full.0_1[idx.0_1], na.rm=T))
  S.max_theo <- sum(S.max_theo.bg[idx.obs])
  
  #Return
  c(
    S.agree,
    S.disagree,
    S.max_theo,
    #Get transmission score: scale by maximally unlikely agreement
    (S.agree - S.disagree) / S.max_theo
  )
}

#Define function to Z transform a distribution
z.transform <- function(x) {(x - mean(x, na.rm=T)) / sd(x, na.rm=T)}

#Iterate through subject.timepoints and calculate transmission scores
#=> using the subject-specific background allele frequencies, both within and across cohorts
for (s_t in curr.s_t.oi) {
  ###################################
  #Get current subject ID
  subj <- s_t
  ###################################
  
  ###################################
  #Select alleles which have non-trivial observations in current cohort
  #=> observed in ≥1 OR and ≥1 ST sample w/in cohort
  keep.alleles <- is.finite(S.max_theo.bg)
  
  #Skip current subject.timepoint if too few alleles with finite background frequencies are left
  if (sum(keep.alleles) < 100) {next()}
  ###################################
  
  ###################################
  #Calculate scores across all pairs (THIS TAKES TIME!)
  #=> subsample to a (lower) maximum number of comparisons
  #
  #Obtain background score distributions
  #=> pre-require minimum shared observations
  #=> ignore intra-subject comparisons
  #=> but consider all comparisons to current subject.timepoint
  ###################################
  #=> against a full background
  idx.pairs <- cp.subject[,1] != cp.subject[,2]
  if (sum(idx.pairs) > PARAM$max.background_size) {
    curr.S <- t(apply(curr.sample_pairs[sample(which(idx.pairs), size=PARAM$max.background_size, replace=F), ], 1, get.S, s_t=s_t, keep=keep.alleles))
  } else {
    curr.S <- t(apply(curr.sample_pairs[idx.pairs, ], 1, get.S, s_t=s_t, keep=keep.alleles))
  }
  colnames(curr.S) <- c("S.agree", "S.disagree", "S.max_theo", "S")
  S.bg[[s_t]] <- na.omit(curr.S[, "S"])
  ###################################
  #Get current intra-individual (intra-timepoint) comparison
  ###################################
  collect.S.intra <- get.S(curr.sample_pairs[cp.subject[,1] == s_t & cp.subject[,2] == s_t, ], s_t=s_t, keep=keep.alleles)
  names(collect.S.intra) <- c("S.agree", "S.disagree", "S.max_theo", "S")
  ###################################
  
  ###################################
  #Re-scale (Z-transform) background distributions
  S.bg.Z[[s_t]] <- z.transform(S.bg[[s_t]])
  #Z-transform current intra-individual comparisons
  S.z.intra <- (collect.S.intra["S"] - mean(S.bg[[s_t]], na.rm=T)) / sd(S.bg[[s_t]], na.rm=T)
  ###################################
  
  ###################################
  #Store data on background distributions
  ###################################
  #Full background
  data.bg <- rbind(data.bg, data.frame(
    s_t = s_t, n.comparisons = sum(idx.pairs), S.mean = mean(S.bg[[s_t]], na.rm=T), S.sd = sd(S.bg[[s_t]], na.rm=T)
  ))
  ###################################
  
  ###################################
  #Store data on intra-individual transmission
  #=> metadata (coverages, etc.)
  ###################################
  #Get current oral and gut sample id
  intra.oral <- rownames(data.sample)[data.sample$subject_id == s_t & data.sample$body_site == "saliva"][1]
  intra.gut <- rownames(data.sample)[data.sample$subject_id == s_t & data.sample$body_site == "feces"][1]
  ###################################
  #Update data.s_t.oi
  ###################################
  data.s_t.oi[s_t, "status"] <- as.character(data.sample[intra.oral, "status"])
  data.s_t.oi[s_t, "sample.oral"] <-intra.oral
  data.s_t.oi[s_t, "sample.gut"] <-intra.gut
  #Allele observation and count statistics
  data.s_t.oi[s_t, "n.obs"] <- sum(rowSums(obs.allele[, c(intra.oral, intra.gut)]) >= 1)
  data.s_t.oi[s_t, "n.allele"] <- sum(rowSums(inc.allele[, c(intra.oral, intra.gut)]) >= 1)
  data.s_t.oi[s_t, "n.obs.oral"] <- sum(obs.allele[, intra.oral])
  data.s_t.oi[s_t, "n.allele.oral"] <- sum(inc.allele[, intra.oral])
  data.s_t.oi[s_t, "n.obs.gut"] <- sum(obs.allele[, intra.gut])
  data.s_t.oi[s_t, "n.allele.gut"] <- sum(inc.allele[, intra.gut])
  data.s_t.oi[s_t, "n.obs.shared"] <- sum(rowSums(obs.allele[, c(intra.oral, intra.gut)]) == 2)
  data.s_t.oi[s_t, "n.allele.shared"] <- sum(rowSums(inc.allele[, c(intra.oral, intra.gut)]) == 2)
  #Jaccard overlap of observed alleles
  data.s_t.oi[s_t, "jac.shared"] <- data.s_t.oi[s_t, "n.allele.shared"] / data.s_t.oi[s_t, "n.obs.shared"]
  #Scores on full background
  data.s_t.oi[s_t, "S.agree"] <- collect.S.intra["S.agree"]
  data.s_t.oi[s_t, "S.disagree"] <- collect.S.intra["S.disagree"]
  data.s_t.oi[s_t, "S.max_theo"] <- collect.S.intra["S.max_theo"]
  data.s_t.oi[s_t, "S"] <- collect.S.intra["S"]
  #Z-transformed transmission score
  data.s_t.oi[s_t, "transmission.score"] <- S.z.intra
  ###################################
}

###################################
#Save current data
save(data.bg, data.s_t.oi, file=paste0(PARAM$folder.results, "detect.transmission/", tax.id, ".", print.taxon,  ".detect_transmission.RData"))
save(S.bg, file=paste0(PARAM$folder.results, "detect.transmission/", tax.id, ".", print.taxon,  ".background_scores.RData"))
save(S.bg.Z, file=paste0(PARAM$folder.results, "detect.transmission/", tax.id, ".", print.taxon,  ".background_scores.z_transformed.RData"))
###################################

#Report
writeLines(paste(date(), "=> Done with", tax.id, taxon))
################################################################################
################################################################################



q(save="no")



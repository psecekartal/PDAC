#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Pancreatic Cancer Microbiom, Oral-Fecal Transmission
#
#Prepare SNV data for further analyses
#
#Per each tax ID, as read from command line...
#=> load sample metadata
#=> load SNV table
#=> filter SNV table
#=> save and export individual datafiles in R format
#
#2018-10-05
#sebastian.schmidt@embl.de
################################################################################


################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("gtools", warn.conflicts=F, quietly=T)
################################################################################
################################################################################


################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.base <- paste0("/g/bork1/tschmidt/pancreatic_cancer_microbiome/")
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.metadata <- paste0(PARAM$folder.base, "metadata/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "metaSNV_output/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Read current taxon of interest from command line
args <- commandArgs(TRUE)
tax.id <- args[2]
################################################################################
################################################################################


################################################################################
################################################################################
#Load sample metadata
load(paste0(PARAM$folder.metadata, "metaG.oral_fecal.data_sample.Rdata"))
samples <- rownames(data.sample)
n.samples <- length(samples)

#Load specI abundance data
load(paste0(PARAM$folder.data, "metaG.oral_fecal.data_mOTU_rel.Rdata"))

#Load taxonomy data
load(paste0(PARAM$folder.data, "metaG.oral_fecal.data_taxonomy.Rdata"))
################################################################################
################################################################################


################################################################################
################################################################################
#Load and process SNV count tables per taxon
#=> process population SNVs and individual SNVs separately
#=> skip all taxa which have not been observed in at least one oral and gut sample
################################################################################
################################################################################
#Get current taxon name
curr.taxon <- data.taxonomy$refOTUs[data.taxonomy$repID == tax.id]
print.taxon <- gsub(" ", "_", curr.taxon)
print.taxon <- gsub("/", "_", print.taxon)

#Get current file names for current taxon
curr.file.pop <- paste0(PARAM$folder.SNV_data, "filtered/pop/", tax.id, ".filtered.freq.gz")
curr.file.indiv <- paste0(PARAM$folder.SNV_data, "filtered/indiv/", tax.id, ".filtered.freq.gz")

#Report
writeLines(paste(date(), "=> Processing", tax.id, curr.taxon))
###################################

###################################
#Read data
###################################
#Check if population-level SNV file exists
if (! file.exists(curr.file.pop)) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> No population-level SNV file available.\n")); q(save="no")}

#Read data
#=> first, read only the first few lines to check if data comes from ≥1 habitat
tmp.dat <- read.table(curr.file.pop, sep="\t", header=T, nrows=2, row.names=1)

#Process sample names
tmp.samples <- colnames(tmp.dat)
curr.samples <- curr.samples <- gsub("[.]", "_", tmp.samples)
curr.samples <- gsub("_ULRepGenomesv11_unique_sorted_bam", "", curr.samples)
curr.samples <- curr.samples[curr.samples %in% rownames(data.sample)]

#Check if current taxon was observed in more than one habitat
curr.habitats <- as.character(data.sample[curr.samples, "body_site"])
n.hab <- length(unique(curr.habitats))

#Skip current taxon if only samples from one habitat are covered
if (n.hab == 1) {
  writeLines(paste0("Skip: only one habitat covered (", unique(curr.habitats), ")\n"))
  q(save="no")
}

#Report
writeLines(paste0(date(), " => ", tax.id, " ", curr.taxon, " => Both habitats covered: oral cavity (", length(which(curr.habitats == "saliva")), "), stool (", length(which(curr.habitats == "feces")), ")."))
###################################

###################################
#Read data
#=> for real, this time
tmp.dat <- read.table(curr.file.pop, sep="\t", header=T, row.names=1)

#Report
writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", nrow(tmp.dat), "raw pop-level SNV positions."))
###################################

###################################
#Check whether individual SNPs exist for current taxon
if (file.exists(curr.file.indiv)) {
  #Read individual SNP data
  tmp.dat.indiv <- read.table(curr.file.indiv, sep="\t", header=T, row.names=1)
  tmp.mat.pos.indiv <- as.matrix(tmp.dat.indiv)
  #Report
  writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", nrow(tmp.dat.indiv), "raw indiv-level SNP positions."))
  
  #Concatenate
  tmp.tmp.dat <- rbind(tmp.dat, tmp.dat.indiv)
  #Re-order
  new.order <- mixedorder(rownames(tmp.tmp.dat))
  tmp.dat <- tmp.tmp.dat[new.order, ]
  rm(tmp.tmp.dat); invisible(gc())
}

#Skip if too few positions were encountered
if (nrow(tmp.dat) < 10) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Skip: too few positions.\n")); q(save="no")}
###################################

###################################
#Process sample names
tmp.samples <- colnames(tmp.dat)
curr.samples <- curr.samples <- gsub("[.]", "_", tmp.samples)
curr.samples <- gsub("_.+", "", curr.samples)
colnames(tmp.dat) <- curr.samples

#Filter to include only characterised samples (with metadata)
curr.samples <- curr.samples[curr.samples %in% rownames(data.sample)]
tmp.dat <- tmp.dat[, curr.samples]

#"Vertically" filter out samples which have no information (-1) across all positions
#=> moreover, remove samples that have previously been dropped based on abundance/coverage data
#=> this is also a sanity check
tmp.mat.pos <- as.matrix(tmp.dat)
keep.samples <- (colSums(tmp.mat.pos) > (-nrow(tmp.mat.pos))) & (curr.samples %in% rownames(data.sample))
curr.samples <- curr.samples[keep.samples]
tmp.mat.pos <- tmp.mat.pos[, keep.samples]

#Report
writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", nrow(tmp.dat), "SNP positions across", length(curr.samples), "samples."))

if (length(curr.samples) < 2) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Skip: too few samples covered appropriately.\n")); q(save="no")}

#Report if a sample ID does not match the original table
if (length(which(! curr.samples %in% rownames(data.sample))) > 0) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Unknown samples:", curr.samples[! curr.samples %in% rownames(data.sample)]))}

#Get current habitats (from fixed sample names)
curr.habitats <- as.character(data.sample[curr.samples, "body_site"])
n.hab <- length(unique(curr.habitats))
###################################

###################################
#Extract current information per position
tmp.pos_info <- do.call(rbind, strsplit(rownames(tmp.dat), ":"))
tmp.ref_base <- substr(tmp.pos_info[,4], start=1, stop=1)
tmp.base <- substr(tmp.pos_info[,4], start=3, stop=3)
tmp.subst_effect <- substr(tmp.pos_info[,5], start=1, stop=1)
data.pos <- data.frame(Contig=tmp.pos_info[,1], Gene=tmp.pos_info[,2], Position=as.numeric(tmp.pos_info[,3]), Substitution=tmp.pos_info[,4], Ref_Base=tmp.ref_base, Base=tmp.base, Effect=tmp.subst_effect, Effect_Full=tmp.pos_info[,5])
rownames(data.pos) <- paste0(data.pos$Contig, ":", data.pos$Position, ":", data.pos$Base)
rm(tmp.dat); invisible(gc())
###################################

###################################
#Get clean positional frequency matrix
mat.pos <- tmp.mat.pos[, curr.samples]
rownames(mat.pos) <- rownames(data.pos)
rm(tmp.mat.pos); invisible(gc())

#Add a minute pseudo-count to ref positions and transform "-1"s to "0s"
mat.pos[mat.pos == 0] <- 10^-12
mat.pos[mat.pos == -1] <- 0

#Make matrix sparse
#=> but only if that won't create overflow trouble for the Matrix package
if (length(mat.pos) < (2^31 / 2)) {mat.pos <- Matrix(mat.pos, sparse=T)}

#"Horizontally" filter positions
#=> to select only those with positive non-reference frequencies
#=> this is a sanity check, as well
tmp.freq.sum <- rowSums(mat.pos > 0)
keep.pos <- tmp.freq.sum > 0
if (length(which(keep.pos)) <= 10) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Skip: too few positions after filtering for non-reference.")); q(save="no")}

#Report
writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", length(which(keep.pos)), "filtered non-reference SNP positions."))

#Filter
mat.pos <- mat.pos[keep.pos, ]
data.pos <- data.pos[keep.pos, ]
###################################

###################################
#Turn mat.pos into per-allele table
mat.tmp <- 1 - as.matrix(mat.pos)
#Account for observations of other alleles only ("1" in the original matrix => everything non-reference)
mat.tmp[mat.tmp == 0] <- 10^-12
#Account for missing observations (zeroes in the original allele freq matrix)
mat.tmp[mat.tmp == 1] <- 0
#Account for positions that were all ref (10^-12 in the original allele freq matrix)
mat.tmp[mat.tmp == 1 - 10^-12] <- 1
#Re-set rownames
rownames(mat.tmp) <- paste(data.pos$Contig, data.pos$Position, data.pos$Ref_Base, sep=":")
#Merge matrices and re-sparsify
if ((length(mat.pos) + length(mat.tmp)) < (2^31 / 2)) {mat.allele <- Matrix(rbind(mat.pos, mat.tmp), sparse = T)} else {mat.allele <- rbind(mat.pos, mat.tmp)}
keep.alleles <- mixedsort(unique(rownames(mat.allele)))     #use gtools::mixedsort to sort naturally
mat.allele <- mat.allele[keep.alleles, ]
rm(mat.tmp); invisible(gc())
#Get per-allele data
tmp.allele_info <- do.call(rbind, strsplit(rownames(mat.allele), ":"))
data.allele <- data.frame(Contig=tmp.allele_info[,1], Position=as.numeric(tmp.allele_info[,2]), Base=tmp.allele_info[,3])
rownames(data.allele) <- rownames(mat.allele)
###################################

###################################
#Filter positions
#=> keep only positions which have been called in at least one oral and at least one stool sample
###################################
#Get summed frequencies per position
# if (length(which(curr.habitats == "saliva")) == 1) {
#   sums.oral.pos <- mat.pos[, curr.habitats == "saliva"]
#   sums.oral.allele <- mat.allele[, curr.habitats == "saliva"]
# } else {
sums.oral.pos <- rowSums(mat.pos[, curr.habitats == "saliva", drop=F], na.rm=T)
sums.oral.allele <- rowSums(mat.allele[, curr.habitats == "saliva", drop=F], na.rm=T)

# if (length(which(curr.habitats == "stool")) == 1) {
#   sums.gut.pos <- mat.pos[, curr.habitats == "stool"]
#   sums.gut.allele <- mat.allele[, curr.habitats == "stool"]
# } else {
sums.gut.pos <- rowSums(mat.pos[, curr.habitats == "feces", drop=F], na.rm=T)
sums.gut.allele <- rowSums(mat.allele[, curr.habitats == "feces", drop=F], na.rm=T)


#Annotate positions and alleles covered in...
#=> at least one oral sample
data.pos$covered.oral <- sums.oral.pos > 0
data.allele$covered.oral <- sums.oral.allele > 0
#=> at least one stool sample
data.pos$covered.gut <- sums.gut.pos > 0
data.allele$covered.gut <- sums.gut.allele > 0

#Annotate positions that are an SNV in...
#=> at least one oral sample
data.pos$snv.oral <- sums.oral.pos > sum(curr.habitats == "saliva") * 10^-12
#=> at least one stool sample
data.pos$snv.gut <- sums.gut.pos > sum(curr.habitats == "feces") * 10^-12
###################################

###################################
#Store data on current taxon
curr.data.taxa <- data.frame(
  SNV_data.available = TRUE,
  SNV_data.oral_gut = TRUE,
  n.pos.total = nrow(data.pos),
  n.allele.total = nrow(data.allele),
  #Oral data
  n.pos.oral = sum(data.pos$covered.oral),
  n.snv.oral = sum(data.pos$snv.oral),
  n.allele.oral = sum(data.allele$covered.oral),
  #Oral <-> gut data
  n.pos.shared.oral_gut = sum(data.pos$covered.oral & data.pos$covered.gut),
  n.snv.shared.oral_gut= sum(data.pos$snv.oral & data.pos$snv.gut),
  n.allele.shared.oral_gut = sum(data.allele$covered.oral & data.allele$covered.gut)
)
###################################

###################################
#Store
save(data.allele, mat.allele, file=paste0(PARAM$folder.data, "SNV.data/", tax.id, ".allele_data.filtered.RData"))
save(data.pos, mat.pos, file=paste0(PARAM$folder.data, "SNV.data/", tax.id, ".SNV_data.filtered.RData"))
save(curr.data.taxa, file=paste0(PARAM$folder.data, "SNV.data/", tax.id, ".taxa_data.RData"))
#Beautify output
writeLines("\n")
################################################################################
################################################################################







q(save="no")



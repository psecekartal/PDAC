#!/usr/bin/Rscript
################################################################################
#Project: Pancreatic Cancer Microbiome
#
#Calculate pairwise community (dis)similarities for 16S data.
#
#2017-09-26
#sebastian.schmidt@embl.de
################################################################################


################################################################################
################################################################################
#Source convenience functions
source("https://raw.githubusercontent.com/defleury/Schmidt_et_al_2016_community_similarity/master/functions.community_similarity.R")

#Set basic parameters
PARAM <- list()
PARAM$folder.R <- getwd()
PARAM$folder.base <- gsub("R", "", PARAM$folder.R)
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$use.cores <- 40
################################################################################
################################################################################


################################################################################
################################################################################
#Load data
load(paste0(PARAM$folder.data, "data.count_table.16S.cleaned.Rdata"))
ot <- ano.t
#Normalize to relative abundances
ot <- t(t(ot) / colSums(ot))
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate sample dissimilarities
#-> note: all metrics are calculated as similarities, not dissimilarities
beta.div <- list();

#Bray-Curtis
beta.div[["Bray_Curtis"]] <- community.similarity.par(ot, distance="bray_curtis", use.cores=4)

#Bray-Curtis on sqrt-transformed data
beta.div[["Bray_Curtis.sqrt"]] <- community.similarity.par(sqrt(ot), distance="bray_curtis", use.cores=4)

#Weighted Jaccard
beta.div[["Jaccard_w"]] <- community.similarity.par(ot, distance="jaccard.abd.frac", use.cores=4)

#Get pairwise (raw) SparCC correlations and derived "interaction matrix"
ot.sparcc <- sparcc(ot, size.thresh=0, pseudocount=10^-6, nblocks=4, use.cores=PARAM$use.cores)
ot.sparcc.S <- 0.5 * (cor.par(ot.sparcc, nblocks=100, use.cores=PARAM$use.cores) + 1)

#Unweighted TINA
beta.div[["TINA_uw"]] <- community.similarity.corr.par(ot, S=ot.sparcc.S, distance="jaccard.corr.uw.norm", blocksize=100, use.cores=PARAM$use.cores)
beta.div[["TINA_uw"]][beta.div[["TINA_uw"]] < 0] <- 0

#Weighted TINA
beta.div[["TINA_w"]] <- community.similarity.corr.par(ot, S=ot.sparcc.S, distance="jaccard.corr.w.norm", blocksize=100, use.cores=PARAM$use.cores)
beta.div[["TINA_w"]][beta.div[["TINA_w"]] < 0] <- 0

#Store in R format
save(ot.sparcc, file=paste0(PARAM$folder.data, "16S.sparcc.Rdata"))
save(ot.sparcc.S, file=paste0(PARAM$folder.data, "16S.sparcc.S.Rdata"))
save(beta.div, file=paste0(PARAM$folder.data, "16S.beta_div.Rdata"))
################################################################################
################################################################################

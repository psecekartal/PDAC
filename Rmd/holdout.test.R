# Application of model on the holdout datasets
holdoutTest <- function(studyDataframe, 
                        filename, 
                        metaTable, 
                        expTable, 
                        status, 
                        models, 
                        modelname){
  print(paste0(filename, status))
  # calculate relative abundance 
  feat.holdout <- prop.table(as.matrix(studyDataframe),2)
  # remove unmapped reads
  feat.holdout <- feat.holdout[!row.names(feat.holdout) %in% 
                                 ("-1"), ]
  # get metadata
  # select only metagenomes
  expTable= expTable[expTable$experiment_type=="metagenome" & 
                       expTable$library_strategy=="WGS",]
  meta.holdout <- unique(merge(metaTable, 
                               expTable, 
                               by="sample_alias", 
                               all.x=TRUE))
  meta.holdout=meta.holdout[!duplicated(meta.holdout$ena_ers_sample_id), ]
  print(dim(meta.holdout))
  # fix all rownames
  rownames(meta.holdout) <- meta.holdout$ena_ers_sample_id
  # subset only specific group
  meta.holdout <- meta.holdout[meta.holdout$subject_disease_status %in% 
                                 status, ]
  #choose only fecal samples
  meta.holdout <- meta.holdout[meta.holdout$environment_material=="feces [ENVO:00002003]",]
  
  # get feat table
  feat.holdout <- feat.holdout[, match(rownames(meta.holdout), 
                                       colnames(feat.holdout))]
  # remove unnecessary row character
  rownames(feat.holdout) <- str_replace_all(rownames(feat.holdout),
                                            c("ef_mOTU_v25" = "",
                                              "incertae sedis" = "",
                                              "eta_mOTU_v25" = ""))
  # remove NA columns
  feat.holdout= feat.holdout[, !apply(is.na(feat.holdout), 2, all)]
  meta.holdout <- meta.holdout[match(colnames(feat.holdout),
                                     rownames(meta.holdout)),]
  # create siamcat object
  siamcat.holdout <- siamcat(feat=feat.holdout, 
                             meta=meta.holdout,
                             verbose=0) 
  # make predictions
  siamcat.holdout <- make.predictions(siamcat = models, 
                                      siamcat.holdout = siamcat.holdout)
  pred <- rowMeans(pred_matrix(siamcat.holdout))
  #hist(pred)
  # find threshold
  threshold <- models@eval_data$roc$thresholds[
    which(models@eval_data$roc$specificities > 0.90)[1]]
  # calculate FPR
  df.fpr <- tibble(pred=pred, Group=(label(siamcat.holdout))[[1]])
  df.fpr <- df.fpr %>% 
    group_by(Group) %>% 
    dplyr::summarise(n=n(), 
                     pred.pos=sum(pred>threshold), 
                     pred.neg=sum(pred < threshold)) %>% 
    ungroup() %>% 
    mutate(fpr=pred.pos/n)
  print(df.fpr)
  # export predictions per study
  pred=as.data.frame(pred)
  rownames(pred) <- colnames(siamcat.holdout@phyloseq@otu_table)
  pred = merge(pred, meta.holdout,
               by.x="row.names",
               by.y="ena_ers_sample_id",
               all.x=TRUE)
  pred$study <- paste0(filename)
  pred$model <- paste0(modelname)
  # save predictions for each study 
  write.table(pred,
              paste0(PARAM$folder.results, 
                     'pred.0.90.spec.', 
                     modelname, ".", 
                     filename, ".", 
                     status, '.tsv'),
              row.names=TRUE, sep="\t")
  return(df.fpr)
} 
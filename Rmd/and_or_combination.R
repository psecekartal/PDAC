PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder <- gsub("Rmd/", "", PARAM$folder.R)
PARAM$folder.files <- paste0(PARAM$folder, "files/")
PARAM$folder.Rdata <- paste0(PARAM$folder, "Rdata/")
PARAM$folder.results <- paste0(PARAM$folder, "results/")
# load siamcat objects
load(paste0(PARAM$folder.Rdata, "2020-05-03.pc.st.cancer.3mg.lasso_ll.log.clr.prev0.siamcat.Rdata"))
siamcat.pc <- siamcat

load(paste0(PARAM$folder.Rdata, "/pos.pc.st.cancer.3mg.lasso_ll.log.clr.prev0.50gFC.Rdata"))
siamcat.pos <- siamcat
#reassign ca19 
a=meta(siamcat.pc)
meta.st.cancer$CA19.new <- as.numeric(gsub("<","",meta.st.cancer$CA19))
meta(siamcat.pc) <- a
meta(siamcat.pos) <- meta(siamcat.pc)
# function for AND and OR combination
.f_and_or_combination <- function(x){
  require("tidyverse")
  require("pROC")
  require(SIAMCAT)
  stopifnot(class(x) == 'siamcat')
  stopifnot('coded_ca19' %in% colnames(meta(x)))
  stopifnot(!is.null(eval_data(x, verbose=0)))
  temp <- eval_data(x)$roc

  # compute for AND combination and OR combination separately
  pred <- enframe(rowMeans(pred_matrix(x)),
                  name='Sample_ID', value="pred_status") %>% 
    full_join(tibble(Sample_ID=rownames(meta(x)),
                     ca_status=meta(x)$coded_ca), by='Sample_ID') %>% 
    full_join(enframe(label(x)$label, 
                      name='Sample_ID', value='Group'), by='Sample_ID')
  pred$ca_status <- as.numeric(pred$ca_status)
  # AND combo
  combo <- pred %>% 
    filter(ca_status!=0) %>% 
    mutate(ca_status=case_when(ca_status==-1~0, TRUE~ca_status)) %>% 
    mutate(combo=pred_status*ca_status)
  roc.and.red <- roc(predictor=combo$combo, response=combo$Group, 
                     levels = c(-1, 1), direction = '<')
  
  combo <- pred %>% 
    mutate(ca_status=case_when(ca_status==-1~0, TRUE~ca_status)) %>% 
    mutate(combo=pred_status*ca_status)
  roc.and <- roc(predictor=combo$combo, response=combo$Group, 
                 levels = c(-1, 1), direction = '<')
  
  # OR combo
  combo <- pred %>% 
    filter(ca_status!=0) %>% 
    mutate(ca_status=case_when(ca_status==-1~0, TRUE~ca_status)) %>% 
    mutate(combo=pred_status+ca_status)
  roc.or.red <- roc(predictor=combo$combo, response=combo$Group, 
                    levels = c(-1, 1), direction = '<')
  
  combo <- pred %>% 
    mutate(ca_status=case_when(ca_status==-1~0, TRUE~ca_status)) %>% 
    mutate(combo=pred_status+ca_status)
  roc.or <- roc(predictor=combo$combo, response=combo$Group, 
                levels = c(-1, 1), direction = '<')
  
  df.res <- list('or'=roc.or, 'or.red'=roc.or.red, 
                 'and'=roc.and, 'and.red'=roc.and.red)
  return(df.res)
}
# combine Ca19-9 with microbiome
sc.ca19=.f_and_or_combination(siamcat.pc)
sc.ca19.pos=.f_and_or_combination(siamcat.pos)
save(sc.ca19.pos,sc.ca19, file = paste0(PARAM$folder.results, 'ca19.combined.Rdata'))
# plot ROCs
get.roc <- function(sc, name){
  roc.obj <- sc$or
  df.return <- tibble(sensitivity=roc.obj$sensitivities,
                      specificitiy=roc.obj$specificities,
                      name=name,
                      auroc=roc.obj$auc)
}

df.plot <- bind_rows(get.roc(sc.ca19, 'OR'),
                     get.roc(sc.ca19.pos, 'OR.pos'))

p <- df.plot %>% 
  arrange(sensitivity) %>% 
  mutate(label=paste0(name, ' ', sprintf(fmt='%.2f', auroc))) %>% 
  ggplot(aes(x=specificitiy, y=sensitivity, col=label)) + 
  geom_abline(intercept = 1, slope = 1, col='black', lty=3) +
  geom_line() + 
  #scale_color_tableau(name='') + 
  xlab('Specificity') + 
  ylab('Sensitivity') + 
  theme_bw() + 
  scale_x_reverse() +
  theme(panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.65, 0.2), legend.title = element_blank())
# save the plots
ggsave(p,  filename = paste0(PARAM$folder.results, 'roc.ca19.and.or.pdf'), width = 5.5, height = 4.5)

# TRP rate for stages
siamcat=siamcat.pc

pred <- rowMeans(pred_matrix(siamcat))
threshold <- siamcat@eval_data$roc$thresholds[
  which(siamcat@eval_data$roc$specificities > 0.90)[1]]

pred=as.data.frame(pred)
rownames(pred) <- colnames(siamcat@phyloseq@otu_table)
pred = merge(pred, meta.st.cancer, by="row.names", all.x=TRUE)
table(pred$stage)

pred <- pred %>%  
  mutate(stage.n=dplyr::recode(stage,  "1"=1,"2"=1,"3"=3,"4"=3))

stage=pred$stage.n
df.tpr <- tibble(pred=pred$pred, Group=(label(siamcat))[[1]], stage=stage)
df.tpr <- df.tpr %>% 
  group_by(stage) %>% 
  dplyr::summarise(n=n(), 
                   pred.pos=sum(pred>threshold), 
                   pred.neg=sum(pred < threshold)) %>% 
  ungroup() %>% 
  mutate(tpr=pred.pos/n)

df.tpr.m1=df.tpr
df.tpr.m1$model="Model 1"

df.tpr.m2 =df.tpr
df.tpr.m2$model="Model 2"

# repeat for combinations
siamcat=sc.ca19.pos

pred <- siamcat$or$predictor
threshold <- siamcat$or$thresholds[which(siamcat$or$specificities > 0.90)[1]]

pred=as.data.frame(pred)
rownames(pred) <- colnames(siamcat.pc@phyloseq@otu_table)
pred = merge(pred, meta.st.cancer, by="row.names", all.x=TRUE)
table(pred$stage)

pred <- pred %>%  
  mutate(stage.n=dplyr::recode(stage,  "1"=1,"2"=1,"3"=3,"4"=3))

stage=pred$stage.n

df.tpr <- tibble(pred=pred$pred, Group=(label(siamcat.pc))[[1]], stage=stage)
df.tpr <- df.tpr %>% 
  group_by(stage) %>% 
  dplyr::summarise(n=n(), 
                   pred.pos=sum(pred>threshold), 
                   pred.neg=sum(pred < threshold)) %>% 
  ungroup() %>% 
  mutate(tpr=pred.pos/n)

 df.tpr.ca19.m1=df.tpr
 df.tpr.ca19.m1$model="Model 1 OR CA19-9"

df.tpr.c19.m2=df.tpr
df.tpr.c19.m2$model = "Model 2 OR CA19-9"

#combine all tpr rate for plotting
df= rbind(df.tpr.m1, df.tpr.m2, df.tpr.ca19.m1, df.tpr.c19.m2)
# add CP cases
df[13:14,] <- ""
df$stage[13:14] <- c("CP", "CP")
df$tpr[13:14] <- c(0.14, 0.17)
df$model[13:14] <- c("Model 1","Model 2")
df$tpr<-as.numeric(df$tpr)
# plot tpr
p.tpr <- ggplot(data=df, aes(x=stage, y=tpr, fill=stage)) + 
  geom_col(alpha=0.8) + 
  theme_classic()  +  facet_grid(cols = vars(model), scales = "free") 


ggsave(p.tpr, filename=paste0(PARAM$folder.results, "stage.tpr.barplot.cpincluded.pdf"), 
       width = 5, height=6)


# get total TPR
siamcat=sc.ca19.pos
pred <- siamcat$or$predictor
threshold <- siamcat$or$thresholds[which(siamcat$or$specificities > 0.90)[1]]
pred=as.data.frame(pred)
df.fpr <- tibble(pred=pred$pred, Group=(label(siamcat.pc))[[1]])

df.fpr <- df.fpr %>%
  group_by(Group) %>%
  dplyr::summarise(n=n(),
                   pred.pos=sum(pred > threshold),
                   pred.neg=sum(pred < threshold)) %>%
  ungroup() %>%
  mutate(fpr=pred.pos/n)



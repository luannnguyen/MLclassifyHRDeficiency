#========= Load packages =========#
packages <- c(#'glmnet',
              'ggplot2',
              'ggfortify',
              'dplyr',
              'reshape2',
              'randomForest',
              'ROCR')

ipak <- function(pkg){
   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
   if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
   sapply(pkg, require, character.only = TRUE)
}

ipak(packages)

options(stringsAsFactors = F)

#========= Misc. functions =========#
#--------- Normalize data according to HRDetect paper ---------#
normalize_data <- function(df, na.replace = T){
   df_norm <- apply(df,2,function(v){
      #v <- as.numeric(v)
      ln_v <- log(v+1, exp(1))
      norm_v <- (ln_v - mean(ln_v, na.rm = T) ) / sd(ln_v, na.rm = T)
      
      if(na.replace == T){
         norm_v[is.na(norm_v)] <- median(norm_v, na.rm = T) ## replace NA with median
      }
      
      return(norm_v)
   })
   return( as.data.frame(df_norm) )
}

#--------- Calculate false positive and negative rate ---------#
get_fPosNegRate <- function(v_logical_expected,v_logical_predicted)
{
   # v_logical_expected <- rf_vs_annotation$Gene_simple
   # v_logical_predicted <- rf_vs_annotation$RFpBRCA_umcuData_sangerTrained > 0.62
   df <- cbind(v_logical_expected,v_logical_predicted) %>% as.data.frame()
   
   fpos <- nrow(df[df$v_logical_expected == 1 & df$v_logical_predicted == 0,])
   fneg <- nrow(df[df$v_logical_expected == 0 & df$v_logical_predicted == 1,])
   
   tpos <- nrow(df[df$v_logical_expected == 1 & df$v_logical_predicted == 1,])
   tneg <- nrow(df[df$v_logical_expected == 0 & df$v_logical_predicted == 0,])
   
   fpos_rate <- fpos/(fpos+tneg)
   fneg_rate <- fneg/(fneg+tpos)
   
   tpos_rate <- tpos/(tpos+fpos)
   tneg_rate <- tneg/(tneg+fneg)
   
   fPosNegRate <- c(fpos_rate,fneg_rate,tpos_rate,tneg_rate)
   names(fPosNegRate) <- c('fpos_rate','fneg_rate','tpos_rate','tneg_rate')
   
   return(fPosNegRate)
}

#========= Random forest functions =========#
#--------- Select mtry where out of bag (OOB) error is minimal ---------#
## Each decision tree node is split using a number of features (mtry)
tuneRF_iterate <- function(df, colname_response, tuneRF_iterations = 20)
{
   ## run random forest with multiple mtry's
   # df <- normdata_ss
   # colname_response <- 'brca_deficiency'
   v_test_mtry <- sapply(1:tuneRF_iterations, function(i)
   {
      mtry <- tuneRF(x = df[,-which(colnames(df) == colname_response)], #df of features/observations
                     y = df[,colname_response], ## vector of expected response
                     ntreeTry=500, #number of trees to create
                     stepFactor=1.3, #inflation rate of mtry for each iteration
                     improve=0.01, #relative improvement in OOB error must be by this much for the search to continue
                     trace=F,
                     plot=F
      )
      
      mtrys_minOOBerror <- mtry[mtry[,2] == min(mtry[,2]),1]
      mtrys_minOOBerror[length(mtrys_minOOBerror)] ## always pick the higher mtry value if two rows have the same minOOBerror
      #mtry[(mtry[,2] == min(mtry[,2])) %>% which() %>% min(),][1]
   })
   
   ## get rounded average
   mtry_best <- mean(v_test_mtry) %>% round(0)
   return(mtry_best)
}

#--------- Get bagged training set ---------#
# df <- normdata_ss
# colname_response <- 'brca_deficiency'

get_baggedTrainingSet <- function(df){
   train_set <- sample(1:nrow(df), (nrow(df) * 0.63212056) %>% round(0) ) ## 0.632 rule; choose ~2/3 of data set for bagging
   df_ss <- df[train_set,]
   df_ss <- df_ss[rownames(df_ss) %>% order,]
   
   return(df_ss)
}

#--------- Up sample low occurrence observations ---------#
df_ss <- get_baggedTrainingSet(df)

balanceClasses <- function(df)

#--------- Create and combine multiple random forests ---------#
#set.seed(Sys.time())
randomForest_iterate <- function(df, colname_response, mtry_best, randomForest_iterations = 100)
{
   # df <- normdata_ss
   # colname_response <- 'brca_deficiency'
   # mtry_best <- mtry_best_sanger
   
   ## run multiple random forests
   l_RF <- lapply(1:randomForest_iterations,function(i)
   {
      # train_set <- sample(1:nrow(df), (nrow(df) * 0.63212056) %>% round(0) ) ## 0.632 rule; choose ~2/3 of data set for bagging
      # df_ss <- df[train_set,]
      
      df_ss <- get_baggedTrainingSet(df)
      
      randomForest(formula = as.formula(paste0('df_ss$',colname_response,' ~ .')),
                   data = df_ss,
                   ntree = 500,
                   mtry = mtry_best,
                   importance = T)
   })
   
   ## combine random forests into a consensus forest
   RF_combined <-
      paste0('l_RF[[',1:randomForest_iterations,']]', collapse = ',') %>% ## paste 'l_RF[[1]],l_RF[[2]],...'
      paste0('combine(',.,')') %>% ## paste 'combine(l_RF[[1]],l_RF[[2]],...)'
      parse(text = .) %>% eval() ## evaluate string as expression
   
   return(RF_combined)
}

#========= Data prep =========#
#--------- Prepare umcu data ---------#
umcu_signature_matrices_path <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/matrices/umcu_matrices/'

umcu_snv_tsv <- paste0(umcu_signature_matrices_path,'final_table_SNV.txt')
umcu_indel_tsv <- paste0(umcu_signature_matrices_path,'final_table_indels.txt')
umcu_sv_tsv <- paste0(umcu_signature_matrices_path,'final_table_SV.txt')

read.tsv <- function(file){ 
   read.table(file,sep='\t',header=T)
}

df_umcu_snv <- read.tsv(umcu_snv_tsv)
df_umcu_indel <- t(read.tsv(umcu_indel_tsv)) %>% as.data.frame()
df_umcu_sv <- read.tsv(umcu_sv_tsv)

## rm 'a' from donor id
rownames(df_umcu_snv) <- rownames(df_umcu_snv) %>% gsub('a','',.)
rownames(df_umcu_indel) <- rownames(df_umcu_indel) %>% gsub('a','',.)

## intersect donors
l_rownames <- list(rownames(df_umcu_snv),rownames(df_umcu_indel),rownames(df_umcu_sv))
common_donors <- Reduce(intersect, l_rownames)

## only keep common donors
df_umcu_snv_ss <- df_umcu_snv[rownames(df_umcu_snv) %in% common_donors,]
df_umcu_indel_ss <- df_umcu_indel[rownames(df_umcu_indel) %in% common_donors,]
df_umcu_sv_ss <- df_umcu_sv[rownames(df_umcu_sv) %in% common_donors,]

## replace snv signature colnames with sanger format
colnames(df_umcu_snv_ss) <- colnames(df_umcu_snv_ss) %>% gsub('Signature','e',.)

## replace sv signature colnames with sanger format
colnames(df_umcu_sv_ss) <- colnames(df_umcu_sv_ss) %>% gsub('SV_Signature_','SV',.)

## combine repeat insertion/deletion 1,2,3 contexts; calculate relative contribution for all contexts
# df_umcu_indel_ss %>% head()
# colnames(rawdata)

df_umcu_indel_ss_rowsums <- df_umcu_indel_ss %>% apply(.,1,sum)

del.rep.prop <- (df_umcu_indel_ss[,c(1,3,5)] %>% apply(.,1,sum)) / df_umcu_indel_ss_rowsums
del.mh.prop <- df_umcu_indel_ss$micro_hom_deletion / df_umcu_indel_ss_rowsums

ins.rep.prop <- (df_umcu_indel_ss[,c(2,4,6)] %>% apply(.,1,sum)) / df_umcu_indel_ss_rowsums
ins.mh.prop <- df_umcu_indel_ss$micro_hom_insertion / df_umcu_indel_ss_rowsums

ins <- ( ( df_umcu_indel_ss[,c(2,4,6)] %>% apply(.,1,sum) ) + df_umcu_indel_ss$micro_hom_insertion ) / df_umcu_indel_ss_rowsums

indel.none.prop <- df_umcu_indel_ss$no_context / df_umcu_indel_ss_rowsums

df_umcu_indel_ss2 <- cbind(df_umcu_indel_ss[,1], del.rep.prop, del.mh.prop,ins.rep.prop,ins.mh.prop,indel.none.prop)
df_umcu_indel_ss2 <- df_umcu_indel_ss2[,-1]

df_umcu_indel_ss2_sangerFormat <- cbind(df_umcu_indel_ss[,1], del.rep.prop, del.mh.prop, ins, indel.none.prop)
df_umcu_indel_ss2_sangerFormat <- df_umcu_indel_ss2_sangerFormat[,-1]

## final normalized signature martrix
umcu_normdata <- cbind(df_umcu_snv_ss,df_umcu_indel_ss2,df_umcu_sv_ss) %>% normalize_data()
umcu_normdata_indelSangerFormat <- cbind(df_umcu_snv_ss,df_umcu_indel_ss2_sangerFormat,df_umcu_sv_ss) %>% normalize_data()

#--------- Sanger data ---------#
normdata_ss <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/matrices/signature_matrices/sanger_pipeline_signature_matrix_norm.tsv',
                       sep = '\t', header = T)

rownames(normdata_ss) <- rownames(normdata_ss) %>% 
   gsub('a2','',.) %>% 
   gsub('a','',.) %>%
   gsub('b','',.)

rawdata <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/nm.4292-S2_raw_data.txt',
                      sep = '\t', header = T)

rawdata$Sample <- rawdata$Sample %>% 
   gsub('a2','',.) %>% 
   gsub('a','',.) %>%
   gsub('b','',.)

#--------- Assign brca deficiency columns ---------#
v_brca_deficiency <- as.integer(nchar(as.vector(rawdata$Gene)) > 1) %>% as.factor()
names(v_brca_deficiency) <- rawdata$Sample

normdata_ss$brca_deficiency <- v_brca_deficiency

common_donors <- intersect(rownames(normdata_ss),rownames(umcu_normdata))

umcu_normdata$brca_deficiency <- v_brca_deficiency[names(v_brca_deficiency) %in% common_donors]
umcu_normdata_indelSangerFormat$brca_deficiency <- v_brca_deficiency[names(v_brca_deficiency) %in% common_donors]

#========= Training random forest on sanger data =========#
## Test run random forest
# ## Determine in-bag samples and run random forest
# rf_train <- sample(1:nrow(normdata_ss), (nrow(normdata_ss) * 0.63212056) %>% round(0) ) ## 0.632 rule; choose ~2/3 of data set for training
# rf_normdata_ss <- randomForest(formula = brca_deficiency ~ ., data = normdata_ss, subset = rf_train, ntree = 500)

## Get best mtry
mtry_best_sanger <- tuneRF_iterate(normdata_ss, 'brca_deficiency')

## training and prediction with all available features in sanger data set
# rf2_normdata_ss <- randomForest_iterate(normdata_ss, 'brca_deficiency',randomForest_iterations=10)
# rf_pBRCA <- predict(object = rf2_normdata_ss, newdata = normdata_ss[-23], type = "prob")

## training and prediction without features unavailable in umcu pipeline
rf2_normdata_ss <- randomForest_iterate(normdata_ss[,-c(5,22)],'brca_deficiency',mtry_best_sanger,randomForest_iterations=100)

RFpBRCA_sangerData_sangerDataTrained <- predict(object = rf2_normdata_ss, newdata = normdata_ss[,-23], type = "prob")
RFpBRCA_umcuData_sangerDataTrained <- predict(object = rf2_normdata_ss, newdata = umcu_normdata_indelSangerFormat[,colnames(normdata_ss[,-c(5,22,23)])], type = "prob")

#========= Train with umcu data =========#
## Get best mtry's
mtry_best_umcu <- tuneRF_iterate(umcu_normdata, 'brca_deficiency')

## Train RF
rf_umcu_normdata <- randomForest_iterate(umcu_normdata,'brca_deficiency',mtry_best_umcu,randomForest_iterations=100)

## Predict
RFpBRCA_umcuData_umcuDataTrained <- predict(object = rf_umcu_normdata,newdata = umcu_normdata[,-42],type = "prob")

#========= Train with 50:50 sanger:umcu data =========#
## Data prep
#umcu_normdata
normdata_ss2 <- normdata_ss[rownames(normdata_ss) %in% common_donors,]

s1 <- sample(1:length(common_donors),round(length(common_donors)/2))
s2 <- which(!(1:length(common_donors) %in% s1))

common_cols <- intersect(colnames(normdata_ss2),colnames(umcu_normdata))

mix_normdata <- rbind(umcu_normdata[s1,common_cols], normdata_ss2[s2,common_cols])
mix_normdata <- mix_normdata[order(rownames(mix_normdata)),]

## Get best mtry's
mtry_best_mix <- tuneRF_iterate(mix_normdata, 'brca_deficiency')

## Train RF
rf_mix_normdata <- randomForest_iterate(mix_normdata,'brca_deficiency',mtry_best_mix,randomForest_iterations=100)

## Predict
RFpBRCA_mixData_mixDataTrained <- predict(object = rf_mix_normdata,newdata = mix_normdata,type = "prob")

RFpBRCA_sangerData_mixDataTrained <- predict(object = rf_mix_normdata,newdata = normdata_ss[,common_cols[-20]],type = "prob")
RFpBRCA_umcuData_mixDataTrained <- predict(object = rf_mix_normdata,newdata = umcu_normdata[,common_cols[-20]],type = "prob")

#========= Train with umcu data; use expanded indel contexts =========#
## data prep
# Load expanded indel matrices
expandedMatricesFiles <- list.files('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_ICGC/Breast_560/Indels/',pattern = '*.txt',full.names = T)
df_umcu_indelExpanded <- lapply(expandedMatricesFiles,function(i){ 
   df <- read.table(i,header = T,sep='\t')
   rownames(df) <- df[,1]
   df[,1] <- NULL
   return(df)
})
df_umcu_indelExpanded <- do.call(cbind,df_umcu_indelExpanded) %>% t() %>% as.data.frame()
rownames(df_umcu_indelExpanded) <- rownames(df_umcu_indelExpanded) %>% 
   gsub('a2','',.) %>% 
   gsub('a','',.) %>%
   gsub('b','',.)

df_umcu_indelExpanded_ss <- df_umcu_indelExpanded[rownames(df_umcu_indelExpanded) %in% common_donors,]
colnames(df_umcu_indelExpanded_ss) <- paste0(c('del.rep.','ins.rep.','del.mh.','ins.mh.','indel.none.') %>% rep(.,each=5),1:5,'.prop')

umcu_normdata_expandedIndel <- cbind(df_umcu_snv_ss,df_umcu_indelExpanded_ss/apply(df_umcu_indelExpanded_ss,1,sum),df_umcu_sv_ss) %>% normalize_data()
umcu_normdata_expandedIndel$brca_deficiency <- v_brca_deficiency[names(v_brca_deficiency) %in% common_donors]

## Get best mtry's
mtry_best_umcu_normdata_expandedIndel <- tuneRF_iterate(umcu_normdata_expandedIndel, 'brca_deficiency')

## Train RF
rf_umcu_normdata_expandedIndel <- randomForest_iterate(umcu_normdata_expandedIndel,
                                                       'brca_deficiency',
                                                       mtry_best_umcu_normdata_expandedIndel,
                                                       randomForest_iterations=100)

## Predict
RFpBRCA_umcuDataExpandedIndel_umcuDataExpandedIndelTrained <- predict(object = rf_umcu_normdata_expandedIndel,
                                                                      newdata = umcu_normdata_expandedIndel[,-62],
                                                                      type = "prob")

## Combine del.mh.prop into one feature
# umcu_normdata_expandedIndel_topFeat <- importance(rf_umcu_normdata_expandedIndel) %>% .[. > median(.),] %>% names()
# umcu_normdata_expandedIndelFeatSel <- umcu_normdata_expandedIndel[,colnames(umcu_normdata_expandedIndel) %in% umcu_normdata_expandedIndel_topFeat]
# umcu_normdata_expandedIndelFeatSel$brca_deficiency <- umcu_normdata_expandedIndel$brca_deficiency

mtry_best_umcu_normdata_expandedIndelFeatSel <- tuneRF_iterate(umcu_normdata_expandedIndel, 'brca_deficiency')

umcu_normdata_expandedIndelFeatSel <- cbind(umcu_normdata_expandedIndel[,-(41:45)],df_umcu_indel_ss2 %>% as.data.frame() %>% .$del.mh.prop)
colnames(umcu_normdata_expandedIndelFeatSel)[58] <- 'del.mh.prop'

rf_umcu_normdata_expandedIndelFeatSel <- randomForest_iterate(umcu_normdata_expandedIndelFeatSel,
                                                       'brca_deficiency',
                                                       mtry_best_umcu_normdata_expandedIndelFeatSel,
                                                       randomForest_iterations=100)

#========= Create rf_vs_annotation for comparison between sanger and umcu data/training =========#
rf_vs_annotation <- rawdata[rawdata$Sample %in% common_donors,c("Sample","Gene")]
rf_vs_annotation$Gene_simple <- as.integer(nchar(as.vector(rf_vs_annotation$Gene)) > 0)

rf_vs_annotation$RF1_sangerData_sangerDataTrained <- RFpBRCA_sangerData_sangerDataTrained[rownames(RFpBRCA_sangerData_sangerDataTrained) %in% common_donors,2]

rf_vs_annotation$RF2_umcuData_sangerDataTrained <- RFpBRCA_umcuData_sangerDataTrained[rownames(RFpBRCA_umcuData_sangerDataTrained) %in% common_donors,2]
rf_vs_annotation$RF3_umcuData_umcuDataTrained <- RFpBRCA_umcuData_umcuDataTrained[,2]

rf_vs_annotation$RF4_mixData_mixDataTrained <- RFpBRCA_mixData_mixDataTrained[,2]
rf_vs_annotation$RF5_sangerData_mixDataTrained <- RFpBRCA_sangerData_mixDataTrained[rownames(RFpBRCA_sangerData_mixDataTrained) %in% common_donors,2]
rf_vs_annotation$RF6_umcuData_mixDataTrained <- RFpBRCA_umcuData_mixDataTrained[,2]

rf_vs_annotation$RF7_umcuDataExpandedIndel_umcuDataExpandedIndelTrained <- RFpBRCA_umcuDataExpandedIndel_umcuDataExpandedIndelTrained[,2]

#========= Performance metrics =========#
#--------- Variable importance; MeanDecreaseAccuracy ---------#
pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/random_forest_training/MeanDecreaseGini_rfTraining_sangerData_umcuData.pdf', 10, 6)
varImpPlot(rf2_normdata_ss,type=1, main = 'rfTraining_sangerData')
varImpPlot(rf_umcu_normdata,type=1, main = 'rfTraining_umcuData')
varImpPlot(rf_mix_normdata,type=1, main = 'rfTraining_mixData')
# varImpPlot(rf_umcu_normdata_expandedIndel,type=2, main = 'rfTraining_umcuData_expandedIndel')
# varImpPlot(rf_umcu_normdata_expandedIndelFeatSel,type=2, main = 'rfTraining_umcuData_expandedIndelFeatSel')
dev.off()

# #--------- ROC curve and AUC ---------#
# sapply(4:ncol(rf_vs_annotation),function(i){
#    pred <- prediction(rf_vs_annotation[,i], rf_vs_annotation$Gene_simple)
#    perf <- performance(pred,'auc')
#    perf@y.values[[1]]
# })

# pred <- prediction(rf_vs_annotation$RF3_umcuData_umcuDataTrained, rf_vs_annotation$Gene_simple)
# perf <- performance(pred,"tpr","fpr")
# performance(pred, "auc")
# plot(perf,colorize=TRUE)

#--------- Plot true positive and true negative rate across p values ---------#
p_range <- seq(from=0, to=1.0, by=0.01)

## for random forest
l_fpnRates <- lapply(4:ncol(rf_vs_annotation),function(i){
   fpnRates <- lapply(p_range, function(j){
      get_fPosNegRate(rf_vs_annotation$Gene_simple,rf_vs_annotation[,i] > j)
   }) %>% do.call(rbind,.) %>% as.data.frame()
   
   fpnRates$p <- p_range
   fpnRates$analysis <- colnames(rf_vs_annotation)[i]
   
   return(fpnRates[3:6])
   #return(fpnRates)
})

abs(l_fpnRates[[1]][,1]-l_fpnRates[[1]][,2]) %>% min(.)

fpnRates_melt <- do.call(rbind,l_fpnRates)

## for HRDetect logistic regression
brca_annotation_vs_pbrca_ss <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/plots_compare_hrdetect_pipelines/brca_annotation_vs_pbrca_ss_multi_datasets.tsv',
           sep = '\t',header=T)[1:5]

fpnRate_predictorProb <- lapply(p_range, function(j){
   get_fPosNegRate(brca_annotation_vs_pbrca_ss$brca_annotation,brca_annotation_vs_pbrca_ss[,'predictorProb'] > j)
}) %>% do.call(rbind,.) %>% as.data.frame() %>% .[,3:4]
fpnRate_predictorProb$p <- p_range
fpnRate_predictorProb$analysis <- 'LRpBRCA_HRDetect'

## plot true positive and true negative rates
# fpnRates_melt <- melt(fpnRates[,c('p','tpos_rate')], id.vars = 'p')
fpnRates_melt <- rbind(fpnRates_melt,fpnRate_predictorProb)

## separate tpos and tneg rate plots
# tposRates_plot <- ggplot(data=fpnRates_melt,aes(x=p, y=tpos_rate, colour=analysis)) +
#    geom_line() +
#    ggtitle('True positive rates at varying probability cutoffs') +
#    xlab('Probability') +
#    ylab('True positive rate')
# 
# tnegRates_plot <- ggplot(data=fpnRates_melt,aes(x=p, y=tneg_rate, colour=analysis)) +
#    geom_line(linetype = 2) +
#    ggtitle('True negative rates at varying probability cutoffs') +
#    xlab('Probability') +
#    ylab('True negative rate')
# 
# pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/random_forest_training/truePosNegRates_bagged.pdf', 8, 4)
# plot(tposRates_plot)
# plot(tnegRates_plot)
# dev.off()

## combine tpos and tneg plots
tPosPegRates_include <- fpnRates_melt$analysis %>% factor %>% levels %>% .[c(3,4,5,6,7)]
tPosPegRates_plot <- ggplot(data=fpnRates_melt[fpnRates_melt$analysis %in% tPosPegRates_include,]) +
   geom_line(aes(x=p, y=tpos_rate, colour=analysis)) +
   geom_line(aes(x=p, y=tneg_rate, colour=analysis),linetype = 2) +
   ggtitle('Solid: true positive rate; Dashed: true negative rate') +
   xlab('Probability') +
   ylab('True positive/negative rate')

pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/random_forest_training/truePosNegRates_mixData.pdf', 10, 4)
plot(tPosPegRates_plot)
dev.off()


#::::::::: Distinguishing BRCA1 and BRCA1 deficiency :::::::::#
#========= Train with umcu data =========#
#--------- Data prep ---------#
v_brca1_deficiency <- as.integer(as.vector(rawdata$Gene) == 'BRCA1') %>% as.factor()
v_brca2_deficiency <- as.integer(as.vector(rawdata$Gene) == 'BRCA2') %>% as.factor()
# 
# names(v_brca1_deficiency) <- rawdata$Sample
# names(v_brca2_deficiency) <- rawdata$Sample
# 
# umcu_normdata_brca1.2 <- umcu_normdata
# 
# umcu_normdata_brca1.2$brca_deficiency <- NULL
# umcu_normdata_brca1.2$brca1_deficiency <- v_brca1_deficiency[names(v_brca1_deficiency) %in% common_donors]
# umcu_normdata_brca1.2$brca2_deficiency <- v_brca2_deficiency[names(v_brca2_deficiency) %in% common_donors]

v_brca1.2_deficiency <- rawdata[,c('Sample','Gene')]
v_brca1.2_deficiency[nchar(v_brca1.2_deficiency$Gene) == 0,'Gene'] <- 0

v_brca1.2_deficiency$Gene <- v_brca1.2_deficiency$Gene %>% gsub('BRCA','',.)
rownames(v_brca1.2_deficiency) <- v_brca1.2_deficiency$Sample
v_brca1.2_deficiency$Sample <- NULL

umcu_normdata_brca1.2 <- umcu_normdata
umcu_normdata_brca1.2$brca_deficiency <- v_brca1.2_deficiency[rownames(v_brca1.2_deficiency) %in% common_donors,] %>% as.factor

#--------- Test run random forest ---------#
# ## Determine in-bag samples and run random forest
# rf_train <- sample(1:nrow(umcu_normdata_brca1.2), (nrow(umcu_normdata_brca1.2) * 0.63212056) %>% round(0) ) ## 0.632 rule; choose ~2/3 of data set for training
# rf_umcu_normdata_brca1.2 <- randomForest(formula = brca_deficiency ~ ., data = umcu_normdata_brca1.2, subset = rf_train, ntree = 500,importance =T)
# predict(object = rf_umcu_normdata_brca1.2,newdata = umcu_normdata_brca1.2,type = "prob") %>%
#    cbind(., umcu_normdata_brca1.2$brca_deficiency %>% as.integer() - 1)



#--------- Get best mtry's ---------#
mtry_best_umcu_brca1.2 <- tuneRF_iterate(umcu_normdata_brca1.2, 'brca_deficiency')

#--------- Train RF ---------#
rf_umcu_normdata_brca1.2 <- randomForest_iterate(umcu_normdata_brca1.2,'brca_deficiency',mtry_best_umcu_brca1.2,randomForest_iterations=100)

RFpBRCA_umcuDataB1B2_umcuDataB1B2Trained <- predict(object = rf_umcu_normdata_brca1.2,newdata = umcu_normdata_brca1.2,type = "prob")

pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/random_forest_training/MeanDecreaseAcc_rfTraining_BRCA1and2.pdf', 10, 6)
varImpPlot(rf_umcu_normdata_brca1.2,type=1,class = 1,main = 'BRCA1 deficiency')
varImpPlot(rf_umcu_normdata_brca1.2,type=1,class = 2,main = 'BRCA2 deficiency')
dev.off()

# importance(rf_umcu_normdata_brca1.2)
# varImpPlot(rf_umcu_normdata_brca1.2,type=1,class = 1)
# varImpPlot(rf_umcu_normdata_brca1.2,type=1,class = 2)

#--------- RF vs annotation ---------#
rf_vs_annotation2 <- rawdata[rawdata$Sample %in% common_donors,c("Sample","Gene")]
rf_vs_annotation2$B1_exp <- as.integer(rf_vs_annotation2$Gene == 'BRCA1')
rf_vs_annotation2$B2_exp <- as.integer(rf_vs_annotation2$Gene == 'BRCA2')

rf_vs_annotation2 <- cbind(rf_vs_annotation2,RFpBRCA_umcuDataB1B2_umcuDataB1B2Trained[,2:3])
colnames(rf_vs_annotation2)[5:6] <- c('B1_pred', 'B2_pred')

# ## check if there are any predictions of both brca1/2
# rf_vs_annotation2_ss <- rf_vs_annotation2[nchar(rf_vs_annotation2$Gene) > 0,]
# 
# apply(rf_vs_annotation2_ss[3:6],1,function(i){
#    if(i[1] == 1 & i[3] > i[4]){
#       return(1)
#    } else if(i[2] == 1 & i[4] > i[3]){
#       return(1)
#    } else {
#       return(0)
#    }
#    #return(i)
# })

if(rf_vs_annotation2_ss$B1_exp == 1 & rf_vs_annotation2_ss$B1_pred > rf_vs_annotation2_ss$B2_pred){1} else {0}

#--------- True pos/neg rates
p_range <- seq(from=0, to=1.0, by=0.01)

B1_tPosNegRates <- lapply(p_range,function(j){
   get_fPosNegRate(rf_vs_annotation2$B1_exp,rf_vs_annotation2$B1_pred > j)
}) %>% do.call(rbind,.) %>% as.data.frame() %>% .[,c(3,4)]
B1_tPosNegRates$Gene <- 'BRCA1'
B1_tPosNegRates$p <- p_range

B2_tPosNegRates <- lapply(p_range,function(j){
   get_fPosNegRate(rf_vs_annotation2$B2_exp,rf_vs_annotation2$B2_pred > j)
}) %>% do.call(rbind,.) %>% as.data.frame() %>% .[,c(3,4)]
B2_tPosNegRates$Gene <- 'BRCA2'
B2_tPosNegRates$p <- p_range

B1_B2_tPosNegRates <- rbind(B1_tPosNegRates,B2_tPosNegRates)

B1_B2_tPosNegRates_plot <- ggplot(data=B1_B2_tPosNegRates) +
   geom_line(aes(x=p, y=tpos_rate, colour=Gene)) +
   geom_line(aes(x=p, y=tneg_rate, colour=Gene),linetype = 2) +
   ggtitle('Solid: true positive rate; Dashed: true negative rate') +
   xlab('Probability') +
   ylab('True positive/negative rate')

pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/random_forest_training/truePosNegRates_BRCA1and2.pdf', 8, 4)
plot(B1_B2_tPosNegRates_plot)
dev.off()





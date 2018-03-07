#========= Load packages =========#
library(ggplot2)
library(dplyr)
library(reshape2)
library(randomForest)
library(ROCR)
library(stringr)
library(gridExtra)

#========= Global options =========#
options(stringsAsFactors = F)

###############################
#          Functions          #
###############################

#========= Misc. functions =========#
#--------- Convert multiclass response to binary response ---------#
get_simpleResponse <- function(chrVector, Tclasses, Treturn, Fclasses, Freturn, asFactor = TRUE)
{
   
   chrVector <- as.character(chrVector)
   
   chrVector[chrVector %in% Tclasses] <- Treturn
   chrVector[chrVector %in% Fclasses] <- Freturn
   
   if(asFactor == TRUE){
      return(as.factor(chrVector))
   } else {
      return(chrVector)
   }
}

#--------- Integration by trapezoidal rule from data frame giving discrete x and y points ---------#
integrateDiscrete <- function(df, x='x', y='y')
{
   sapply(1:(nrow(df)-1),function(i){
      trapez <- df[i:(i+1),]
      area <- (trapez[2,x] - trapez[1,x])*((trapez[2,y] + trapez[1,y])/2)
      return(area)
   }) %>% sum()
}

#========= Random forest functions =========#
#--------- Get train/test sets ---------#
get_kfoldTrainTest <- function(df, k=10)
{
   ## Total samples
   totalSamples <- nrow(df)
   
   ## Number of desired TEST samples per fold
   nSamples <- totalSamples %/% k
   nRemainder <- totalSamples %% k
   
   ## Create vector of nSamples per fold (to take into account if totalSamples is not a factor of k)
   v_nSamples <- rep(nSamples,k)
   v_nSamples[1:nRemainder] <- v_nSamples[1:nRemainder] + 1
   
   ## Create list with first start/end row number for each fold
   l_StartEnd <- list( c(1,v_nSamples[1]) ) ## initialize start/end row numbers list with first fold
   for(i in 2:length(v_nSamples))
   {
      prevStartEnd <- l_StartEnd[[(i-1)]]
      l_StartEnd[[i]] <- c(prevStartEnd[2]+1, prevStartEnd[2]+v_nSamples[i])
   }
   
   ## Get train and test sets
   l_TrainTest <- lapply(l_StartEnd, function(i)
   {
      Start <- i[1]
      End <- i[2]
      
      testSet <- df[Start:End,]
      trainSet <- df[-(Start:End),]
      
      list(train = trainSet, test = testSet)
   })
   
   return(l_TrainTest)
}

## Depends on: get_kfoldTrainTest()
get_cvTrainTestSets <- function(df, k=10, stratifyByCol = NULL)
{
   ## Shuffle data
   df_shuffled <- df[sample(1:nrow(df)),]
   
   ## No stratification
   if( is.null(stratifyByCol) ){
      return( get_kfoldTrainTest(df_shuffled) )
   }
   
   ## Stratification
   else {
      ## Split data by stratifyByCol
      responseLvls <- df[,stratifyByCol] %>% as.factor() %>% levels()
      l_df_strat <- lapply(responseLvls,function(lvl){ 
         df_shuffled %>% .[.[,stratifyByCol] == lvl ,]
      })
      
      ## Get train/test set per response group
      ## List structure: l_TrainTestStrat[[responseGroup]][[kfold]][[train/test]]
      l_TrainTestStrat <- lapply(l_df_strat, function(df_responseGroup){
         get_kfoldTrainTest(df_responseGroup, k=k)
      })
      
      ## Merge response groups per k-fold
      nResponses <- 1:length(responseLvls)
      
      l_TrainTest <- lapply(1:k,function(k){
         ## Avoided loop for train/set list levels for readability
         trainSet <-
            paste0('l_TrainTestStrat[[',nResponses,']][[k]][[\'train\']]', collapse=',') %>% 
            paste0('rbind(',.,')') %>%
            parse(text = .) %>%
            eval()
         testSet <-
            paste0('l_TrainTestStrat[[',nResponses,']][[k]][[\'test\']]', collapse=',') %>% 
            paste0('rbind(',.,')') %>%
            parse(text = .) %>%
            eval()
         list(train = trainSet, test = testSet)
      })
      
      return(l_TrainTest)
   }
}

#--------- Balance classes ---------#
balanceClasses <- function(df, colname_response, scaling = 'up')
{
   # df <- df_umcuNormData
   
   classFreq <- df[,colname_response] %>% table() %>% as.data.frame()
   
   if(scaling == 'up'){ 
      targetFreq <- classFreq$Freq %>% max() 
   } else if(scaling == 'down') {
      targetFreq <- classFreq$Freq %>% min() 
   }
   
   targetFreqClass <- classFreq[classFreq$Freq == targetFreq,'.']
   
   df_BalancedClasses <- lapply( as.character(classFreq$.), function(i){
      df_ss <- df[df[,colname_response] == i,] 
      
      if(scaling == 'up'){
         if(i != targetFreqClass){ 
            sample(1:nrow(df_ss), targetFreq, replace = T) %>% df_ss[.,]
         } else {
            df_ss
         }
      
      } else if(scaling == 'down'){
         sample(1:nrow(df_ss), targetFreq, replace = F) %>% df_ss[.,]
      }
      
   }) %>% do.call(rbind,.) %>% .[order(rownames(.)),]

   return(df_BalancedClasses)
}

#--------- Random forest training and validation ---------# 
## Train RF on a train set
## Predict on a test set
## Outputs list(RF object, predictions)

## Depends on: balanceClasses()

randomForest_trainAndTest <- function(df_train, df_test, colname_response, balance=F, bootstrap=F,inclExpectedResponse=T,
                                      ntreeTry=500, stepFactor=1.2, improve=0.001, plot=F, trace=F, ## tuneRF() default args
                                      randomForest_ntree = 500, importance = T, ## randomForest() default args
                                      ...)
{
   #--------- Up/down balance classes ---------#
   if(balance == T | balance == 'up'){
      df_train <- balanceClasses(df_train, colname_response, scaling = 'up')
   } else if (balance == 'down'){
      df_train <- balanceClasses(df_train, colname_response, scaling = 'down')
   }
   
   #--------- bootstrap aggregation; prevents overfitting of model ---------#
   if(bootstrap == T){
      inBag_samples <- sample(1:nrow(df_train), (nrow(df_train) * 0.63212056) %>% round(0) ) ## 0.632 rule; choose ~2/3 of data set for bagging
      
      df_train <- df_train[inBag_samples,]
      df_train <- df_train[rownames(df_train) %>% order,]
   }
   
   #--------- Get mtry where OOBE is min ---------#
   mtryTune <- tuneRF(x = df_train %>% .[,colnames(.) != colname_response], #df of features/observations
                      y = df_train %>% .[,colname_response], ## vector of expected response
                      ntreeTry=ntreeTry, ## number of trees to create
                      stepFactor=stepFactor, ## inflation rate of mtry for each iteration
                      improve=improve, ## relative improvement in OOB error must be by this much for the search to continue
                      plot=plot, 
                      trace=trace,
                      ...)
   
   mtryBest <- 
      mtryTune %>% 
      .[.[,2] == min(.[,2]),1] %>% 
      .[length(.)] ## always select the highest mtry
   
   #--------- Fit RF model ---------#
   RF <- randomForest(x = df_train %>% .[,colnames(.) != colname_response],
                      y = df_train %>% .[,colname_response],
                      ntree = randomForest_ntree,
                      importance = importance,
                      mtry = mtryBest,
                      ...)
   
   #--------- Predict ---------#
   pred <- predict(object = RF, newdata = df_test[,colnames(df_test) != 'response'], type = "prob")
   pred <- pred %>% as.data.frame()
   
   if(inclExpectedResponse == T){
      
      pred$response <- df_test[,colname_response]
   }
   
   #--------- Return RF and prediction object ---------#
   RF_CV <- list()
   RF_CV$RF <-  RF
   RF_CV$pred <- pred
   return(RF_CV)
}

#--------- Random forest k-fold cross validation ---------#
## Depends on: get_cvTrainTestSets()
randomForestCV <- function(df,colname_response,k=10,...)
{
   l_TrainTestSets <- get_cvTrainTestSets(df,k=k)
   
   l_RF_CV <- lapply(1:k, function(i){
      message(paste0('< CV: round ', i,' >'))
      df_train <- l_TrainTestSets[[i]]$train
      df_test <- l_TrainTestSets[[i]]$test
      
      randomForest_trainAndTest(df_train, df_test, colname_response,...)
   })
   
   return(l_RF_CV)
}

#--------- Extract information from randomForest_CV ---------#
aggregate_rfCV <- function(rfCV_object, var)
{
   ## MeanDecreaseAccuracy
   if(var == 'importance'){
      output <- lapply(rfCV_object, function(i){ i$RF %>% importance(.,type=1) }) %>% do.call(cbind,.) %>% apply(.,1,mean)
      output <- data.frame(feature = names(output), MeanDecreaseAccuracy=output)
      rownames(output) <- NULL
   }
   
   ## Prediction probabilities
   else if (var == 'prediction'){
      output <- lapply(rfCV_object, function(i){ i$pred }) %>% do.call(rbind,.) %>% .[order(rownames(.)),]
   }
   
   return(output)
}

#--------- Plot MeanDecreaseAccuracy from aggregate_rfCV() output ---------#
plot_MeanDecreaseAccuracy <- function(df_importance, title='Feature importance', sort=T, topn=10)
{
   #df_importance <- agg_imp
   
   ## Sort features by MeanDecreaseAccuracy
   if(sort == T){
      df_importance <- df_importance %>% .[order(.$MeanDecreaseAccuracy, decreasing = T),]
   }
   
   ## Take top n
   if(topn != F){
      df_importance <- df_importance[1:topn,]
   }
   
   ## Reverse MeanDecreaseAccuracy order again. ggplot2 coord_flip orders the MeanDecreaseAccuracy in the reverse order
   df_importance <- df_importance %>% .[order(.$MeanDecreaseAccuracy, decreasing = F),]
   
   ## Plot id as x to ensure proper sorting order 
   df_importance$id <- 1:nrow(df_importance)
   
   MeanDecAccPlot <- ggplot(df_importance,aes(id,MeanDecreaseAccuracy)) + 
      geom_point(shape=1,size=2) + 
      coord_flip() +
      scale_x_discrete(labels=df_importance$feature,limits=1:length(df_importance$feature)) +
      
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size=10),
            plot.title = element_text(hjust=0.5))
   
   if(title == 'Feature importance'){
      MeanDecAccPlot <- MeanDecAccPlot + ggtitle('Feature importance')
   } else {
      MeanDecAccPlot <- MeanDecAccPlot + ggtitle(title)
   }
}

#========= Stats functions =========#
#--------- Get true/false positive/negative rates ---------#
perfAsDf <- function(v_prob_predicted, v_logical_expected, metrics=c('tpr','tnr','fpr','fnr'))
{
   
   # v_prob_predicted <- aggPred[,'BRCA']
   # v_logical_expected <- df$response == 'BRCA'
   
   ROCRPred <- ROCR::prediction(v_prob_predicted, v_logical_expected)
   
   df_tfpnr <- lapply(metrics,function(i){
      df_metric <- cbind(
         performance(ROCRPred, measure = i)@x.values %>% unlist(),
         performance(ROCRPred, measure = i)@y.values %>% unlist()
      ) %>% as.data.frame()
      colnames(df_metric) <- c('x','y')
      
      df_metric$metric <- i
      
      return(df_metric)
   }) %>% do.call(rbind,.)
   
   return(df_tfpnr)
}

#--------- Plot true positive/negative rate ---------#
## Depends on: perfAsDf()
plot_tfpnRates <- function(v_prob_predicted=NULL, v_logical_expected=NULL, df_melt=NULL,
                           metrics=c('tpr','tnr'), title='True positive vs. true negative rate', showIntersection=T)
{
   ## metrics check
   if(!all(metrics %in% c('tpr','tnr','fpr','fnr'))){
      stop('Metrics must be: tpr, tnr, fpr, fnr')
   }
   
   ## By default, expect v_logical_expected and v_prob_predicted
   ## Also accepts df output from get_tfpnRates
   if( !(is.null(v_prob_predicted)) & !(is.null(v_logical_expected)) == T ){
      df_melt <- perfAsDf(v_prob_predicted, v_logical_expected)
   }
   
   df_melt <- df_melt[df_melt$metric %in% metrics,]
   
   tfpnRatesPlot <- ggplot(data=df_melt, aes(x=x, y=y, colour=metric)) +
      geom_line() +
      
      ggtitle(title) +
      xlab('Probability cutoff') +
      ylab('Fractional value') +
      scale_color_discrete(name='',labels = metrics) +
      
      theme(plot.title = element_text(hjust = 0.5))
   
   if(title == 'True positive vs. true negative rate'){
      tfpnRatesPlot <- tfpnRatesPlot + ggtitle('True positive vs. true negative rate')
   } else {
      tfpnRatesPlot <- tfpnRatesPlot + ggtitle(title)
   }
   
   if(showIntersection == T){
      
      diff_tpnRate <- abs( df_melt[df_melt$metric == 'tpr','y'] - df_melt[df_melt$metric == 'tnr','y'] )
      intersection <- which(diff_tpnRate == min(diff_tpnRate)) %>% df_melt[.,'x'] %>% signif(.,4)
      
      tfpnRatesPlot <- tfpnRatesPlot +
         geom_vline(xintercept = intersection,linetype = 3,colour = 'black') +
         annotate('text', x=intersection*0.95, y=0.07, hjust = 1, vjust = 0.5, label=paste0('P = ',intersection))
   }
   
   return(tfpnRatesPlot)
}

#--------- Plot ROC curve ---------#
## Depends on: perfAsDf()
plot_ROC <- function(v_prob_predicted=NULL, v_logical_expected=NULL, title='ROC', showAUC=T)
{
   df_melt <- perfAsDf(v_prob_predicted, v_logical_expected, metrics=c('tpr','fpr'))
   
   df_tpr_fpr <- cbind(
      df_melt[df_melt$metric == 'fpr','y'],
      df_melt[df_melt$metric == 'tpr','y']
   ) %>% as.data.frame()
   colnames(df_tpr_fpr) <- c('x','y')
   
   ROCPlot <- ggplot(data=df_tpr_fpr, aes(x=x, y=y)) +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype=3) + 
      
      xlab('False positive rate') +
      ylab('True positive rate') +
      theme(plot.title = element_text(hjust = 0.5))
   
   if(title == 'ROC'){
      ROCPlot <- ROCPlot + ggtitle('ROC')
   } else {
      ROCPlot <- ROCPlot + ggtitle(title)
   }
   
   if(showAUC == T){
      auc <- prediction(v_prob_predicted, v_logical_expected) %>% performance(.,'auc') %>% .@y.values %>% unlist() 
      ROCPlot <- ROCPlot + annotate('text', x=0.5, y=min(df_tpr_fpr$y), hjust = 0.5, vjust = 0.5, label=paste0('AUROC = ', auc %>% round(.,4)))
   }
   
   return(ROCPlot)
}

#--------- Plot precision recall ---------#
## Precision/Positive predictive value: probability that subjects with a positive screening test truly have the disease
## Recall/True positive rate
calc_AUPRC <- function(v_prob_predicted, v_logical_expected)
{
   df_melt <- perfAsDf(v_prob_predicted, v_logical_expected, metrics=c('ppv','tpr'))
   
   ## Prepare dataframe
   df_ppv_tpr <- cbind(
      df_melt[df_melt$metric == 'tpr','y'],
      df_melt[df_melt$metric == 'ppv','y']
   ) %>% .[,] %>% as.data.frame()
   colnames(df_ppv_tpr) <- c('x','y')
   
   ## Assign point (0,1) as first point
   df_ppv_tpr[df_ppv_tpr$x==0,]$y <- 1
   
   ## AUC calculation
   df_ppv_tpr_aucCalc <- df_ppv_tpr[df_ppv_tpr$x != 1,] %>% rbind(.,c(1,0)) ## Set (1,0) as last coordinate (set last y value to 0)
   auc <- integrateDiscrete(df_ppv_tpr_aucCalc) ## No need to divide by max(x)*max(y); lims are already (1,1)
   return(auc)
}

plot_PRC <- function(v_prob_predicted=NULL, v_logical_expected=NULL, title='PRC', showAUC=T)
{
   #df_melt <- perfAsDf(agg_pred$BRCA,agg_pred$response, metrics=c('ppv','tpr'))
   #title='test'
   
   df_melt <- perfAsDf(v_prob_predicted, v_logical_expected, metrics=c('ppv','tpr'))
   
   ## Prepare dataframe
   df_ppv_tpr <- cbind(
      df_melt[df_melt$metric == 'tpr','y'],
      df_melt[df_melt$metric == 'ppv','y']
   ) %>% .[,] %>% as.data.frame()
   colnames(df_ppv_tpr) <- c('x','y')
   
   ## Assign point (0,1) as first point
   df_ppv_tpr[df_ppv_tpr$x==0,]$y <- 1
   
   PRCPlot <- ggplot(data=df_ppv_tpr, aes(x=x, y=y)) +
      geom_line() +
      geom_hline(yintercept = 0.5, linetype=3) +
      
      ggtitle(title) +
      ylab('Positive predictive value') +
      xlab('True positive rate') +
      
      theme(plot.title = element_text(hjust = 0.5))
   
   if(title == 'PRC'){
      PRCPlot <- PRCPlot + ggtitle('PRC')
   } else {
      PRCPlot <- PRCPlot + ggtitle(title)
   }
   
   ## Calculate AUC
   if(showAUC == T){
      df_ppv_tpr_aucCalc <- df_ppv_tpr[df_ppv_tpr$x != 1,] %>% rbind(.,c(1,0)) ## Set (1,0) as last coordinate (set last y value to 0)
      
      auc <- integrateDiscrete(df_ppv_tpr_aucCalc) ## No need to divide by max(x)*max(y); lims are already (1,1)
      PRCPlot <- PRCPlot + annotate('text', x=0.5, y=min(df_ppv_tpr$y), hjust = 0.5, vjust = 0.5, label=paste0('AUPRC = ', auc %>% round(.,4)))
   }
   
   return(PRCPlot)
}

#--------- Plot feature AUPRC or AUROC ---------#
plot_AUCFeatSel <- function(df_summaryFeatSel, metric='auprc')
{
   if(!(metric %in% c('auprc','auroc'))){
      stop('Metric must be: auprc, auroc')
   }
   
   AUCFeatSelPlot <- ggplot(df_summaryFeatSel, aes_string(x='iter', y=metric, label='feature')) +
      geom_point(shape=1) +
      stat_smooth(method = 'loess', span = 1, fill = NA, size = 0.5, colour = 'red') +
      
      xlab('Minimum meanDecreaseAccuracy feature |  Iteration') +
      ylab(toupper(metric)) +
      
      ggtitle('Iterative feature selection by MeanDecreaseAccuracy') +
      
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text=element_text(size=11),
            axis.title=element_text(size=12)) +
      
      scale_x_continuous(breaks = df_summaryFeatSel$iter,
                         labels = paste0(df_summaryFeatSel$feature, ' | ',
                                         str_pad(df_summaryFeatSel$iter,2,pad='0') ) )
   
   return(AUCFeatSelPlot)
}

####################################
#          Functions: End          #
####################################

#========= Paths =========#
localPathPrefix <- '/Users/lnguyen/'
hpcPath <- '/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/random_forest_training/feature_selection_umcuData_whitelist_contextMatrix/'
if( dir.exists(localPathPrefix) ){
   base_dir <- paste0(localPathPrefix, hpcPath)
} else {
   base_dir <- hpcPath
}

#========= Donor white list =========#
## donors that are CERTAIN to be BRCA proficient or deficient
donorWhitelist <- read.table( '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/ML_data/donorWhitelist.tsv', sep='\t',header=T)
donorWhitelistNames <- rownames(donorWhitelist)

#========= Input: UMCU data =========#
## Norm signature data
df_umcuNormData <- read.table( paste0(base_dir,'df_umcuNormData_BRCAsimple.tsv'),sep='\t',header=T)
df_umcuNormData$response <- as.factor(df_umcuNormData$response)

df_umcuNormData_whitelist <- df_umcuNormData %>% .[rownames(.) %in% donorWhitelistNames,]

## Context data
df_umcuDataContext <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/ML_data/df_umcuDataContext_BRCAannot.tsv',sep='\t',header=T,check.names=F)

df_umcuDataContext$response <- df_umcuDataContext$response %>% get_simpleResponse(.,c('BRCA1','BRCA2'),'BRCA','none','none')
df_umcuDataContext_whitelist <- df_umcuDataContext %>% .[rownames(.) %in% donorWhitelistNames,]

#========= Set input =========#
#df <- df_umcuNormData_whitelist
df <- df_umcuDataContext_whitelist

#========= CV and plots =========#
rfCV <- randomForestCV(df, colname_response = 'response', inclExpectedResponse = T)

## Aggregate information
agg_pred <- aggregate_rfCV(rfCV, 'prediction')
agg_pred$response <- get_simpleResponse(agg_pred$response,'BRCA',1,'none',0)

agg_imp <- aggregate_rfCV(rfCV, 'importance')

## Plotting
rfCV_summaryPlots <- list(
   plot_MeanDecreaseAccuracy(agg_imp,topn = 50),
   plot_tfpnRates(agg_pred$BRCA,agg_pred$response),
   plot_ROC(agg_pred$BRCA,agg_pred$response),
   plot_PRC(agg_pred$BRCA,agg_pred$response)
)

## Export
pdf(paste0(base_dir,'rfCV_summaryPlots.pdf'),10,10)
grid.arrange(grobs = rfCV_summaryPlots, nrow = 2, ncol = 2)
dev.off()

#========= Feature selection =========#
nFeatures <- ncol(df)-1
featSelIters <- 1:(nFeatures-1) ## End when only 1 feature(s) are left

f_report <- paste0(base_dir,'summaryFeatureSelection.tsv')
write.table(paste('iter','minImpVar','feature','auroc','auprc','optimalCutOff',sep = '\t'), f_report, col.names = F, row.names = F, quote = F)

v_minImpVar <- c()
for(iter in featSelIters){
   message(paste0('<< Feature selection: round ', iter,' >>'))
   iterName <- str_pad(iter, 2, pad = "0")
   
   df_ss <- df[,!(colnames(df) %in% names(v_minImpVar))]
   rfCV <- randomForestCV(df_ss, colname_response = 'response', inclExpectedResponse = T)
   
   ## Aggregate information
   message('Aggregating CV information...')
   agg_pred <- aggregate_rfCV(rfCV, 'prediction')
   agg_pred$response <- get_simpleResponse(agg_pred$response,'BRCA',1,'none',0)
   
   agg_imp <- aggregate_rfCV(rfCV, 'importance')
   
   ## Plotting
   message('Making summary plots...')
   rfCV_summaryPlots <- list(
      plot_MeanDecreaseAccuracy(agg_imp, topn = F),
      plot_tfpnRates(agg_pred$BRCA, agg_pred$response),
      plot_ROC(agg_pred$BRCA, agg_pred$response),
      plot_PRC(agg_pred$BRCA, agg_pred$response)
   )
   
   pdf(paste0(base_dir,'summaryPlotsFeatSel_iter',iterName,'.pdf'),10,10)
   grid.arrange(grobs = rfCV_summaryPlots, nrow = 2, ncol = 2)
   dev.off()
   
   ## Append to summary table
   message('Writing to summary table...')
   
   df_minImpVar <- agg_imp %>% .[.$MeanDecreaseAccuracy == min(.$MeanDecreaseAccuracy),]
   minImpVar <- df_minImpVar$MeanDecreaseAccuracy
   minImpVarFeature <- df_minImpVar$feature
   
   auroc <- prediction(agg_pred$BRCA, agg_pred$response) %>% performance(.,'auc') %>% .@y.values %>% unlist()
   auprc <- calc_AUPRC(agg_pred$BRCA, agg_pred$response)
   
   df_tfpnRates <- perfAsDf(agg_pred$BRCA, agg_pred$response)
   diff_tpnRate <- abs( df_tfpnRates[df_tfpnRates$metric == 'tpr','y'] - df_tfpnRates[df_tfpnRates$metric == 'tnr','y'] )
   intersection <- which(diff_tpnRate == min(diff_tpnRate)) %>% df_tfpnRates[.,'x']
   
   write(paste(iter,minImpVar,minImpVarFeature,auroc,auprc,intersection,sep='\t'), f_report, append = T)
   
   ## Append to feature exclusion vector
   names(minImpVar) <- minImpVarFeature
   v_minImpVar <- c(v_minImpVar,minImpVar)
}

#========= Plot feature selection summary =========#
df_summaryFeatSel <- read.table(paste0(base_dir,'summaryFeatureSelection.tsv'), header = T, sep = '\t')

pdf(paste0(base_dir,'AUCvsFeatSel.pdf'),9,4)
plot_AUCFeatSel(df_summaryFeatSel, metric = 'auroc')
plot_AUCFeatSel(df_summaryFeatSel, metric = 'auprc')
dev.off()

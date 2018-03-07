#========= Load packages =========#
packages <- c('glmnet',
              'ggplot2',
              'magrittr',
              'reshape2',
              'ROCR')

ipak <- function(pkg){
   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
   if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
   sapply(pkg, require, character.only = TRUE)
}

ipak(packages)

#========= Global options =========#
options(stringsAsFactors = F)

#========= Misc. functions =========#
#--------- Normalize data using standard deviation ---------#
normalizeVector <- function(v, na.replace = T)
{
   ln_v <- log(v+1, exp(1))
   norm_v <- (ln_v - mean(ln_v, na.rm = T) ) / sd(ln_v, na.rm = T)

   if(na.replace == T){
      norm_v[is.na(norm_v)] <- median(norm_v, na.rm = T) ## replace NA with median
   }

   return(norm_v)
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

#--------- Remove donor suffix ---------#
rmDonorSuffix <- function(v)
{
   v %>% str_extract_all(.,'PD\\d+') %>% unlist()
}

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
#========= Logistic regression functions =========#
#--------- Balance classes ---------#
balanceClasses <- function(df, colnameResponse, scaling = 'up')
{
   # df <- df_umcuNormData
   
   classFreq <- df[,colnameResponse] %>% table() %>% as.data.frame()
   
   if(scaling == 'up'){ 
      targetFreq <- classFreq$Freq %>% max() 
   } else if(scaling == 'down') {
      targetFreq <- classFreq$Freq %>% min() 
   }
   
   targetFreqClass <- classFreq[classFreq$Freq == targetFreq,'.']
   
   df_BalancedClasses <- lapply( as.character(classFreq$.), function(i){
      df_ss <- df[df[,colnameResponse] == i,] 
      
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

#--------- Logistic regression function ---------#
logisticRegression <- function(v_featureValues, intercept, v_coef)
{
   P <- 1 / (1 + exp(1) ^ -(intercept + sum(v_featureValues * v_coef) ) )
   return( P %>% unname() )
}

#--------- Nested cross validation ---------#
logRegCVnested <- function(df = NULL, colnameResponse = NULL, stratify = T, kOuter=10, balance = T,
                           kInner = 10, alpha = 1, type.measure = 'class', nlambda = 100, lower.limits = 0,
                           predictLambda = 'lambda.min', ...)
{
   ## Get outer fold train and test sets
   if(stratify == T){
      l_TrainTestSets <- get_cvTrainTestSets(df, k=kOuter, stratifyByCol = colnameResponse)
   } else {
      l_TrainTestSets <- get_cvTrainTestSets(df, k=kOuter)
   }
   
   ## Nested cross validation
   LR_CV_nested <- lapply(1:length(l_TrainTestSets), function(outerFold)
   {
      message(paste0('<< Outer fold CV: ', outerFold, ' >>'))
      
      ## Get train and test for current outer fold
      df_train <- l_TrainTestSets[[outerFold]]$train
      df_test <- l_TrainTestSets[[outerFold]]$test
      
      ## Balance classes
      if(balance == T | balance == 'up'){
         df_train <- balanceClasses(df_train, colnameResponse, scaling = 'up')
      } else if (balance == 'down'){
         df_train <- balanceClasses(df_train, colname_response, scaling = 'down')
      }
      
      ## Inner fold CV: Determining best lambda
      LR_CV_inner <- 
         cv.glmnet(x = as.matrix(df_train %>% .[, colnames(.) != colnameResponse ]), ## signature matrix
                   y = df_train[,'response'], ## response vector
                   family = 'binomial',
                   nfolds = kInner,
                   alpha = alpha,
                   type.measure = type.measure,
                   nlambda = nlambda,
                   lower.limits = lower.limits,
                   ...)
      
      ## Predict on test set with best lambda
      pred <- 
         predict.cv.glmnet(LR_CV_inner,
                           as.matrix(df_test %>% .[, colnames(.) != colnameResponse ]),
                           s = predictLambda,
                           type = 'response'
         )
      
      colnames(pred) <- 'prediction'
      pred <- pred %>% as.data.frame()
      pred$response <- df_test$response
      
      ## Store CV object and prediction probabilities
      LR_summary <- list(LR_CV_inner = LR_CV_inner, pred = pred)
      return(LR_summary)
   })
   
   return(LR_CV_nested)
}

#--------- Summarize nested cross validation ---------#
aggregate_logRegCVnested <- function(logRegCVnested_object = NULL, var = NULL)
{
   ## lambda
   if(var == 'lambda.1se'){
      output <- sapply(logRegCVnested_object, function(fold){ fold$LR_CV_inner[['lambda.1se']] })
   } else if (var == 'lambda.min'){
      output <- sapply(logRegCVnested_object, function(fold){ fold$LR_CV_inner[['lambda.min']] })
   }
   
   ## coef
   else if(var == 'coef'){
      output <- lapply(logRegCVnested_object, function(fold){ fold$LR_CV_inner %>% coef() }) %>% do.call(cbind,.)
   }
   
   ## pred
   else if(var == 'pred'){
      output <- lapply(logRegCVnested_object, function(outerFold){ outerFold$pred }) %>% do.call(rbind,.) %>% .[order(rownames(.)),]
   }
   
   return(output)
}

#--------- Plot weights from logistic regression CV ---------#
plot_logRegWeights <- function(df, title = 'Feature weights')
{
   ## Convert to df
   df <- df %>% .[apply(.,1,mean) != 0,] %>% as.matrix() %>% as.data.frame() 
   
   ## Rename rows and cols
   colnames(df) <- 1:ncol(df)
   rownames(df) <- rownames(df) %>% gsub('[()]','',.)
   
   ## Coef means
   mean_coef <- apply(df,1,mean) %>% .[. != 0] %>% sort(.,decreasing = F)
   mean_coef <- c(mean_coef %>% .[names(.) == 'Intercept'], mean_coef %>% .[names(.) != 'Intercept'])
   
   ## Reorder
   df <- df %>% .[names(mean_coef),]
   
   df$feature <- rownames(df)

   df$id <- 
      1:nrow(df) %>% 
      as.character() %>% 
      str_pad(., width = nchar(.) %>% max, pad = '0',side = 'left')
   df_melt <- melt(df, c('id','feature'))
   
   ## Annotation of mean
   mean_label <- signif(mean_coef,3) %>% as.character()
   value_range <- max(df_melt$value) - min(df_melt$value)
   ylim_min <- min(df_melt$value) - value_range*0.18
   ylim_max <- max(df_melt$value) + value_range*0.1
   
   ## Plot
   weightsPlot <- ggplot(df_melt, aes(x=id,y=value)) + 
      
      geom_hline(yintercept = 0, color = 'dark gray') +
      geom_boxplot() + 
      geom_jitter() +
      
      stat_summary(fun.y="mean", geom = "point", shape=4, size=5, color="red") +
      annotate("text", x = 1:nrow(df), y = ylim_min, label = mean_label, size = 3.2, hjust = 0, color = "red4") +
      
      scale_x_discrete(labels = df$feature) +
      scale_y_continuous(expand = c(0, 0), limits = c(ylim_min - value_range*0.01, ylim_max) ) +
      
      xlab('Feature') +
      ylab('Coefficient') +
      
      theme(
         axis.ticks.y = element_blank(),
         panel.background = element_blank(),
         panel.border = element_rect(fill = NA),
         plot.title = element_text(hjust = 0.5)
      ) +
      
      coord_flip()
   
   if(title == 'Feature weights'){
      weightsPlot <- weightsPlot + ggtitle('Feature weights')
   } else {
      weightsPlot <- weightsPlot + ggtitle(title)
   }
   
   return(weightsPlot)
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

#========= Working dir =========#
localPath <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/umcuDataContextRel/'
if( dir.exists(localPath) ){
   base_dir <- localPath
} else {
   base_dir <- sub('/Users/lnguyen/','/',localPath)
}

#========= Donor white list =========#
## donors that are CERTAIN to be BRCA proficient or deficient
donorWhitelist <- read.table( '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/ML_data/donorWhitelist.tsv', sep='\t',header=T)
donorWhitelistNames <- rownames(donorWhitelist)

#========= Input: UMCU data =========#
## Norm signature data
df_umcuNormData <- read.table( '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/ML_data/df_umcuNormData_BRCAannot.tsv',sep='\t',header=T)
df_umcuNormData$response <- df_umcuNormData$response %>% get_simpleResponse(.,c('BRCA1','BRCA2'),1,'none',0)

df_umcuNormData_whitelist <- df_umcuNormData %>% .[rownames(.) %in% donorWhitelistNames,]

#========= Input: UMCU data context matrices =========#
df_umcuDataContextRel <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/ML_data/df_umcuDataContextRel_BRCAannot.tsv',sep='\t',header=T,check.names=F)

df_umcuDataContextRel$response <- df_umcuDataContextRel$response %>% get_simpleResponse(.,c('BRCA1','BRCA2'),1,'none',0)
df_umcuDataContextRel_whitelist <- df_umcuDataContextRel %>% .[rownames(.) %in% donorWhitelistNames,]

## Remove 'all zero' columns and normalize
df_umcuDataContextRelNorm_whitelist <- 
   df_umcuDataContextRel_whitelist %>% 
   .[,colnames(.) != 'response'] %>% 
   .[,apply(.,2,sum) !=0] %>% ## rm zero columns
   apply(.,2,normalizeVector) %>% ## normalize
   as.data.frame()
df_umcuDataContextRelNorm_whitelist$response <- df_umcuDataContextRel_whitelist$response

# #========= Input: Sanger data =========#
# #------ Self normalized ------#
# df_sangerRawData <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/nm.4292-S2_raw_data.txt',
#                       sep = '\t', header = T)
# 
# rownames(df_sangerRawData) <- df_sangerRawData$Sample
# 
# ## only take feature columns
# df_sangerNormData <- df_sangerRawData[,9:33]
# rownames(df_sangerNormData) <- rownames(df_sangerNormData) %>% rmDonorSuffix()
# 
# ## Normalize data
# df_sangerNormData <- df_sangerNormData %>% apply(.,2,normalizeVector) %>% as.data.frame()
# 
# ## Attach brca annotation
# df_brcaAnnotation <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/brca_annotation_simple.tsv')
# colnames(df_brcaAnnotation) <- c('donor','response')
# 
# df_brcaAnnotation$donor <- df_brcaAnnotation$donor %>% rmDonorSuffix()
# df_sangerNormData$response <- df_brcaAnnotation$response
# 
# ## Training data set:
# df_sangerNormData_whitelist <- df_sangerNormData %>% .[rownames(.) %in% donorWhitelistNames,]

#========= Logistic regression training =========#
df <- df_umcuNormData_whitelist

# df <- df_umcuDataContextRel_whitelist
# df <- df_umcuDataContextRelNorm_whitelist

lrCV <- logRegCVnested(df, 'response')

#========= Summarize training =========#
#agg_lambda <- aggregate_logRegCVnested(lrCV,'lambda.min')
agg_coef <- aggregate_logRegCVnested(lrCV,'coef')
agg_pred <- aggregate_logRegCVnested(lrCV,'pred')

lrCV_summaryPlots <- list(
   plot_logRegWeights(agg_coef),
   plot_tfpnRates(agg_pred$prediction, agg_pred$response),
   plot_ROC(agg_pred$prediction, agg_pred$response),
   plot_PRC(agg_pred$prediction, agg_pred$response)
)

pdf(paste0(base_dir,'lrCV_summaryPlots.pdf'),10,10)
grid.arrange(grobs = lrCV_summaryPlots, nrow = 2, ncol = 2)
dev.off()


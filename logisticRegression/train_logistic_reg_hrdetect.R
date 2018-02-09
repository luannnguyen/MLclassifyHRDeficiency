packages <- c('glmnet',
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


#========= Data prep =========#
rawdata <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/nm.4292-S2_raw_data.txt',
                      sep = '\t', header = T)
rawdata_ss <- rawdata[,9:33]
rownames(rawdata_ss) <- rawdata$Sample
rownames(rawdata_ss) <- 
   rownames(rawdata_ss) %>% 
   gsub('a2','',.) %>% 
   gsub('a','',.) %>%
   gsub('b','',.)

normalize_data <- function(df, na.replace = T){
   df_norm <- apply(df,2,function(v){
      ln_v <- log(v+1, exp(1))
      norm_v <- (ln_v - mean(ln_v, na.rm = T) ) / sd(ln_v, na.rm = T)
      
      if(na.replace == T){
         norm_v[is.na(norm_v)] <- median(norm_v, na.rm = T) ## replace NA with median
      }
      
      return(norm_v)
   })
   return( as.data.frame(df_norm) )
}

normdata <- normalize_data(rawdata_ss)

## Load brca status annotation
df_brca_annotation <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/brca_annotation_simple.tsv')
colnames(df_brca_annotation) <- c('donor','brca_deficiency')

df_brca_annotation$donor <- 
   df_brca_annotation$donor %>% 
   gsub('a2','',.) %>% 
   gsub('a','',.) %>%
   gsub('b','',.)

## Training data set
normdata_ss <- normdata[,1:22] ## only take feature columns
normdata_ss$brca_deficiency <- as.factor(df_brca_annotation$brca_deficiency)

#========= Train logistic regression by cross-validation with multiple iterations =========#
## Per iteration:
## - determine optimal lambda (i.e. regularization parameter. Need to find balanced lamdba; i.e. not underfit, but also not over fit)
## - determine coefficients with the corresponding lambda

cv_iterations <- 10

l_CV <- lapply(1:cv_iterations, function(i){ #do cross-validation 10 times
   CV <- cv.glmnet(x = as.matrix(normdata_ss[,1:22]), #signature matrix
                   y = as.factor(normdata_ss$brca_deficiency), #brca status
                   family = 'binomial', #i.e. brca status is either 0 or 1
                   type.measure = "class", #cross validate by misclassification error (i.e. wrong prediction of brca status)
                   alpha = 1, #alpha = 1: only lasso regularization, alpha = 0: only ridge regularization
                   nlambda = 100, #number of lambda values to try
                   nfolds = 10,
                   lower.limits =  rep(0,ncol(normdata_ss)) #force no negative coefs
                   )
   return(CV)
})

# ## Plot log(lambda) vs missclassification error for each iteration
# # 1. left dotted line: minimum lambda (i.e lambda for overfitted model)
# # 2. right dotted line: lambda that is within 1 standard error of minimum lambda (i.e. lambda for regularized model; the lambda value that is chosen)
# pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/logLambda_vs_misclassError.pdf', 5*(16/9), 5)
# for(i in l_CV){
#    plot(i, main = paste0('Lambda 1SE: ', i$lambda.1se))
# }
# dev.off()

#plot(l_CV[[1]])

## Get coefs from cross-validation iterations
l_coefs <- lapply(l_CV, function(i){
   #coef(i,s = "lambda.min") %>% as.matrix()
   coef(i, s = "lambda.1se") %>% as.matrix()
})

df_coefs <- do.call(cbind, l_coefs)
colnames(df_coefs) <- paste0('CV.', 1:length( colnames(df_coefs) ))
rownames(df_coefs)[1] <- 'intercept'

df_coefs

predict.cv.glmnet(l_CV[[1]], as.matrix(normdata_ss),s="lambda.1se",type='response') %>% head 

coef(CV, s = "lambda.min")
coef(CV, s = "lambda.1se")

#--------- Boxplot of coefs ---------#
df_coefs_ss <- df_coefs[apply(df_coefs,1,sum) != 0,] ## remove features where coefs always equals 0; cleaner boxplot

df_coefs_ss_melt <- melt(df_coefs_ss)
colnames(df_coefs_ss_melt) <- c('feature','CV.iteration','weight')

CV_iter_coefs <- ggplot(df_coefs_ss_melt, aes(x = feature, y = weight)) +
   geom_boxplot() +
   xlab('') +
   ylab('Coefficient value') +
   ggtitle(paste0('Distribution of coefficients for ',cv_iterations,' cross-validation iterations'))

pdf('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/CV_iter_coefs_nonzero.pdf', 8, 8)
plot(CV_iter_coefs)
dev.off()

## Calculate mean per coefficient
weight_vector_trained <- apply(df_coefs,1,mean)

#--------- Compare HRDetect weights vs. trained weights ---------#
df_hrdetect_weights <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/nm.4292-S4_hrdetect_weights.txt',
                                  sep='\t', header = T, stringsAsFactors = F)

df_weights_hrdetect_vs_trained <- cbind(df_hrdetect_weights$Weight,weight_vector_trained[df_hrdetect_weights$Feature_label])
colnames(df_weights_hrdetect_vs_trained) <- c('hrdetect','trained')

write.table(df_weights_hrdetect_vs_trained,'/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/weights_hrdetect_vs_trained.tsv',
            sep = '\t', quote = F, col.names = F)

#--------- Run BRCA predictor function on trained weights ---------#
# p_brca <- function(v_features_donor, intercept, weight_vector){
#    p <- 1 / (1 + exp(1) ^ -(intercept + sum(v_features_donor * weight_vector) ) )
#    return(p)
# }
# Pbrca_sanger_trained <- apply(normdata_ss[,1:22],1, function(v){ p_brca(v,weight_vector_trained[1],weight_vector_trained[-1]) })

#Pbrca_sanger_trained <- predict.cv.glmnet(CV, as.matrix(normdata_ss),s=0.00369,type='response') 

#--------- Comparing pipelines ---------#
annotated_donors_ss <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/plots_compare_hrdetect_pipelines/annotated_donors_ss_multi_datasets.tsv',
                                  sep = '\t')

annotated_donors_ss$Pbrca_sanger_trained <- Pbrca_sanger_trained[rownames(annotated_donors_ss)]
annotated_donors_ss$Pbrca_sanger_trained_gt_0.7 <- as.integer(annotated_donors_ss$Pbrca_sanger_trained > 0.7)

annotated_donors_ss <- annotated_donors_ss[,c(1:5,10,6:9,11)]
annotated_donors_ss$Pbrca_sanger_trained <- annotated_donors_ss$Pbrca_sanger_trained %>% format(scientific = T)

# cbind(annotated_donors_ss[,c('brca_annotation','predictorProb')],Pbrca_sanger_trained[rownames(annotated_donors_ss)])


# write.table(annotated_donors_ss[c(1,3,8,6,11)],'/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/predictorProb_vs_Pbrca_sanger_trained.tsv',
#             sep = '\t', quote = F)

# write.table(annotated_donors_ss,'/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/comparing_pipelines.tsv',
#             sep = '\t', quote = F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# #========= Load norm data from paper =========#
# normdata_sanger <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/nm.4292-S4_normalized data.txt',
#                               sep = '\t', header = T)
# 
# normdata_sanger$normalised.hrd[is.na(normdata_sanger$normalised.hrd)] <- median(normdata_sanger$normalised.hrd,na.rm = T)
# 
# # apply(normdata_sanger,2,function(col){
# #    col[is.na(col)]
# # })
# 
# # normdata_sanger_ss <- normdata_sanger[,3:15][,-7]
# # 
# # CV <- cv.glmnet(x = as.matrix(normdata_sanger_ss), ## signature matrix
# #                 y = as.factor(df_brca_annotation$brca_deficiency), ## brca status
# #                 family = 'binomial', ## i.e. brca status is either 0 or 1
# #                 type.measure = "class", ## cross validate by misclassification error (i.e. wrong prediction of brca status)
# #                 alpha = 1, ## alpha = 1: only lasso regularization, alpha = 0: only ridge regularization
# #                 nlambda = 100, ## number of lambda values to try
# #                 lower.limits = 0
# # )
# # coef(CV)
# 
# fit <- glmnet(x = as.matrix(normdata_sanger_ss), ## signature matrix
#               y = as.factor(df_brca_annotation$brca_deficiency), ## brca status
#               family = 'binomial', ## i.e. brca status is either 0 or 1
#               alpha = 1#, ## alpha = 1: only lasso regularization, alpha = 0: only ridge regularization
#               #lambda = 0.000480
# )
# coef(fit)

# #========= train on umcu data =========#
# #--------- Data prep ---------#
# ## umcu data
# df_umcu_data_full_norm <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/matrices/signature_matrices/umcu_pipeline_full_signature_matrix_norm.tsv',
#                                      sep='\t')
# 
# ## Load brca status annotation
# df_brca_annotation <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/supplementary_data/brca_annotation_simple.tsv')
# colnames(df_brca_annotation) <- c('donor','brca_deficiency')
# 
# df_brca_annotation$donor <- 
#    df_brca_annotation$donor %>% 
#    gsub('a2','',.) %>% 
#    gsub('a','',.) %>%
#    gsub('b','',.)
# 
# 
# df_brca_annotation_ss <- df_brca_annotation[df_brca_annotation$donor %in% rownames(df_umcu_data_full_norm),] ## umcu data
# 
# ## bind brca annotation to umcu data
# df_umcu_data_full_norm$brca_deficiency <- df_brca_annotation_ss$brca_deficiency ## umcu data
# 
# #--------- Train logistic regression ---------#
# ## cross-validation to determine optimal lambda (i.e. regularization value; coefficient penalty)
# CV <- cv.glmnet(x = as.matrix(df_umcu_data_full_norm[,-46]), ## signature matrix
#                 y = as.factor(df_umcu_data_full_norm$brca_deficiency), ## brca status
#                 family = 'binomial', ## i.e. brca status is either 0 or 1
#                 type.measure = "class", ## cross validate by misclassification error (i.e. wrong prediction of brca status)
#                 alpha = 1, ## alpha = 1: only lasso regularization, alpha = 0: only ridge regularization
#                 nlambda = 100 ## number of lambda values to try
#                 )
# 
# ## plot log(lamdba) vs missclassification error
# # 1. left dotted line: minimum lambda (i.e lambda for overfitted model)
# # 2. right dotted line: determine lambda that is within 1 standard error of minimum lambda (i.e. lamdba for regularized model)
# #plot(CV)
# #CV$lambda.1se ## lamdba for regularized model
# 
# ## get coefficients
# coefs <- coef(CV) %>% as.matrix()
# intercept_trained <- coefs[1,]
# 
# weight_vector_trained <- coefs[-1,] 
# weight_vector_trained[weight_vector_trained < 0] <- 0 ## constrain to positive weights as mentioned in the hrdetect paper; set negative weights to 0
# 
# # weight_vector_trained[weight_vector_trained > 0]
# # df_hrdetect_weights
# 
# # write.table(weight_vector_trained,'/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/weight_vector_trained.tsv',
# #             sep = '\t', quote = F, col.names = F)
# 
# ## run BRCA predictor function
# p_brca <- function(v_features_donor, intercept, weight_vector){
#    p <- 1 / (1 + exp(1) ^ -(intercept + sum(v_features_donor * weight_vector) ) )
#    return(p)
# }
# 
# Pbrca_umcu_trained <- apply(df_umcu_data_full_norm[,-46],1, function(v){ p_brca(v,intercept_trained,weight_vector_trained) })
# 
# # write.table(Pbrca_umcu_trained,'/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/results/hr_detect/logistic_regression_training/Pbrca_umcu_trained.tsv',
# #             sep = '\t', quote = F, col.names = F)
# 

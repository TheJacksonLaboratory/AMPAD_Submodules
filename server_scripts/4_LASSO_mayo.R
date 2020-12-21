# Notes
# * Scripts integrated from work by Cai John
# * Data taken from the Amp-AD Repositories

# Goal
# * Use LASSO to avoid overuse of Eigengenes
# * Run for multiple iterations for further analysis

#------------------------------------------------
# Static
#------------------------------------------------

# Clear environment
rm(list = ls())
setwd("../")

# Define cohort (mayo, mssm, rosmap)
COHORT <- "mayo"

# Number of folds
if (COHORT=="mayo") {
  FOLDS <- 20
} else if (COHORT=="mssm") {
  FOLDS <- 150
} else {
  FOLDS <- 125
}

# Number of models to run in LASSO
NUM.MODELS <- 500
LAMBDA.ITER <- 5
TRAIN.RATIO <- 0.9

#------------------------------------------------
# Libraries and Data
#------------------------------------------------

library(dplyr)
library(ggplot2)
library(glmnet)
library(lmtest)
library(mgsub)
library(ROCR)

submodules <- readRDS(paste0("clean_data/2_build_svd_sets/", COHORT, "_submodules.RDS"))
submodules.labels <- readRDS(paste0("clean_data/2_build_svd_sets/", COHORT, "_submodule_labels.RDS"))
eigengenes <- readRDS(paste0("clean_data/2_build_svd_sets/", COHORT, "_eigengenes.RDS"))

#------------------------------------------------
# Run LASSO Models
#------------------------------------------------

# Order data set
eigengenes <- eigengenes %>%
  as_tibble() %>%
  group_by(Diagnosis)

# Create storage objects
metrics <- as.data.frame(matrix(nrow=NUM.MODELS, ncol=5))
names(metrics) <- c("AUC", "ACC", "LRI", "Lambda", "SD.Lambda")

LASSO.models <- list()

# Iterate to produce models
for (count in 1:NUM.MODELS) {
  
  # Create training and validation sets
  cat(paste0("Computing predictions for model number ", count, "\n"))
  train <- sample_frac(eigengenes, TRAIN.RATIO)
  valid <- eigengenes[which(!(eigengenes$Patient %in% train$Patient)),]
  
  # Set up predictors and response variables
  # Sets for training and validation data
  predictors.train <- as.matrix(train[,4:ncol(train)])
  response.train <- as.factor(as.matrix(train[,3])) %>%
    factor(levels(.)[c(2,1)]) %>%
    as.numeric()
  
  predictors.valid <- as.matrix(valid[,4:ncol(valid)])
  response.valid <- as.factor(as.matrix(valid[,3])) %>%
    factor(levels(.)[c(2,1)]) %>%
    as.numeric()
  
  # Determine robust lambda for multinomial regression
  min.lambda <- data.frame(lambda=rep(NA, LAMBDA.ITER))
  for (i in 1:LAMBDA.ITER){
    logis.model <- cv.glmnet(predictors.train, response.train, family="binomial", type.measure="auc", nfolds=20)
    min.lambda$lambda[i]<- logis.model$lambda.min
  }
  
  # Find average lambda value and corresponding model
  cat("Number of unique Lambda values found: ", length(unique(min.lambda$lambda)), "/", paste0(LAMBDA.ITER),"\n")
  lambda <- mean(min.lambda$lambda)
  logis.model <- glmnet(predictors.train, response.train, family="binomial", lambda=lambda)
  
  # Save model members and parameters
  sig.params <- logis.model$beta[which(logis.model$beta!=0)]
  sig.modules <- rownames(logis.model$beta)[which(logis.model$beta!=0)]
  LASSO.models[[count]] <- data.frame(Module=sig.modules, Parameter=sig.params, Num.Modules=rep(length(sig.modules), length(sig.modules)))
  
  # Obtain predictions for validation set
  logis.pred <- predict(logis.model, newx=predictors.valid, type="class")
  
  # Compute metrics
  predicted.val <- mgsub(logis.pred, c("AD", "Control"), c(1, 0)) %>%
    as.numeric()
  true.val <- mgsub(response.valid, c("AD", "Control"), c(1,0)) %>%
    as.numeric()
  pred <- prediction(predicted.val, true.val)
  auc <- performance(pred, measure="auc")
  
  # Compute and save metrics
  auc <- auc@y.values[[1]]
  metrics[count,] <- c(
    as.numeric(auc),
    sum(response.valid == logis.pred) / length(response.valid),
    logis.model$dev.ratio,
    lambda,
    sd(min.lambda$lambda)
  )
}

# Create data frame with important metrics
metrics$Model.Index <- rownames(metrics)
LASSO.metrics <- data.frame(
  AUC=mean(metrics$AUC), 
  AUC.SD=sd(metrics$AUC), 
  ACC=mean(metrics$ACC), 
  ACC.SD=sd(metrics$ACC),
  LRI=mean(metrics$LRI), 
  LRI.SD=sd(metrics$LRI), 
  Lambda.SD=sd(metrics$Lambda)
)

# Store LASSO Metrics and Models
write.table(LASSO.metrics, file=paste0("results/4_LASSO/", COHORT, "_LASSO_metrics.txt"), sep="\t", row.names=FALSE)
saveRDS(metrics, file=paste0("results/4_LASSO/", COHORT, "_LASSO_metrics.RDS"))
saveRDS(LASSO.models, file=paste0("results/4_LASSO/", COHORT, "_LASSO_models.RDS"))

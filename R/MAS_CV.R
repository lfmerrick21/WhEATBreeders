MAS_CV <- function(phenotype,PCA=NULL,CV=NULL, folds = 5)
{

  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  library(dplyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    pheno_train <- phenotype[,2]
    pheno_train[-fold_indices] <- NA

    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]


    if(length(PCA)!=0){
      myPCA_train <- PCA
      myPCA_test <- PCA
    }
    gc()

    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)

      GLM_data <- data.frame(pheno_train, MAS_train)

      names(GLM_data)[1] <- "Y"

      MAS_model <- lm(Y ~ ., data = GLM_data)
      predictions <- predict(MAS_model,MAS_test)

      acc <- cor(predictions[-fold_indices], myY_test, use = "pairwise.complete")
      sacc <- cor(predictions[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions[-fold_indices],obs=myY_test)
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])



      if(!is.null(PCA)){
        MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
        MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))

        GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

        names(GLM_data_PC)[1] <- "Y"
        #Linear model to calculate effects
        #You can run all signficant markers at once to see cumulative effect
        MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
        predictions_PC <- predict(MAS_model_PC, MAS_test_PC)

      acc_PC <- cor(predictions_PC[-fold_indices], myY_test, use = "pairwise.complete")
      sacc_PC <- cor(predictions_PC[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics_PC=postResample(pred=predictions_PC[-fold_indices],obs=myY_test)
      results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
      prediction_PC=data.frame(Fold=rep(i,length(phenotype[-fold_indices,])),phenotype[-fold_indices,],GEBV_PC=predictions_PC[-fold_indices])}

    if(!is.null(PCA)){
    Predictions<-cbind(prediction,prediction_PC[,4])
    BGLR_acc_results[[i]] <- list(results,results_PC)}else{
      Predictions<-prediction
      BGLR_acc_results[[i]] <- list(results)
    }
    Predictions_ALL=rbind(Predictions_ALL,Predictions)

    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  if(!is.null(PCA)){
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  }else{
    model_vect <- c("Pearson","Spearman","RMSE","R2","MAE")
    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
    for (i in 1:length(BGLR_acc_results))
    {
      results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
      BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
    }
    names(BGLR_acc_table) <- c("fold", "model", "r")
    data_wide <- spread(BGLR_acc_table, model, r)
    acc_fold=data.frame(data_wide)
    results=colMeans(acc_fold[2:6], na.rm=TRUE)
  }

  results_ALL=list(Results=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

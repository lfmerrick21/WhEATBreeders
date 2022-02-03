rrBLUP_CV <- function(genotypes, phenotype,PCA=NULL,CV=NULL, folds = 5)
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

    # Split into training and testing data
    pheno_train <- phenotype[fold_indices,]
    pheno_test=phenotype[-fold_indices,]

    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    if(length(PCA)!=0){
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    }
    gc()

    if(length(CV)==0){rrBLUP_model <- mixed.solve(y = myY_train,
                                                  Z = myGD_train)
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)

    if(!is.null(PCA)){

    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = myPCA_train)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Fold=rep(i,length(phenotype[-fold_indices,])),phenotype[-fold_indices,],GEBV_PC=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    }
    }else{

      #CV=impute(as.factor(CV))
      if(length(ncol(CV))==0){
        myCV_train <- CV[fold_indices]
        myCV_test <- CV[-fold_indices]
      }else{
        myCV_train <- CV[fold_indices,]
        myCV_test <- CV[-fold_indices,]
      }

      fix_train <- as.matrix(myCV_train)
      fix_test  <- as.matrix(myCV_test)

      p <- ncol(fix_train)
      XtX <- crossprod(fix_train, fix_train)
      rank.X <- qr(XtX)$rank
      if (rank.X < p) {
        sm2=findLinearCombos(fix_train)$remove
        fix_train=fix_train[,-sm2]
        fix_test=fix_test[,-sm2]}

      gc()
      rrBLUP_model <- mixed.solve(y = myY_train,
                                  Z = myGD_train,
                                  X = fix_train)

      pred_effects <- myGD_test %*% rrBLUP_model$u
      fix_effects <- fix_test  %*% rrBLUP_model$beta
      predictions <- c(pred_effects) + c(fix_effects)

      acc <- cor(predictions, myY_test, use = "pairwise.complete")
      sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions,obs=myY_test)
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)



      if(!is.null(PCA)){
      fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
      fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))

      p <- ncol(fix_train_PC)
      XtX <- crossprod(fix_train_PC, fix_train_PC)
      rank.X <- qr(XtX)$rank
      if (rank.X < p) {
        sm2=findLinearCombos(fix_train_PC)$remove
        fix_train_PC=fix_train_PC[,-sm2]
        fix_test_PC=fix_test_PC[,-sm2]}

      rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                     Z = myGD_train,
                                     X = fix_train_PC)
      pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
      fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
      predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
      acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
      sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
      results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
      prediction_PC=data.frame(Fold=rep(i,length(phenotype[-fold_indices,])),phenotype[-fold_indices,],GEBV_PC=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)}

    }
    if(!is.null(PCA)){
    Predictions<-cbind(prediction,prediction_PC[,4:6])
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

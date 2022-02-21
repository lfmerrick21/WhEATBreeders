MAS_CV <- function(phenotype,PCA=NULL,CV=NULL, folds = 5,transformation=NULL)
{

  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    if(transformation=="sqrt"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
      phenotype[-fold_indices,2] <-sqrt(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-sqrt(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="log"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
      phenotype[-fold_indices,2] <-log(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-log(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="boxcox"){
      phenotype[-fold_indices,2] <-boxcox_t(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-boxcox_t(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="none"){
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }


    if(length(PCA)!=0){
      myPCA_train <- PCA
      myPCA_test <- PCA
    }
    gc()

    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)


    if(!is.null(PCA)){
      MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
      MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))

      GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

      names(GLM_data_PC)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
      predictions <- predict(MAS_model_PC, MAS_test_PC)

      acc <- cor(predictions[-fold_indices], myY_test, use = "pairwise.complete")
      sacc <- cor(predictions[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions[-fold_indices],obs=myY_test)
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])

    }else{
      GLM_data <- data.frame(pheno_train, MAS_train)

      names(GLM_data)[1] <- "Y"

      MAS_model <- lm(Y ~ ., data = GLM_data)
      predictions <- predict(MAS_model,MAS_test)

      acc <- cor(predictions[-fold_indices], myY_test, use = "pairwise.complete")
      sacc <- cor(predictions[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions[-fold_indices],obs=myY_test)
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])
    }


    Predictions<-prediction
    BGLR_acc_results[[i]] <- list(results)
    Predictions_ALL=rbind(Predictions_ALL,Predictions)

    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }

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

  results_ALL=list(Results=results,Predictions=Predictions_ALL)
  return(results_ALL)
}

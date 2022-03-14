
MAS_VS <- function(train_phenotype,train_PCA=NULL,train_CV=NULL,test_phenotype,test_PCA=NULL,test_CV=NULL, transformation=NULL)
{


  # Split into training and testing data

  if(transformation=="sqrt"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
    test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)
    train_phenotype[,2] <-sqrt(train_phenotype[,2])
    test_phenotype[,2] <-sqrt(test_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="log"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
    test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
    train_phenotype[,2] <-log(train_phenotype[,2])
    test_phenotype[,2] <-log(test_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="boxcox"){
    train_phenotype[,2] <-boxcox_t(train_phenotype[,2])
    test_phenotype[,2] <-boxcox_t(test_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="none"){
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

  }
  #Make sure columns have same column names

  if(length(ncol(train_CV))==0){
    CV<-c(train_CV,test_CV)
  }else{
    CV<-rbind(train_CV,test_CV)
  }

    if(!is.null(train_PC)){
      PCA<-rbind(train_PC,test_PC)
    }else{
      PCA=NULL
    }




    if(!is.null(train_PCA)){
      MAS_train_PC <- data.frame(cbind(CV,PCA))
      MAS_test_PC  <- data.frame(cbind(CV,PCA))

      GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

      names(GLM_data_PC)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
      predictions <- predict(MAS_model_PC, MAS_test_PC)

      acc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete")
      sacc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions[-c(1:length(train_phenotype[,2]))],obs=test_phenotype[,2])
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(test_phenotype,GEBV=predictions[-c(1:length(train_phenotype[,2]))])

    }else{
      MAS_train <- data.frame(CV)
      MAS_test  <- data.frame(CV)
      GLM_data <- data.frame(pheno_train, MAS_train)

      names(GLM_data)[1] <- "Y"

      MAS_model <- lm(Y ~ ., data = GLM_data)
      predictions <- predict(MAS_model,MAS_test)

      acc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete")
      sacc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions[-c(1:length(train_phenotype[,2]))],obs=test_phenotype[,2])
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(test_phenotype,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
    }




  names(results)<- c("Pearson","Spearman","RMSE","R2","MAE")

  results_ALL=list(Results=results,Predictions=prediction)
  return(results_ALL)
}

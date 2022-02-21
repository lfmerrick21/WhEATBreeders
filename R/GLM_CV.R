GLM_CV <- function(genotypes, phenotype,markers=NULL,Kernel="Markers",fam="poisson", folds = 5,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL){

  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]

    if(Kernel=="Markers"){
      if(!is.null(markers)){
        samp=sample(1:ncol(genotypes), markers)
        m_samp=genotypes[,samp]
        myGD_train <- m_samp[fold_indices,]
        myGD_test <- m_samp[-fold_indices,]
      }else{
        myGD_train <- genotypes[fold_indices,]
        myGD_test <- genotypes[-fold_indices,]
      }

    }else{
      maf <- calc_maf_apply(genotypes, encoding = c(0, 1, 2))
      mono_indices <- which(maf ==0)
      if(length(mono_indices)!=0){
        genotypes = genotypes[,-mono_indices]
      }

      if(Sparse==TRUE){
        X=as.matrix(genotypes)
        X=apply(genotypes,2,as.numeric)
        #sum(rowSums(is.na(X)))
        XS<-scale(X,center=TRUE,scale=TRUE)
        #sum(rowSums(is.na(XS)))
        library(caret)
        nzv <- nearZeroVar(XS)
        XSZV <- XS[, -nzv]
        #sum(rowSums(is.na(XSZV)))
        X=XSZV
        K=Sparse_Kernel_Construction(m=m,X=X,name=Kernel, degree=degree,nL=nL)
      }else{
        X=as.matrix(genotypes)
        X=apply(genotypes,2,as.numeric)
        #sum(rowSums(is.na(X)))
        XS<-scale(X,center=TRUE,scale=TRUE)
        #sum(rowSums(is.na(XS)))
        library(caret)
        nzv <- nearZeroVar(XS)
        XSZV <- XS[, -nzv]
        #sum(rowSums(is.na(XSZV)))
        X=XSZV
        K=Kernel_computation(X=X,name=Kernel,degree=degree, nL=nL)
      }
      myGD_train <- K[fold_indices,]
      myGD_test <- K[-fold_indices,]
    }


    myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
    myGD_test1=as.matrix(sapply(myGD_test, as.numeric))

    #glmnet.control(mxitnr = 50)
    if(fam=="poisson"){
      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = poisson(),
                                   alpha=1,type.measure="mse", standardize = FALSE,
                                   intercept = FALSE)
    }
    if(fam=="quasipoisson"){
      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = quasipoisson(),
                                   alpha=1,type.measure="mse", standardize = FALSE,
                                   intercept = FALSE)
    }
    if(fam=="negative.binomial"){
      A1_RR=glmnetUtils::cv.glmnet(myGD_train2, myY_train, family = negative.binomial(theta = 5),
                                   alpha=1,type.measure="mse", standardize = FALSE,
                                   intercept = FALSE)
    }
    predictions= as.numeric(predict(A1_RR,newx=myGD_test1,s='lambda.min',type='response'))

    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions)

    Predictions<-prediction
    BGLR_acc_results[[i]] <- list(results)
    Predictions_ALL=rbind(Predictions_ALL,Predictions)

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

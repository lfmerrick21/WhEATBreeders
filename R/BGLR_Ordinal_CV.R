BGLR_Ordinal_CV <- function(genotypes, phenotype,model="BL",Kernel="Markers",markers=NULL,PCA=NULL,CV=NULL, nIter = 15000, burnIn = 5000, folds = 5,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL)
{
  # Make the CV list


  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_metrics<- list()
  Predictions_ALL<-list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])
    if(Kernel=="Markers"){
      if(!is.null(markers)){
        samp=sample(1:ncol(genotypes), markers)
        genotypes=genotypes[,samp]
      }else{
        genotypes <- genotypes
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

    }
    # Split into training and testing data
    pheno_train <- phenotype[,2]
    pheno_train[pheno_train=="NaN"]<-NA
    pheno_train=droplevels(pheno_train)
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using BGLR
    ##Ordinal
    #Without PCs
    if(length(CV)==0){
      if(!is.null(PCA)){
        if(Kernel=="Markers"){
          BO_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
        }else{
          BO_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)

        }

      }else{
        if(Kernel=="Markers"){
          BO_ETA<-list(list(X=as.matrix(genotypes),model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
        }else{
          BO_ETA<-list(list(K=K,model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)

        }
      }
    }else{
      if(!is.null(PCA)){
        fix_PC=as.matrix(cbind(CV,PCA))
        if(Kernel=="Markers"){
          BO_ETA<-list(list(X=fix_PC,model="FIXED"),list(X=as.matrix(genotypes),model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
        }else{
          BO_ETA<-list(list(X=fix_PC,model="FIXED"),list(K=K,model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)

        }

      }else{
        if(Kernel=="Markers"){
          BO_ETA<-list(list(X=CV,model="FIXED"),list(X=as.matrix(genotypes),model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
        }else{
          BO_ETA<-list(list(X=CV,model="FIXED"),list(K=K,model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)

        }
      }

    }
    BO_acc <- cor(as.numeric(phenotype[-fold_indices,2]), BO_predictions[-fold_indices],use = "complete.obs")
    BO_sacc <- cor(as.numeric(phenotype[-fold_indices,2]), BO_predictions[-fold_indices],use = "complete.obs", method = c("spearman"))
    DF=BO_results$probs[-fold_indices,]
    probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    #probs=colnames(DF)[max.col(DF,ties.method="first")]
    BO_acc_cat=cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(probs),use = "complete.obs")
    metrics=postResample(as.factor(probs),phenotype[-fold_indices,2])
    tests=data.frame(Observed=phenotype[-fold_indices,2],BO_results$probs[-fold_indices,],Predicted=factor(probs))
    mets=confusionMatrix(data = tests$pred, reference = tests$obs)
    results=c(ACC=BO_acc,SACC=BO_sacc,C_ACC=BO_acc_cat,metrics)
    #ALL
    Predictions<-data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,1],tests)
    Predictions_ALL[[i]]<-Predictions
    BGLR_acc_results[[i]] <- results
    BGLR_acc_metrics[[i]] <- mets

  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","Categorical","R2","Kappa")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:6], na.rm = TRUE)
  results_ALL=list(Results=results,CM=BGLR_acc_metrics,Predictions=Predictions_ALL)
  return(results_ALL)
}

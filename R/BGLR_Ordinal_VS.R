
BGLR_Ordinal_VS <- function(train_genotypes, train_phenotype,train_PCA=NULL,train_CV=NULL,test_genotypes,test_phenotype,test_PCA=NULL,test_CV=NULL,model="BL",Kernel="Markers",markers=NULL, nIter = 15000, burnIn = 5000,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL)
{

  # Make the CV list
    if(Kernel=="Markers"){
      if(!is.null(markers)){
        genotypes=rbind(train=train_genotypes,test=test_genotypes)
        samp=sample(1:ncol(genotypes), markers)
        genotypes=genotypes[,samp]
      }else{
        genotypes=rbind(train=train_genotypes,test=test_genotypes)
      }


    }else{
      genotypes=rbind(train=train_genotypes,test=test_genotypes)
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
  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[pheno_train=="NaN"]<-NA
  pheno_train=droplevels(pheno_train)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA

    # Calculate the GS model using BGLR
    ##Ordinal
    #Without PCs
    if(length(CV)==0){
      if(!is.null(PCA)){
        PCA<-rbind(train=train_PCA,test=test_PCA)
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
      CV<-rbind(train=train_CV,test=test_CV)
      if(!is.null(PCA)){
        PCA<-rbind(train=train_PCA,test=test_PCA)
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
    BO_acc <- cor(as.numeric(test_phenotype[,2]), BO_predictions[-c(1:length(train_phenotype[,2]))],use = "complete.obs")
    BO_sacc <- cor(as.numeric(test_phenotype[,2]), BO_predictions[-c(1:length(train_phenotype[,2]))],use = "complete.obs", method = c("spearman"))
    DF=BO_results$probs[-c(1:length(train_phenotype[,2])),]
    probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    #probs=colnames(DF)[max.col(DF,ties.method="first")]
    BO_acc_cat=cor(as.numeric(test_phenotype[,2]), as.numeric(probs),use = "complete.obs")
    metrics=postResample(as.factor(probs),test_phenotype[,2])
    tests=data.frame(Observed=test_phenotype[,2],BO_results$probs[-c(1:length(train_phenotype[,2])),],Predicted=factor(probs))
    mets=confusionMatrix(data = tests$pred, reference = tests$obs)
    results=c(ACC=BO_acc,SACC=BO_sacc,C_ACC=BO_acc_cat,metrics)
    #ALL
    Predictions<-data.frame(test_phenotype[,1],tests)

    names(results)<- c("Pearson","Spearman","Categorical","R2","Kappa")

    results_ALL=list(Results=results,CM=mets,Predictions=Predictions)

  return(results_ALL)
}

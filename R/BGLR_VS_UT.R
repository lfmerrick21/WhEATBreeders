
BLGR_VS_UT <- function(train_genotypes, train_phenotype,train_PCA=NULL,train_CV=NULL,test_genotypes,test_phenotype,test_PCA=NULL,test_CV=NULL,model="RKHS",Kernel="Markers",markers=NULL, nIter = 15000, burnIn = 5000,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,transformation=NULL)
{

  # Make the CV list

  if(Kernel=="Markers"){
    if(!is.null(markers)){
      genotypes=rbind(train=train_genotypes,test=test_genotypes)
      samp=sample(1:ncol(genotypes), markers)
      genotypes=genotypes[,samp]
    }else{
      genotypes <- genotypes
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
  if(transformation=="sqrt"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
    train_phenotype[,2] <-sqrt(train_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="log"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
    train_phenotype[,2] <-log(train_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="boxcox"){
    train_phenotype[,2] <-boxcox_t(train_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="none"){
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

  }


  if(length(CV)==0){

    if(!is.null(train_PCA)){
      PCA<-rbind(train=train_PCA,test=test_PCA)
      if(Kernel=="Markers"){
        ETA<-list(list(X=PCA,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)}
      else{
        ETA<-list(list(X=PCA,model="FIXED"),G=list(K=K,model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)
      }
      prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
    }else{
      if(Kernel=="Markers"){
        ETA<-list(G=list(X=as.matrix(genotypes),model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)
      }else{

        ETA<-list(G=list(K=K,model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)

      }
      prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])

    }
  }else{

    CV<-rbind(train=train_CV,test=test_CV)

    if(!is.null(PCA)){
      PCA<-rbind(train=train_PCA,test=test_PCA)

      fix_PC=as.matrix(cbind(CV,PCA))


      if(Kernel=="Markers"){
        ETA<-list(list(X=fix_PC,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)
      }else{
        ETA<-list(list(X=fix_PC,model="FIXED"),G=list(K=K,model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)
      }
      prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
    }else{
      gc()
      if(Kernel=="Markers"){
        ETA<-list(list(X=CV,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)
      }else{
        ETA<-list(list(X=CV,model="FIXED"),G=list(K=K,model=model))
        BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
        predictions <- predict(BGLR_results)

      }
      prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])

    }

  }

  results_ALL=list(Predictions=prediction)
  return(results_ALL)
}

GLM_VS_UT <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,markers=NULL,Kernel="Markers",fam="poisson",Sparse=FALSE,m=NULL,degree=NULL, nL=NULL){

  # Split into training and testing data

  myY_train <- train_phenotype[,2]
  test_phenotype[,2]<-NA

  if(Kernel=="Markers"){
    if(!is.null(markers)){
      samp=sample(1:ncol(train_genotypes), markers)
      myGD_train <- train_genotypes[,samp]
      myGD_test <- test_genotypes[,samp]
      genotypes=rbind(train=myGD_train,test=myGD_test)
    }else{
      myGD_train <- train_genotypes
      myGD_test <- test_genotypes
      genotypes=rbind(train=myGD_train,test=myGD_test)
    }

  }else{
    genotypes=rbind(train=myGD_train,test=myGD_test)
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
    myGD_train <- K[c(1:length(train_phenotype[,2])),]
    myGD_test <- K[-c(1:length(train_phenotype[,2])),]
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

  prediction=data.frame(test_phenotype,GEBV=predictions)

  results_ALL=list(Predictions=prediction)
  return(results_ALL)
}

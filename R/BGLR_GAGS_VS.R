
BGLR_GAGS_CV <- function(train_genotypes, train_phenotype,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,model="RKHS",Kernel="Markers",GWAS="BLINK",alpha=0.05,threshold=NULL, markers=NULL, nIter = 15000, burnIn = 5000,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,QTN=10,PCA.total=3,transformation=NULL)
{

  # Make the CV list

    if(Kernel=="Markers"){
      if(!is.null(markers)){
        genotypes=rbind(train=myGD_train,test=myGD_test)
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

    GWASR<- GAPIT(Y = train_phenotype,
                  GD = train_GD,
                  GM = train_GM,
                  PCA.total=PCA.total,
                  model = GWAS,
                  file.output=F)

    if(threshold=="Bonferonni"){

      GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
      if(length(GWASSM)==0){

        acc <- NA
        sacc <- NA
        metrics=c(RMSE=NA,Rsquared=NA,MAE=NA)
        results=c(ACC=acc,SACC=sacc,metrics)

        prediction=data.frame(test_phenotype,GEBV=NA)

      }else{
        sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

        myCV<-as.matrix(genotypes[,sm])

        if(!is.null(PCA)){
          myPCA_train <- train_PCA
          myPCA_test <- test_PCA
          PCA<-rbind(train_PCA,test_PCA)
          fix_PC=as.matrix(cbind(myCV,PCA))

          if(Kernel=="Markers"){
            ETA<-list(list(X=fix_PC,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
            BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
            predictions <- predict(BGLR_results)
          }else{
            ETA<-list(list(X=fix_PC,model="FIXED"),G=list(K=K,model=model))
            BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
            predictions <- predict(BGLR_results)
          }
          acc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete")
          sacc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete", method = c("spearman"))
          metrics=postResample(pred=predictions[-c(1:length(train_phenotype[,2]))],obs=test_phenotype[,2])
          results=c(ACC=acc,SACC=sacc,metrics)
          prediction=data.frame(test_phenotype,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
        }else{
          if(Kernel=="Markers"){
            ETA<-list(list(X=myCV,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
            BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
            predictions <- predict(BGLR_results)
          }else{
            ETA<-list(list(X=myCV,model="FIXED"),G=list(K=K,model=model))
            BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
            predictions <- predict(BGLR_results)
          }
          acc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete")
          sacc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete", method = c("spearman"))
          metrics=postResample(pred=predictions[-c(1:length(train_phenotype[,2]))],obs=test_phenotype[,2])
          results=c(ACC=acc,SACC=sacc,metrics)
          prediction=data.frame(test_phenotype,GEBV=predictions[-c(1:length(train_phenotype[,2]))])

        }
      }



    }else{
      top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
      GWASSM=top10[1:QTN,]$SNP
      myCV<-as.matrix(genotypes[,GWASSM])

      if(!is.null(PCA)){
        myPCA_train <- train_PCA
        myPCA_test <- test_PCA
        PCA<-rbind(train_PCA,test_PCA)
        fix_PC=as.matrix(cbind(myCV,PCA))

        if(Kernel=="Markers"){
          ETA<-list(list(X=fix_PC,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
          BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
          predictions <- predict(BGLR_results)
        }else{
          ETA<-list(list(X=fix_PC,model="FIXED"),G=list(K=K,model=model))
          BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
          predictions <- predict(BGLR_results)
        }
        acc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete")
        sacc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions[-c(1:length(train_phenotype[,2]))],obs=test_phenotype[,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(test_phenotype,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
      }else{
        if(Kernel=="Markers"){
          ETA<-list(list(X=myCV,model="FIXED"),G=list(X=as.matrix(genotypes),model=model))
          BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
          predictions <- predict(BGLR_results)
        }else{
          ETA<-list(list(X=myCV,model="FIXED"),G=list(K=K,model=model))
          BGLR_results <- BGLR(y = pheno_train, ETA = ETA, nIter=nIter, burnIn=burnIn)
          predictions <- predict(BGLR_results)
        }
        acc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete")
        sacc <- cor(predictions[-c(1:length(train_phenotype[,2]))], test_phenotype[,2], use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions[-c(1:length(train_phenotype[,2]))],obs=test_phenotype[,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(test_phenotype,GEBV=predictions[-c(1:length(train_phenotype[,2]))])

      }

    }
    names(results)<- c("Pearson","Spearman","RMSE","R2","MAE")

    results_ALL=list(Results=results,Predictions=prediction)
  return(results_ALL)
}

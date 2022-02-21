

BGLR_Ordinal_GAGS_VS <- function(train_genotypes, train_phenotype,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,model="BL",Kernel="Markers",GWAS="BLINK",alpha=0.05,threshold=NULL, markers=NULL,nIter = 15000, burnIn = 5000,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,QTN=10,PCA.total=3)
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

  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[pheno_train=="NaN"]<-NA
  pheno_train=droplevels(pheno_train)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA

    GWASR<- GAPIT(Y = train_phenotype,
                  GD = train_GD,
                  GM = train_GM,
                  PCA.total=PCA.total,
                  model = GWAS,
                  file.output=F)

    if(threshold=="Bonferonni"){

      GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
      if(length(GWASSM)==0){

        BO_acc <- NA
        BO_sacc <- NA
        BO_acc_cat<-NA
        metrics=c(Rsquared=NA,Kappa=NA)
        results=c(ACC=BO_acc,SACC=BO_sacc,C_ACC=BO_acc_cat,metrics)
        mets<-NA
        #ALL
        prediction<-data.frame(test_phenotype[,1],Observed=NA,Predicted=NA)

      }else{
        sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

        myCV<-as.matrix(genotypes[,sm])

        if(!is.null(PCA)){
          PCA<-rbind(train=train_PCA,test=test_PCA)
          fix_PC=as.matrix(cbind(myCV,PCA))

          if(Kernel=="Markers"){
            BO_ETA<-list(list(X=fix_PC,model="FIXED"),list(X=as.matrix(genotypes),model=model))
            BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
            BO_predictions <- predict(BO_results)
          }else{
            BO_ETA<-list(list(X=fix_PC,model="FIXED"),list(K=K,model=model))
            BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
            BO_predictions <- predict(BO_results)
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
        }else{
          if(Kernel=="Markers"){
            BO_ETA<-list(list(X=myCV,model="FIXED"),list(X=as.matrix(genotypes),model=model))
            BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
            BO_predictions <- predict(BO_results)
          }else{
            BO_ETA<-list(list(X=myCV,model="FIXED"),list(K=K,model=model))
            BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
            BO_predictions <- predict(BO_results)
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

        }
        Predictions<-data.frame(test_phenotype[,1],tests)
      }



    }else{
      top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
      GWASSM=top10[1:QTN,]$SNP
      myCV<-as.matrix(genotypes[,GWASSM])

      if(!is.null(PCA)){
        PCA<-rbind(train=train_PCA,test=test_PCA)
        fix_PC=as.matrix(cbind(myCV,PCA))

        if(Kernel=="Markers"){
          BO_ETA<-list(list(X=fix_PC,model="FIXED"),list(X=as.matrix(genotypes),model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
        }else{
          BO_ETA<-list(list(X=fix_PC,model="FIXED"),list(K=K,model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA_PC, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
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
      }else{
        if(Kernel=="Markers"){
          BO_ETA<-list(list(X=myCV,model="FIXED"),list(X=as.matrix(genotypes),model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
        }else{
          BO_ETA<-list(list(X=myCV,model="FIXED"),list(K=K,model=model))
          BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
          BO_predictions <- predict(BO_results)
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

      }
      Predictions<-data.frame(test_phenotype[,1],tests)

    }
  names(results)<- c("Pearson","Spearman","Categorical","R2","Kappa")

  results_ALL=list(Results=results,CM=mets,Predictions=Predictions)
  return(results_ALL)
}

BGLR_Ordinal_GAGS_CV <- function(genotypes, phenotype,model="BL",Kernel="Markers",Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05,threshold=NULL, markers=NULL, folds = 5,nIter = 15000, burnIn = 5000,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,QTN=10,PCA.total=3)
{
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_metrics<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
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

    if(!is.null(markers)){
      samp=sample(2:ncol(GD), markers)
      m_samp=GD[,samp]
      myGD_train <- m_samp[fold_indices,]
      myGD_test <- m_samp[-fold_indices,]
    }else{
      myGD_train <- GD[fold_indices,]
      myGD_test <- GD[-fold_indices,]
    }

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    pheno_train <- phenotype[,2]
    pheno_train[-fold_indices] <- NA

    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
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
        prediction<-data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,1],Observed=phenotype[-fold_indices,2],Predicted=NA)
        Predictions<-prediction
        BGLR_acc_results[[i]] <- results
        BGLR_acc_metrics[[i]] <- mets
        Predictions_ALL[[i]]<-Predictions

      }else{
        sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

        myCV<-as.matrix(genotypes[,sm])

        if(!is.null(PCA)){
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

        }
        Predictions<-data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,1],tests)
        Predictions_ALL[[i]]<-Predictions
        BGLR_acc_results[[i]] <- results
        BGLR_acc_metrics[[i]] <- mets
      }



    }else{
      top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
      GWASSM=top10[1:QTN,]$SNP
      myCV<-as.matrix(genotypes[,GWASSM])

      if(!is.null(PCA)){
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

      }
      Predictions<-data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,1],tests)
      Predictions_ALL[[i]]<-Predictions
      BGLR_acc_results[[i]] <- results
      BGLR_acc_metrics[[i]] <- mets

    }


  }
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
  results=colMeans(acc_fold[2:6], na.rm=TRUE)

  results_ALL=list(Results=results,CM=BGLR_acc_metrics,Predictions=Predictions_ALL)
  return(results_ALL)
}

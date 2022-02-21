BGLR_GAGS_CV <- function(genotypes, phenotype,model="RKHS",Kernel="Markers",Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05,threshold=NULL, markers=NULL, folds = 5,nIter = 15000, burnIn = 5000,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,QTN=10,PCA.total=3,transformation=NULL)
{

  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
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

    if(transformation=="sqrt"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
      phenotype[-fold_indices,2] <-sqrt(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-sqrt(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      Y[,2]=replace(Y[,2], Y[,2] < 0, 0)
      Y_train <-Y[fold_indices,c(1,2)]
      Y_test <-Y[-fold_indices,c(1,2)]
      Y_train[,2] <-sqrt(Y_train[,2])
      Y_test[,2] <-sqrt(Y_train[,2])
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="log"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
      phenotype[-fold_indices,2] <-log(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-log(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      Y[,2]=replace(Y[,2], Y[,2] <= 0, 0.000001)
      Y_train <-Y[fold_indices,c(1,2)]
      Y_test <-Y[-fold_indices,c(1,2)]
      Y_train[,2] <-log(Y_train[,2])
      Y_test[,2] <-log(Y_train[,2])
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA

    }

    if(transformation=="boxcox"){
      phenotype[-fold_indices,2] <-boxcox_t(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-boxcox_t(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      Y_train <-Y[fold_indices,c(1,2)]
      Y_test <-Y[-fold_indices,c(1,2)]
      Y_train[,2] <-boxcox_t(Y_train[,2])
      Y_test[,2] <-boxcox_t(Y_train[,2])
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA

    }

    if(transformation=="none"){
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      Y_train <-Y[fold_indices,c(1,2)]
      Y_test <-Y[-fold_indices,c(1,2)]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
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

        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=NA)

        Predictions<-prediction
        BGLR_acc_results[[i]] <- list(results)
        Predictions_ALL=rbind(Predictions_ALL,Predictions)

      }else{
        sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

        myCV<-as.matrix(genotypes[,sm])

        if(!is.null(PCA)){
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
          acc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete")
          sacc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete", method = c("spearman"))
          metrics=postResample(pred=predictions[-fold_indices],obs=phenotype[-fold_indices,2])
          results=c(ACC=acc,SACC=sacc,metrics)
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])
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
          acc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete")
          sacc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete", method = c("spearman"))
          metrics=postResample(pred=predictions[-fold_indices],obs=phenotype[-fold_indices,2])
          results=c(ACC=acc,SACC=sacc,metrics)
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])

        }
        Predictions<-prediction
        BGLR_acc_results[[i]] <- list(results)
        Predictions_ALL=rbind(Predictions_ALL,Predictions)
      }



    }else{
      top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
      GWASSM=top10[1:QTN,]$SNP
      myCV<-as.matrix(genotypes[,GWASSM])

      if(!is.null(PCA)){
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
        acc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete")
        sacc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions[-fold_indices],obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])
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
        acc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete")
        sacc <- cor(predictions[-fold_indices], phenotype[-fold_indices,2], use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions[-fold_indices],obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions[-fold_indices])

      }
      Predictions<-prediction
      BGLR_acc_results[[i]] <- list(results)
      Predictions_ALL=rbind(Predictions_ALL,Predictions)

    }


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

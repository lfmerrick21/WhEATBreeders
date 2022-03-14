

rrBLUP_VS <- function(train_genotypes, train_phenotype,train_PCA=NULL,train_CV=NULL,test_genotypes, test_phenotype,test_PCA=NULL,test_CV=NULL,Kernel="Markers",markers=NULL, Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,transformation=NULL)
{

  # Make the CV list

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
      # Calculate the GS model using rrBLUP
      myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
      myGD_train=apply(myGD_train,2,as.numeric)
      myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
      myGD_test=apply(myGD_test,2,as.numeric)
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

    }
    # Split into training and testing data
  if(transformation=="sqrt"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
    test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)
    train_phenotype[,2] <-sqrt(train_phenotype[,2])
    test_phenotype[,2] <-sqrt(test_phenotype[,2])
    myY_train <- train_phenotype[,2]
    myY_test <- test_phenotype[,2]
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="log"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
    test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
    train_phenotype[,2] <-log(train_phenotype[,2])
    test_phenotype[,2] <-log(test_phenotype[,2])
    myY_train <- train_phenotype[,2]
    myY_test <- test_phenotype[,2]
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

  }

  if(transformation=="boxcox"){
    train_phenotype[,2] <-boxcox_t(train_phenotype[,2])
    test_phenotype[,2] <-boxcox_t(test_phenotype[,2])
    myY_train <- train_phenotype[,2]
    myY_test <- test_phenotype[,2]
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="none"){
    myY_train <- train_phenotype[,2]
    myY_test <- test_phenotype[,2]
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }



    if(length(train_PCA)!=0){
      myPCA_train <- train_PCA
      myPCA_test <- test_PCA
      PCA<-rbind(train=train_PCA,test=test_PCA)
    }

    if(length(train_CV)==0){

      if(!is.null(PCA)){
        if(Kernel=="Markers"){
          rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                         Z = myGD_train,
                                         X = myPCA_train)
          pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
          fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
          predictions <- c(pred_effects_PC) + c(fix_effects_PC)
        }else{
          gBLUP_model <- mixed.solve(y = pheno_train,
                                     K = K,
                                     X=PCA)
          pred_effects <- gBLUP_model$u[-c(1:length(train_phenotype[,2]))]
          fix_effects <- as.matrix(PCA[-c(1:length(train_phenotype[,2])),])  %*% gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)
        }
        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
      }else{
        if(Kernel=="Markers"){
          rrBLUP_model <- mixed.solve(y = myY_train,
                                      Z = myGD_train)
          pred_effects <- myGD_test %*% rrBLUP_model$u
          fix_effects <- rrBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)}
        else{
          gBLUP_model <- mixed.solve(y = pheno_train,
                                     K = K)
          pred_effects <- gBLUP_model$u[-c(1:length(train_phenotype[,2]))]
          fix_effects <- gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)

        }
        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)

      }
    }else{

        myCV_train <- train_CV
        myCV_test <- test_CV
        CV<-rbind(train=train_CV,test=test_CV)

      if(!is.null(train_PCA)){
        if(Kernel=="Markers"){
          fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
          fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))

          #p <- ncol(fix_train_PC)
          #XtX <- crossprod(fix_train_PC, fix_train_PC)
          #rank.X <- qr(XtX)$rank
          #if (rank.X < p) {
          #sm2=findLinearCombos(fix_train_PC)$remove
          #fix_train_PC=fix_train_PC[,-sm2]
          #fix_test_PC=fix_test_PC[,-sm2]}

          fix_train_PC=make_full_rank(fix_train_PC)
          fix_test_PC=fix_test_PC[,colnames(fix_train_PC)]

          rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                         Z = myGD_train,
                                         X = fix_train_PC)
          pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
          fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
          predictions <- c(pred_effects_PC) + c(fix_effects_PC)}
        else{
          fix_PC=as.matrix(cbind(CV,PCA))
          fix_PC=make_full_rank(fix_PC)
          gBLUP_model <- mixed.solve(y = pheno_train,
                                     K = K,
                                     X=fix_PC)
          pred_effects_PC <- gBLUP_model$u[-c(1:length(train_phenotype[,2]))]
          fix_effects_PC <- as.matrix(fix_PC[-c(1:length(train_phenotype[,2])),])  %*% gBLUP_model$beta
          predictions <- c(pred_effects_PC) + c(fix_effects_PC)
        }
        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
      }else{
        if(Kernel=="Markers"){
          fix_train <- as.matrix(myCV_train)
          fix_test  <- as.matrix(myCV_test)

          #p <- ncol(fix_train)
          #XtX <- crossprod(fix_train, fix_train)
          #rank.X <- qr(XtX)$rank
          #if (rank.X < p) {
          #sm2=findLinearCombos(fix_train)$remove
          #fix_train=fix_train[,-sm2]
          #fix_test=fix_test[,-sm2]}

          fix_train=make_full_rank(fix_train)
          fix_test=fix_test_PC[,colnames(fix_train)]
          rrBLUP_model <- mixed.solve(y = myY_train,
                                      Z = myGD_train,
                                      X = fix_train)

          pred_effects <- myGD_test %*% rrBLUP_model$u
          fix_effects <- fix_test  %*% rrBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)
        }else{
          myCV_fix  <- as.matrix(CV)
          myCV_fix=make_full_rank(myCV_fix)

          gBLUP_model <- mixed.solve(y = pheno_train,
                                     K = K,
                                     X=myCV_fix)
          pred_effects <- gBLUP_model$u[-c(1:length(train_phenotype[,2]))]
          fix_effects <- as.matrix(myCV_fix[-c(1:length(train_phenotype[,2])),])  %*% gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)

        }

        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)

      }

    }


    names(results)<- c("Pearson","Spearman","RMSE","R2","MAE")


  results_ALL=list(Results=results,Predictions=prediction)
  return(results_ALL)
}

rrBLUP_CV <- function(genotypes, phenotype,Kernel="Markers",PCA=NULL,CV=NULL,markers=NULL, folds = 5,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,transformation=NULL)
{

  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])
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
      # Calculate the GS model using rrBLUP
      myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
      myGD_train=apply(myGD_train,2,as.numeric)
      myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
      myGD_test=apply(myGD_test,2,as.numeric)
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


    if(transformation=="sqrt"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
      phenotype[-fold_indices,2] <-sqrt(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-sqrt(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="log"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
      phenotype[-fold_indices,2] <-log(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-log(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="boxcox"){
      phenotype[-fold_indices,2] <-boxcox_t(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-boxcox_t(phenotype[fold_indices,2])
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }

    if(transformation=="none"){
      myY_train <- phenotype[fold_indices,2]
      myY_test <- phenotype[-fold_indices,2]
      pheno_train <- phenotype[,2]
      pheno_train[-fold_indices] <- NA
    }


    if(length(PCA)!=0){
      myPCA_train <- PCA[fold_indices,]
      myPCA_test <- PCA[-fold_indices,]
    }
    gc()

    if(length(CV)==0){

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
          pred_effects <- gBLUP_model$u[-fold_indices]
          fix_effects <- as.matrix(PCA[-fold_indices,])  %*% gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)
        }
        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)
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
          pred_effects <- gBLUP_model$u[-fold_indices]
          fix_effects <- gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)

        }
        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)

      }
    }else{

      #CV=impute(as.factor(CV))
      if(length(ncol(CV))==0){
        myCV_train <- CV[fold_indices]
        myCV_test <- CV[-fold_indices]
      }else{
        myCV_train <- CV[fold_indices,]
        myCV_test <- CV[-fold_indices,]
      }



      if(!is.null(PCA)){

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
          if(ncol(data.frame(fix_test_PC))==1){
            fix_test_PC=matrix(fix_test_PC)
          }
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
          pred_effects_PC <- gBLUP_model$u[-fold_indices]
          fix_effects_PC <- as.matrix(fix_PC[-fold_indices,])  %*% gBLUP_model$beta
          predictions <- c(pred_effects_PC) + c(fix_effects_PC)
        }
        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)
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
          fix_test=fix_test[,colnames(fix_train)]
          if(ncol(data.frame(fix_test))==1){
            fix_test=matrix(fix_test)
          }

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
          pred_effects <- gBLUP_model$u[-fold_indices]
          fix_effects <- as.matrix(myCV_fix[-fold_indices,])  %*% gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)

        }

        acc <- cor(predictions, myY_test, use = "pairwise.complete")
        sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=predictions,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)

      }

    }

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

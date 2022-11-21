rrBLUP_GAGS_CV <- function(genotypes, phenotype,Kernel="Markers",Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05,threshold=NULL, markers=NULL, folds = 5,Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,QTN=10,PCA.total=3,transformation=NULL)
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
      if(!is.null(markers)){
        samp=sample(2:ncol(GD), markers)
        m_samp=GD[,samp]
        train_GM<-train_GM[samp,]
        train_GD <- m_samp[fold_indices,]
        #myGD_test <- m_samp[-fold_indices,]
      }else{
        train_GD <- GD[fold_indices,]
        train_GM<-GM
        #myGD_test <- GD[-fold_indices,]
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
        prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=NA,RE=NA,FE=NA)

        Predictions<-prediction
        BGLR_acc_results[[i]] <- list(results)

        Predictions_ALL=rbind(Predictions_ALL,Predictions)

      }else{
        sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

        myCV_train <- myGD_train[,sm]
        myCV_test <- myGD_test[,sm]
        myCV<-as.matrix(genotypes[,sm])
        if(!is.null(PCA)){
          myPCA_train <- PCA[fold_indices,]
          myPCA_test <- PCA[-fold_indices,]
        }

        if(!is.null(PCA)){

          if(Kernel=="Markers"){
            fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
            fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
            #fix_PC=as.matrix(cbind(myCV,PCA))
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
            predictions <- c(pred_effects_PC) + c(fix_effects_PC)
          }else{
            fix_PC=as.matrix(cbind(myCV,PCA))
            #p <- ncol(fix_PC)
            #XtX <- crossprod(fix_PC, fix_PC)
            #rank.X <- qr(XtX)$rank
            #if (rank.X < p) {
              #sm2=findLinearCombos(fix_PC)$remove
              #fix_PC=fix_PC[,-sm2]}

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

            #fix_myCV  <- as.matrix(myCV)
            #p <- ncol(fix_myCV)
            #XtX <- crossprod(fix_myCV, fix_myCV)
            #rank.X <- qr(XtX)$rank
            #if (rank.X < p) {
            #sm2=findLinearCombos(fix_myCV)$remove
            #myCV=myCV[,-sm2]}

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
            fix_myCV  <- as.matrix(myCV)
            fix_myCV=make_full_rank(fix_myCV)
            gBLUP_model <- mixed.solve(y = pheno_train,
                                       K = K,
                                       X=fix_myCV)
            pred_effects <- gBLUP_model$u[-fold_indices]
            fix_effects <- as.matrix(fix_myCV[-fold_indices,])  %*% gBLUP_model$beta
            predictions <- c(pred_effects) + c(fix_effects)
          }
          acc <- cor(predictions, myY_test, use = "pairwise.complete")
          sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
          metrics=postResample(pred=predictions,obs=myY_test)
          results=c(ACC=acc,SACC=sacc,metrics)
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)

        }
        Predictions<-prediction
        BGLR_acc_results[[i]] <- list(results)
        Predictions_ALL=rbind(Predictions_ALL,Predictions)
      }



    }else{
      top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
      GWASSM=top10[1:QTN,]$SNP
      myCV_train <- myGD_train[,GWASSM]
      myCV_test <- myGD_test[,GWASSM]
      myCV<-as.matrix(genotypes[,GWASSM])



      if(!is.null(PCA)){

        gc()
        if(Kernel=="Markers"){
          myPCA_train <- PCA[fold_indices,]
          myPCA_test <- PCA[-fold_indices,]
          fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
          fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
          #fix_PC=as.matrix(cbind(myCV,PCA))

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
          predictions <- c(pred_effects_PC) + c(fix_effects_PC)
        }else{
          fix_PC=as.matrix(cbind(myCV,PCA))
          #p <- ncol(fix_PC)
          #XtX <- crossprod(fix_PC, fix_PC)
          #rank.X <- qr(XtX)$rank
          #if (rank.X < p) {
            #sm2=findLinearCombos(fix_PC)$remove
            #fix_PC=fix_PC[,-sm2]}
          fix_PC=make_full_rank(fix_PC)
          gBLUP_model <- mixed.solve(y = pheno_train,
                                     K = K,
                                     X=fix_PC)
          pred_effects <- gBLUP_model$u[-fold_indices]
          fix_effects <- as.matrix(fix_PC[-fold_indices,])  %*% gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)
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
          #myCV_fix  <- as.matrix(myCV)
          #p <- ncol(myCV_fix)
          #XtX <- crossprod(myCV_fix, myCV_fix)
          #rank.X <- qr(XtX)$rank
          #if (rank.X < p) {
          #sm2=findLinearCombos(myCV_fix)$remove
          #myCV=myCV[,-sm2]}
          rrBLUP_model <- mixed.solve(y = myY_train,
                                      Z = myGD_train,
                                      X = fix_train)

          pred_effects <- myGD_test %*% rrBLUP_model$u
          fix_effects <- fix_test  %*% rrBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)
        }else{
          myCV_fix  <- as.matrix(myCV)
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

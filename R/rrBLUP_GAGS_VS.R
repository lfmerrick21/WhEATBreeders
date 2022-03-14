
rrBLUP_GAGS_VS <- function(train_genotypes, train_phenotype,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,Kernel="Markers",Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,GWAS="BLINK",alpha=0.05,threshold=NULL, QTN=10,markers=NULL,PCA.total=3,transformation=NULL)
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
        prediction=data.frame(test_phenotype,GEBV=NA,RE=NA,FE=NA)


      }else{
        sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

        myCV_train <- myGD_train[,sm]
        myCV_test <- myGD_test[,sm]
        myCV<-as.matrix(genotypes[,sm])

        if(!is.null(train_PCA)){

          if(Kernel=="Markers"){
            myPCA_train <- train_PCA
            myPCA_test <- test_PCA
            PCA<-rbind(train_PCA,test_PCA)
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

            #fix_myCV  <- as.matrix(myCV)
            #p <- ncol(fix_myCV)
            #XtX <- crossprod(fix_myCV, fix_myCV)
            #rank.X <- qr(XtX)$rank
            #if (rank.X < p) {
            #sm2=findLinearCombos(fix_myCV)$remove
            #myCV=myCV[,-sm2]}

            fix_train=make_full_rank(fix_train)
            fix_test=fix_test_PC[,colnames(fix_train)]
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
            pred_effects <- gBLUP_model$u[-c(1:length(train_phenotype[,2]))]
            fix_effects <- as.matrix(fix_myCV[-c(1:length(train_phenotype[,2])),])  %*% gBLUP_model$beta
            predictions <- c(pred_effects) + c(fix_effects)
          }
          acc <- cor(predictions, myY_test, use = "pairwise.complete")
          sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
          metrics=postResample(pred=predictions,obs=myY_test)
          results=c(ACC=acc,SACC=sacc,metrics)
          prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)

        }
      }



    }else{
      top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
      GWASSM=top10[1:QTN,]$SNP
      myCV_train <- myGD_train[,GWASSM]
      myCV_test <- myGD_test[,GWASSM]
      myCV<-as.matrix(genotypes[,GWASSM])



      if(!is.null(train_PCA)){

        gc()
        if(Kernel=="Markers"){
          myPCA_train <- train_PCA
          myPCA_test <- test_PCA
          PCA<-rbind(train_PCA,test_PCA)
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
          pred_effects <- gBLUP_model$u[-c(1:length(train_phenotype[,2]))]
          fix_effects <- as.matrix(fix_PC[-c(1:length(train_phenotype[,2])),])  %*% gBLUP_model$beta
          predictions <- c(pred_effects) + c(fix_effects)
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
          fix_test=fix_test[,colnames(fix_train)]

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

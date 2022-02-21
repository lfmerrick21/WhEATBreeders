Caret_Models_CV <- function(genotypes, phenotype,model="svmRadial",type="Regression", folds = 5,markers=NULL,Kernel="Gaussian",Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,transformation=NULL,sampling="up",repeats=5,method="repeatedcv"){
  # for fitting SVMs
  #library(DMwR)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  svm_results <- list()
  svm_results_metrics<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP


    if(Kernel=="Markers"){

      if(type=="Classification"){
        phenotype[,2]=as.factor(phenotype[,2])
        phenotype=droplevels(phenotype)
        myY_train <- phenotype[fold_indices,2]
        myY_test <- phenotype[-fold_indices,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)
      }else{
        if(transformation=="sqrt"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
          myY_train <- sqrt(phenotype[fold_indices,2])
          myY_test <- sqrt(phenotype[-fold_indices,2])
          pheno_train <- sqrt(phenotype[,2])
          pheno_train[-fold_indices] <- NA
        }

        if(transformation=="log"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
          myY_train <- log(phenotype[fold_indices,2])
          myY_test <- log(phenotype[-fold_indices,2])
          pheno_train <- log(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(phenotype[fold_indices,2])
          myY_test <- boxcox_t(phenotype[-fold_indices,2])
          pheno_train <- boxcox_t(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="none"){
          myY_train <- phenotype[fold_indices,2]
          myY_test <- phenotype[-fold_indices,2]
          pheno_train <- phenotype[,2]
          pheno_train[-fold_indices] <- NA
        }
      }

      if(!is.null(markers)){
        samp=sample(1:ncol(genotypes), markers)
        m_samp=genotypes[,samp]
        myGD_train <- m_samp[fold_indices,]
        myGD_test <- m_samp[-fold_indices,]
      }else{
        myGD_train <- genotypes[fold_indices,]
        myGD_test <- genotypes[-fold_indices,]
      }

      myGD_train2=as.matrix(sapply(myGD_train, as.numeric))
      myGD_test1=as.matrix(sapply(myGD_test, as.numeric))
      if(type=="Classification"){
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none",sampling = sampling)
      }else{
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none")
      }
      svmFit1 <- train(myGD_train2,myY_train,
                       method = model,
                       preProcess = c("center", "scale","nzv"),
                       trControl = ctrl,
                       tuneLength = 10,
                       allowParallel=TRUE)
      svm.linear_pred1 <- predict(svmFit1, myGD_test1)

      Predictions_ALL[[i]]=list(Fold=rep(i,length(phenotype[-fold_indices,])),Genotypes=phenotype[-fold_indices,1],obs=myY_test,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
        svm_results_metrics[[i]] <- list(mets)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)

      }

    }

    else if(Kernel=="VanRaden"){
      #pheno_train <- phenotype
      #pheno_train[-fold_indices,2] <- NA
      #myY_train <- pheno_train[,2]
      #myY_train=droplevels(myY_train)

      VR=GAPIT.kinship.VanRaden(genotypes)


      if(type=="Classification"){
        phenotype[,2]=as.factor(phenotype[,2])
        phenotype=droplevels(phenotype)
        myY_train <- phenotype[fold_indices,2]
        myY_test <- phenotype[-fold_indices,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)
      }else{
        if(transformation=="sqrt"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
          myY_train <- sqrt(phenotype[fold_indices,2])
          myY_test <- sqrt(phenotype[-fold_indices,2])
          pheno_train <- sqrt(phenotype[,2])
          pheno_train[-fold_indices] <- NA
        }

        if(transformation=="log"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
          myY_train <- log(phenotype[fold_indices,2])
          myY_test <- log(phenotype[-fold_indices,2])
          pheno_train <- log(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(phenotype[fold_indices,2])
          myY_test <- boxcox_t(phenotype[-fold_indices,2])
          pheno_train <- boxcox_t(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="none"){
          myY_train <- phenotype[fold_indices,2]
          myY_test <- phenotype[-fold_indices,2]
          pheno_train <- phenotype[,2]
          pheno_train[-fold_indices] <- NA
        }
      }

      myGD_train2 <- VR[fold_indices,]
      myGD_test1 <- VR[-fold_indices,]

      if(type=="Classification"){
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none",sampling = sampling)
      }else{
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none")
      }
      svmFit1 <- train(myGD_train2,myY_train,
                       method = model,
                       preProcess = c("center", "scale","nzv"),
                       trControl = ctrl,
                       tuneLength = 10,
                       allowParallel=TRUE)
      svm.linear_pred1 <- predict(svmFit1,myGD_test1)


      Predictions_ALL[[i]]=list(Fold=rep(i,length(phenotype[-fold_indices,])),Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
        svm_results_metrics[[i]] <- list(mets)
      }else{
        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)

      }

    }else if(Kernel=="PC"){
      PCA=prcomp(genotypes)
      myPCA=PCA$x


      if(type=="Classification"){
        phenotype[,2]=as.factor(phenotype[,2])
        phenotype=droplevels(phenotype)
        myY_train <- phenotype[fold_indices,2]
        myY_test <- phenotype[-fold_indices,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)
      }else{
        if(transformation=="sqrt"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
          myY_train <- sqrt(phenotype[fold_indices,2])
          myY_test <- sqrt(phenotype[-fold_indices,2])
          pheno_train <- sqrt(phenotype[,2])
          pheno_train[-fold_indices] <- NA
        }

        if(transformation=="log"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
          myY_train <- log(phenotype[fold_indices,2])
          myY_test <- log(phenotype[-fold_indices,2])
          pheno_train <- log(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(phenotype[fold_indices,2])
          myY_test <- boxcox_t(phenotype[-fold_indices,2])
          pheno_train <- boxcox_t(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="none"){
          myY_train <- phenotype[fold_indices,2]
          myY_test <- phenotype[-fold_indices,2]
          pheno_train <- phenotype[,2]
          pheno_train[-fold_indices] <- NA
        }
      }

      myGD_train2 <- myPCA[fold_indices,]
      myGD_test1 <- myPCA[-fold_indices,]

      if(type=="Classification"){
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none",sampling = sampling)
      }else{
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none")
      }
      svmFit1 <- train(myGD_train2,myY_train,
                       method = model,
                       preProcess = c("center", "scale","nzv"),
                       trControl = ctrl,
                       tuneLength = 10,
                       allowParallel=TRUE)
      svm.linear_pred1 <- predict(svmFit1,myGD_test1)


      Predictions_ALL[[i]]=list(Fold=rep(i,length(phenotype[-fold_indices,])),Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
        svm_results_metrics[[i]] <- list(mets)
      }else{
        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)

      }

    }

    else if(Kernel=="Zhang"){
      #pheno_train <- phenotype
      #pheno_train[-fold_indices,2] <- NA
      #myY_train <- pheno_train[,2]
      #myY_train=droplevels(myY_train)
      ZZ=GAPIT.kinship.Zhang(genotypes)


      if(type=="Classification"){
        phenotype[,2]=as.factor(phenotype[,2])
        phenotype=droplevels(phenotype)
        myY_train <- phenotype[fold_indices,2]
        myY_test <- phenotype[-fold_indices,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)
      }else{
        if(transformation=="sqrt"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
          myY_train <- sqrt(phenotype[fold_indices,2])
          myY_test <- sqrt(phenotype[-fold_indices,2])
          pheno_train <- sqrt(phenotype[,2])
          pheno_train[-fold_indices] <- NA
        }

        if(transformation=="log"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
          myY_train <- log(phenotype[fold_indices,2])
          myY_test <- log(phenotype[-fold_indices,2])
          pheno_train <- log(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(phenotype[fold_indices,2])
          myY_test <- boxcox_t(phenotype[-fold_indices,2])
          pheno_train <- boxcox_t(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="none"){
          myY_train <- phenotype[fold_indices,2]
          myY_test <- phenotype[-fold_indices,2]
          pheno_train <- phenotype[,2]
          pheno_train[-fold_indices] <- NA
        }
      }

      myGD_train2 <- ZZ[fold_indices,]
      myGD_test1 <- ZZ[-fold_indices,]

      if(type=="Classification"){
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none",sampling = sampling)
      }else{
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none")
      }
      svmFit1 <- train(myGD_train2,myY_train,
                       method = model,
                       preProcess = c("center", "scale","nzv"),
                       trControl = ctrl,
                       tuneLength = 10,
                       allowParallel=TRUE)
      svm.linear_pred1 <- predict(svmFit1,myGD_test1)

      Predictions_ALL[[i]]=list(Fold=rep(i,length(phenotype[-fold_indices,])),Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
        svm_results_metrics[[i]] <- list(mets)
      }else{
        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)

      }
    }

    else if(Kernel=="Endelman"){
      #pheno_train <- phenotype
      #pheno_train[-fold_indices,2] <- NA
      #myY_train <- pheno_train[,2]
      #myY_train=droplevels(myY_train)
      G=A.mat(genotypes)


      if(type=="Classification"){
        phenotype[,2]=as.factor(phenotype[,2])
        phenotype=droplevels(phenotype)
        myY_train <- phenotype[fold_indices,2]
        myY_test <- phenotype[-fold_indices,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)
      }else{
        if(transformation=="sqrt"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
          myY_train <- sqrt(phenotype[fold_indices,2])
          myY_test <- sqrt(phenotype[-fold_indices,2])
          pheno_train <- sqrt(phenotype[,2])
          pheno_train[-fold_indices] <- NA
        }

        if(transformation=="log"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
          myY_train <- log(phenotype[fold_indices,2])
          myY_test <- log(phenotype[-fold_indices,2])
          pheno_train <- log(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(phenotype[fold_indices,2])
          myY_test <- boxcox_t(phenotype[-fold_indices,2])
          pheno_train <- boxcox_t(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="none"){
          myY_train <- phenotype[fold_indices,2]
          myY_test <- phenotype[-fold_indices,2]
          pheno_train <- phenotype[,2]
          pheno_train[-fold_indices] <- NA
        }
      }

      myGD_train2 <- G[fold_indices,]
      myGD_test1 <- G[-fold_indices,]

      if(type=="Classification"){
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none",sampling = sampling)
      }else{
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none")
      }
      svmFit1 <- train(myGD_train2,myY_train,
                       method = model,
                       preProcess = c("center", "scale","nzv"),
                       trControl = ctrl,
                       tuneLength = 10,
                       allowParallel=TRUE)
      svm.linear_pred1 <- predict(svmFit1,myGD_test1)

      Predictions_ALL[[i]]=list(Fold=rep(i,length(phenotype[-fold_indices,])),Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
        svm_results_metrics[[i]] <- list(mets)
      }else{
        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)

      }
    }
    else{
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




      if(type=="Classification"){
        phenotype[,2]=as.factor(phenotype[,2])
        phenotype=droplevels(phenotype)
        myY_train <- phenotype[fold_indices,2]
        myY_test <- phenotype[-fold_indices,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)
      }else{
        if(transformation=="sqrt"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
          myY_train <- sqrt(phenotype[fold_indices,2])
          myY_test <- sqrt(phenotype[-fold_indices,2])
          pheno_train <- sqrt(phenotype[,2])
          pheno_train[-fold_indices] <- NA
        }

        if(transformation=="log"){
          phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
          myY_train <- log(phenotype[fold_indices,2])
          myY_test <- log(phenotype[-fold_indices,2])
          pheno_train <- log(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(phenotype[fold_indices,2])
          myY_test <- boxcox_t(phenotype[-fold_indices,2])
          pheno_train <- boxcox_t(phenotype[,2])
          pheno_train[-fold_indices] <- NA

        }

        if(transformation=="none"){
          myY_train <- phenotype[fold_indices,2]
          myY_test <- phenotype[-fold_indices,2]
          pheno_train <- phenotype[,2]
          pheno_train[-fold_indices] <- NA
        }
      }

      myGD_train2 <- K[fold_indices,]
      myGD_test1 <- K[-fold_indices,]

      if(type=="Classification"){
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none",sampling = sampling)
      }else{
        ctrl=trainControl(method = method, repeats = repeats, savePredictions="none")
      }
      svmFit1 <- train(myGD_train2,myY_train,
                       method = model,
                       preProcess = c("center", "scale","nzv"),
                       trControl = ctrl,
                       tuneLength = 10,
                       allowParallel=TRUE)
      svm.linear_pred1 <- predict(svmFit1,myGD_test1)

      Predictions_ALL[[i]]=list(Fold=rep(i,length(phenotype[-fold_indices,])),Genotypes=phenotype[-fold_indices,1],obs=phenotype[-fold_indices,2],pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(phenotype[,2])))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        phenotype[,2] <- factor(phenotype[,2], levels = mynames)

        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
        svm_results_metrics[[i]] <- list(mets)
      }else{
        acc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete")
        sacc <- cor(as.numeric(phenotype[-fold_indices,2]), as.numeric(svm.linear_pred1), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=phenotype[-fold_indices,2])
        results=c(ACC=acc,SACC=sacc,metrics)
        svm_results[[i]] <- list(results)

      }
    }

  }
  if(type=="Classification"){
    model_vect <- c("Pearson","Spearman","R2","Kappa")
  }else{
    model_vect <- c("Pearson","Spearman","RMSE","R2","MAE")
  }
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  if(type=="Classification"){
    results=colMeans(acc_fold, na.rm = TRUE)[2:5]
    results_ALL=list(Results=results,CM=svm_results_metrics,Predictions=Predictions_ALL)
  }else{
    results=colMeans(acc_fold, na.rm = TRUE)[2:6]
    results_ALL=list(Results=results,Predictions=Predictions_ALL)
  }
  return(results_ALL)
}

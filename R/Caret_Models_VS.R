Caret_Models_VS <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,model="svmRadial",type="Regression",markers=NULL,Kernel="Gaussian",Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,transformation=NULL,sampling="up",repeats=5,method="repeatedcv"){

  # for fitting SVMs
  #library(DMwR)
  # Make the CV list
    # Split into training and testing data
    # Calculate the GS model using rrBLUP


    if(Kernel=="Markers"){

      if(type=="Classification"){
        train_phenotype[,2] <- as.factor(train_phenotype[,2])
        test_phenotype[,2] <- as.factor(test_phenotype[,2])

        train_phenotype[train_phenotype=="NaN"]<-NA
        train_phenotype=droplevels(train_phenotype)

        test_phenotype[test_phenotype=="NaN"]<-NA
        test_phenotype=droplevels(test_phenotype)

        myY_train <- train_phenotype[,2]
        myY_test <- test_phenotype[,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)

        pheno_test=test_phenotype
        pheno_test[,2]<-NA
        pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

      }else{
        if(transformation=="sqrt"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)

          myY_train <- sqrt(train_phenotype[,2])
          myY_test <- sqrt(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
        }

        if(transformation=="log"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
          myY_train <- log(train_phenotype[,2])
          myY_test <- log(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(train_phenotype[,2])
          myY_test <- boxcox_t(test_phenotype[,2])

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
      }

      if(!is.null(markers)){
        samp=sample(1:ncol(train_genotypes), markers)
        myGD_train <- train_genotypes[,samp]
        myGD_test <- test_genotypes[,samp]
      }else{
        myGD_train <- train_genotypes
        myGD_test <- test_genotypes
      }

      myGD_train2=as.matrix(sapply(train_genotypes, as.numeric))
      myGD_test1=as.matrix(sapply(test_genotypes, as.numeric))
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

      prediction=data.frame(test_phenotype,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)

      }

    }

    else if(Kernel=="VanRaden"){
      #pheno_train <- phenotype
      #pheno_train[-fold_indices,2] <- NA
      #myY_train <- pheno_train[,2]
      #myY_train=droplevels(myY_train)
      genotypes=rbind(train=train_genotypes,test=test_genotypes)
      VR=GAPIT.kinship.VanRaden(genotypes)


      if(type=="Classification"){
        train_phenotype[,2] <- as.factor(train_phenotype[,2])
        test_phenotype[,2] <- as.factor(test_phenotype[,2])
        train_phenotype[train_phenotype=="NaN"]<-NA
        train_phenotype=droplevels(train_phenotype)

        test_phenotype[test_phenotype=="NaN"]<-NA
        test_phenotype=droplevels(test_phenotype)

        myY_train <- train_phenotype[,2]
        myY_test <- test_phenotype[,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)

        pheno_test=test_phenotype
        pheno_test[,2]<-NA
        pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

      }else{
        if(transformation=="sqrt"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)

          myY_train <- sqrt(train_phenotype[,2])
          myY_test <- sqrt(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
        }

        if(transformation=="log"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
          myY_train <- log(train_phenotype[,2])
          myY_test <- log(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(train_phenotype[,2])
          myY_test <- boxcox_t(test_phenotype[,2])

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
      }

      myGD_train2=VR[c(1:length(train_phenotype[,2])),]
      myGD_test1=VR[-c(1:length(train_phenotype[,2])),]

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

      prediction=data.frame(test_phenotype,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)

      }

    }else if(Kernel=="PC"){
      genotypes=rbind(train=train_genotypes,test=test_genotypes)
      PCA=prcomp(genotypes)
      myPCA=PCA$x


      if(type=="Classification"){
        train_phenotype[,2] <- as.factor(train_phenotype[,2])
        test_phenotype[,2] <- as.factor(test_phenotype[,2])
        train_phenotype[train_phenotype=="NaN"]<-NA
        train_phenotype=droplevels(train_phenotype)

        test_phenotype[test_phenotype=="NaN"]<-NA
        test_phenotype=droplevels(test_phenotype)

        myY_train <- train_phenotype[,2]
        myY_test <- test_phenotype[,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)

        pheno_test=test_phenotype
        pheno_test[,2]<-NA
        pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

      }else{
        if(transformation=="sqrt"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)

          myY_train <- sqrt(train_phenotype[,2])
          myY_test <- sqrt(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
        }

        if(transformation=="log"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
          myY_train <- log(train_phenotype[,2])
          myY_test <- log(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(train_phenotype[,2])
          myY_test <- boxcox_t(test_phenotype[,2])

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
      }

      myGD_train2=myPCA[c(1:length(train_phenotype[,2])),]
      myGD_test1=myPCA[-c(1:length(train_phenotype[,2])),]

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

      prediction=data.frame(test_phenotype,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)

      }

    }

    else if(Kernel=="Zhang"){
      #pheno_train <- phenotype
      #pheno_train[-fold_indices,2] <- NA
      #myY_train <- pheno_train[,2]
      #myY_train=droplevels(myY_train)
      genotypes=rbind(train=train_genotypes,test=test_genotypes)
      ZZ=GAPIT.kinship.Zhang(genotypes)


      if(type=="Classification"){
        train_phenotype[,2] <- as.factor(train_phenotype[,2])
        test_phenotype[,2] <- as.factor(test_phenotype[,2])
        train_phenotype[train_phenotype=="NaN"]<-NA
        train_phenotype=droplevels(train_phenotype)

        test_phenotype[test_phenotype=="NaN"]<-NA
        test_phenotype=droplevels(test_phenotype)

        myY_train <- train_phenotype[,2]
        myY_test <- test_phenotype[,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)

        pheno_test=test_phenotype
        pheno_test[,2]<-NA
        pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

      }else{
        if(transformation=="sqrt"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)

          myY_train <- sqrt(train_phenotype[,2])
          myY_test <- sqrt(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
        }

        if(transformation=="log"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
          myY_train <- log(train_phenotype[,2])
          myY_test <- log(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(train_phenotype[,2])
          myY_test <- boxcox_t(test_phenotype[,2])

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
      }

      myGD_train2=ZZ[c(1:length(train_phenotype[,2])),]
      myGD_test1=ZZ[-c(1:length(train_phenotype[,2])),]

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

      prediction=data.frame(test_phenotype,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)

      }
    }

    else if(Kernel=="Endelman"){
      #pheno_train <- phenotype
      #pheno_train[-fold_indices,2] <- NA
      #myY_train <- pheno_train[,2]
      #myY_train=droplevels(myY_train)
      genotypes=rbind(train=train_genotypes,test=test_genotypes)
      G=A.mat(genotypes)


      if(type=="Classification"){
        train_phenotype[,2] <- as.factor(train_phenotype[,2])
        test_phenotype[,2] <- as.factor(test_phenotype[,2])
        train_phenotype[train_phenotype=="NaN"]<-NA
        train_phenotype=droplevels(train_phenotype)

        test_phenotype[test_phenotype=="NaN"]<-NA
        test_phenotype=droplevels(test_phenotype)

        myY_train <- train_phenotype[,2]
        myY_test <- test_phenotype[,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)

        pheno_test=test_phenotype
        pheno_test[,2]<-NA
        pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

      }else{
        if(transformation=="sqrt"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)

          myY_train <- sqrt(train_phenotype[,2])
          myY_test <- sqrt(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
        }

        if(transformation=="log"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
          myY_train <- log(train_phenotype[,2])
          myY_test <- log(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(train_phenotype[,2])
          myY_test <- boxcox_t(test_phenotype[,2])

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
      }

      myGD_train2=G[c(1:length(train_phenotype[,2])),]
      myGD_test1=G[-c(1:length(train_phenotype[,2])),]

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

      prediction=data.frame(test_phenotype,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)

      }
    }
    else{
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




      if(type=="Classification"){
        train_phenotype[,2] <- as.factor(train_phenotype[,2])
        test_phenotype[,2] <- as.factor(test_phenotype[,2])
        train_phenotype[train_phenotype=="NaN"]<-NA
        train_phenotype=droplevels(train_phenotype)

        test_phenotype[test_phenotype=="NaN"]<-NA
        test_phenotype=droplevels(test_phenotype)

        myY_train <- train_phenotype[,2]
        myY_test <- test_phenotype[,2]
        myY_train=as.factor(myY_train)
        myY_train=droplevels(myY_train)
        myY_test=as.factor(myY_test)

        pheno_test=test_phenotype
        pheno_test[,2]<-NA
        pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

      }else{
        if(transformation=="sqrt"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] < 0, 0)

          myY_train <- sqrt(train_phenotype[,2])
          myY_test <- sqrt(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
        }

        if(transformation=="log"){
          train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
          test_phenotype[,2]=replace(test_phenotype[,2], test_phenotype[,2] <= 0, 0.000001)
          myY_train <- log(train_phenotype[,2])
          myY_test <- log(test_phenotype[,2])

          pheno_test=test_phenotype
          pheno_test[,2]<-NA
          pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

        }

        if(transformation=="boxcox"){
          myY_train <- boxcox_t(train_phenotype[,2])
          myY_test <- boxcox_t(test_phenotype[,2])

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
      }

      myGD_train2=K[c(1:length(train_phenotype[,2])),]
      myGD_test1=K[-c(1:length(train_phenotype[,2])),]

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

      prediction=data.frame(test_phenotype,pred=svm.linear_pred1)


      if(type=="Classification"){
        mynames <- unique(c(levels(svm.linear_pred1), levels(myY_test)))
        svm.linear_pred1 <- factor(svm.linear_pred1, levels = mynames)
        myY_test <- factor(myY_test, levels = mynames)
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)
        mets=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
      }else{
        acc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete")
        sacc <- cor(as.numeric(svm.linear_pred1), as.numeric(myY_test), use = "pairwise.complete", method = c("spearman"))
        metrics=postResample(pred=svm.linear_pred1,obs=myY_test)
        results=c(ACC=acc,SACC=sacc,metrics)

      }
    }



  if(type=="Classification"){
    names(results) <- c("Pearson","Spearman","R2","Kappa")
    results_ALL=list(Results=results,CM=mets,Predictions=prediction)
  }else{
    names(results) <- c("Pearson","Spearman","RMSE","R2","MAE")
    results_ALL=list(Results=results,Predictions=prediction)
  }
  return(results_ALL)
}

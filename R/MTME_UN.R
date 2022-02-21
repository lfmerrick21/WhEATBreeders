MTME_UN <- function(matrix,model="MT",trait=c("IT","SEV"),nIter = 80000, burnIn = 10000, folds = 5)
{

  # Make the CV list
  fold_list <- make_CV_sets(length(matrix$BUN$Y[,1]), k = folds)

  if(model=="ST"){
    Predictions_ALL=list()
    BGLR_acc_results_ALL=list()
    # Split into training and testing data
    for(j in 1:length(trait)){
      BGLR_acc_results <- c()
      Predictions<-c()
      for (i in 1:length(fold_list)){
        fold_indices <- which(fold_list[[i]])
        # Split into training and testing data
        phenotype=matrix$BUN$Y[,c("Genotype","ENV",trait)]
        Y <- as.matrix(matrix$BUN$Y[,-c(1,2)])
        Y[-fold_indices,] <- NA
        K <- matrix$BUN$K
        #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
        #K_expanded=ZG%*%K%*%t(ZG)
        Z_L=matrix$BUN$Z_L
        Z_E=matrix$BUN$Z_E
        K_expanded=matrix$BUN$K_expanded
        K_GE=matrix$BUN$K_GE
        #Calculate the GS model using BGLR
        ##RKHS
        ST_IT_ETA<-list(G=list(K=K_expanded,model='RKHS'))
        ST_IT_model_results <- BGLR(y = Y[,j], ETA = ST_IT_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('st_IT_',nIter,'_',burnIn,'_'))
        ST_IT_predictions <- predict(ST_IT_model_results)

        #acc_IT <- cor(phenotype[-fold_indices,j], ST_IT_predictions[-fold_indices], use = "pairwise.complete")
        #sacc_IT <- cor(phenotype[-fold_indices,j], ST_IT_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
        #metrics_IT=postResample(pred=ST_IT_predictions[-fold_indices],obs=phenotype[-fold_indices,1])
        #results_IT=c(ACC=acc_IT,SACC=sacc_IT,metrics_IT)


        #model_results<-list(IT=ST_IT_model_results,SEV=ST_SEV_model_results)
        #results<-list(results_IT)
        #names(results)<-trait[j]
        prediction=data.frame(Fold=rep(i,nrow(phenotype[-fold_indices,])),Trait=rep(trait[j],nrow(phenotype[-fold_indices,])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$ENV[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=ST_IT_predictions[-fold_indices])
        #prediction=prediction %>% mutate_if(is.character,as.factor)
        #prediction$Fold=as.factor(prediction$Fold)
        results=prediction%>%
          group_by(Fold,Env)%>%
          summarise(ACC=cor(Obs, Pred, use = "pairwise.complete"),
                    SACC=cor(Obs, Pred, use = "pairwise.complete", method = c("spearman")),
                    Rsquare=R2(pred=Pred,obs=Obs,na.rm=TRUE),
                    RMSE=RMSE(pred=Pred,obs=Obs,na.rm=TRUE),
                    MAE=MAE(pred=Pred,obs=Obs,na.rm=TRUE))
        Predictions=rbind(Predictions,prediction)
        BGLR_acc_results=rbind(BGLR_acc_results,results)
      }
      Predictions_ALL[[j]]=Predictions
      BGLR_acc_results_ALL[[j]]=BGLR_acc_results
    }
    names(Predictions_ALL)=trait
    names(BGLR_acc_results_ALL)=trait
    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 4))
    for(j in 1:length(trait)){
      try=BGLR_acc_results_ALL[[j]]
      results_long <- data.frame(try[,names(try) %in% c("Fold")],Trait=rep(trait[j], nrow(try)),try[,!names(try) %in% c("Fold")] )%>%pivot_longer(!c(Trait,Fold,Env), names_to = "Model", values_to = "Accuracy")
      BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
    }

    names(BGLR_acc_table) <- c("Fold","Trait","Env","Model", "r")
    BGLR_acc_table=BGLR_acc_table %>% mutate_if(is.character,as.factor)
    results_ALL=BGLR_acc_table%>%
      group_by(Trait,Env,Model)%>%
      summarise(Accuracy=mean(r,na.rm = TRUE))
    results_ST=list(Results=results_ALL,Predictions=Predictions_ALL)
    return(results_ST)
  }

  if(model=="MT"){
    Predictions_ALL=list()
    BGLR_acc_results_ALL=list()
    for (i in 1:length(fold_list)){
      fold_indices <- which(fold_list[[i]])
      phenotype=matrix$BUN$Y[,c("Genotype","ENV",trait)]
      Y <- as.matrix(matrix$BUN$Y[,-c(1,2)])
      Y[-fold_indices,] <- NA
      K <- matrix$BUN$K
      #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
      #K_expanded=ZG%*%K%*%t(ZG)
      Z_L=matrix$BUN$Z_L
      Z_E=matrix$BUN$Z_E
      K_expanded=matrix$BUN$K_expanded
      K_GE=matrix$BUN$K_GE
      # Calculate the GS model using BGLR
      ##RKHS
      gc()
      MT_ETA<-list(list(K=K_expanded,model="RKHS"))
      MT_model_results <- Multitrait(y = Y, ETA = MT_ETA, nIter=nIter, burnIn=burnIn)
      MT_predictions <- MT_model_results$ETA[[1]]$u
      Predictions=list()
      BGLR_acc_results=list()
      for(j in 1:length(trait)){
        prediction=data.frame(Fold=rep(i,nrow(phenotype[-fold_indices,])),Trait=rep(trait[j],nrow(phenotype[-fold_indices,])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$ENV[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=MT_predictions[-fold_indices,j])
        results=prediction%>%
          group_by(Fold,Env)%>%
          summarise(ACC=cor(Obs, Pred, use = "pairwise.complete"),
                    SACC=cor(Obs, Pred, use = "pairwise.complete", method = c("spearman")),
                    Rsquare=R2(pred=Pred,obs=Obs,na.rm=TRUE),
                    RMSE=RMSE(pred=Pred,obs=Obs,na.rm=TRUE),
                    MAE=MAE(pred=Pred,obs=Obs,na.rm=TRUE))
        Predictions[[j]]=prediction
        BGLR_acc_results[[j]]=results
      }
      names(Predictions)=trait
      names(BGLR_acc_results)=trait
      Predictions_ALL[[i]]=Predictions
      BGLR_acc_results_ALL[[i]]=BGLR_acc_results
      #names(Predictions_ALL)=trait
      #names(BGLR_acc_results_ALL)=trait
    }

    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 4))
    for (i in 1:length(BGLR_acc_results_ALL)){
      for(j in 1:length(trait)){
        try=BGLR_acc_results_ALL[[i]][[j]]
        results_long <- data.frame(try[,names(try) %in% c("Fold")],Trait=rep(trait[j], nrow(try)),try[,!names(try) %in% c("Fold")] )%>%pivot_longer(!c(Trait,Fold,Env), names_to = "Model", values_to = "Accuracy")
        BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
      }
    }

    names(BGLR_acc_table) <- c("Fold","Trait","Env","Model", "r")
    BGLR_acc_table=BGLR_acc_table %>% mutate_if(is.character,as.factor)
    results_ALL=BGLR_acc_table%>%
      group_by(Trait,Env,Model)%>%
      summarise(accuracy=mean(r,na.rm = TRUE))
    results_MT=list(Results=results_ALL,Predictions=Predictions_ALL)
    return(results_MT)
  }

  if(model=="BME"){
    Predictions_ALL=list()
    BGLR_acc_results_ALL=list()
    # Split into training and testing data
    for(j in 1:length(trait)){
      BGLR_acc_results <- c()
      Predictions<-c()
      for (i in 1:length(fold_list)){
        fold_indices <- which(fold_list[[i]])
        # Split into training and testing data
        phenotype=matrix$BUN$Y[,c("Genotype","ENV",trait)]
        Y <- as.matrix(matrix$BUN$Y[,-c(1,2)])
        Y[-fold_indices,] <- NA
        K <- matrix$BUN$K
        #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
        #K_expanded=ZG%*%K%*%t(ZG)
        Z_L=matrix$BUN$Z_L
        Z_E=matrix$BUN$Z_E
        K_expanded=matrix$BUN$K_expanded
        K_GE=matrix$BUN$K_GE
        #Calculate the GS model using BGLR
        ##RKHS

        #ST_IT_ETA<-list(G=list(K=K_expanded,model='RKHS'))

        #Z_E=matrix$MTME$Z_E
        #K_expanded=matrix$MTME$K_expanded
        #K_GE=matrix$MTME$K_GE

        BME_IT_ETA=list(Env=list(model="FIXED",X=Z_E),Lines=list(model="RKHS",K=K_expanded))

        nEnv=length(levels(as.factor(phenotype$ENV)))

        for(k in 1:nEnv){
          BME_IT_ETA[[k+2]]<-list(K=matrix$BUN$KUN[[k]],model='RKHS')
        }

        BME_IT_ETA <- BGLR(y = Y[,j], ETA = BME_IT_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('st_IT_',nIter,'_',burnIn,'_'))
        BME_IT_predictions <- predict(BME_IT_ETA)

        #acc_IT <- cor(phenotype[-fold_indices,j], BME_IT_predictions[-fold_indices], use = "pairwise.complete")
        #sacc_IT <- cor(phenotype[-fold_indices,j], BME_IT_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
        #metrics_IT=postResample(pred=BME_IT_predictions[-fold_indices],obs=phenotype[-fold_indices,1])
        #results_IT=c(ACC=acc_IT,SACC=sacc_IT,metrics_IT)


        #model_results<-list(IT=ST_IT_model_results,SEV=ST_SEV_model_results)
        #results<-list(results_IT)
        #names(results)<-trait[j]
        prediction=data.frame(Fold=rep(i,nrow(phenotype[-fold_indices,])),Trait=rep(trait[j],nrow(phenotype[-fold_indices,])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$ENV[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=BME_IT_predictions[-fold_indices])
        #prediction=prediction %>% mutate_if(is.character,as.factor)
        #prediction$Fold=as.factor(prediction$Fold)
        results=prediction%>%
          group_by(Fold,Env)%>%
          summarise(ACC=cor(Obs, Pred, use = "pairwise.complete"),
                    SACC=cor(Obs, Pred, use = "pairwise.complete", method = c("spearman")),
                    Rsquare=R2(pred=Pred,obs=Obs,na.rm=TRUE),
                    RMSE=RMSE(pred=Pred,obs=Obs,na.rm=TRUE),
                    MAE=MAE(pred=Pred,obs=Obs,na.rm=TRUE))
        Predictions=rbind(Predictions,prediction)
        BGLR_acc_results=rbind(BGLR_acc_results,results)
      }
      Predictions_ALL[[j]]=Predictions
      BGLR_acc_results_ALL[[j]]=BGLR_acc_results
    }
    names(Predictions_ALL)=trait
    names(BGLR_acc_results_ALL)=trait
    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 4))
    for(j in 1:length(trait)){
      try=BGLR_acc_results_ALL[[j]]
      results_long <- data.frame(try[,names(try) %in% c("Fold")],Trait=rep(trait[j], nrow(try)),try[,!names(try) %in% c("Fold")] )%>%pivot_longer(!c(Trait,Fold,Env), names_to = "Model", values_to = "Accuracy")
      BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
    }

    names(BGLR_acc_table) <- c("Fold","Trait","Env","Model", "r")
    BGLR_acc_table=BGLR_acc_table %>% mutate_if(is.character,as.factor)
    results_ALL=BGLR_acc_table%>%
      group_by(Trait,Env,Model)%>%
      summarise(Accuracy=mean(r,na.rm = TRUE))
    results_BME=list(Results=results_ALL,Predictions=Predictions_ALL)
    return(results_BME)
  }

  if(model=="BMTME"){
    Predictions_ALL=list()
    BGLR_acc_results_ALL=list()
    for (i in 1:length(fold_list)){
      fold_indices <- which(fold_list[[i]])
      phenotype=matrix$BUN$Y[,c("Genotype","ENV",trait)]
      Y <- as.matrix(matrix$BUN$Y[,-c(1,2)])
      Y[-fold_indices,] <- NA
      K <- matrix$BUN$K
      #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
      #K_expanded=ZG%*%K%*%t(ZG)
      Z_L=matrix$BUN$Z_L
      Z_E=matrix$BUN$Z_E
      K_expanded=matrix$BUN$K_expanded
      K_GE=matrix$BUN$K_GE
      # Calculate the GS model using BGLR
      ##RKHS
      gc()

      #Z_E=matrix$MTME$Z_E
      #K_expanded=matrix$MTME$K_expanded
      #K_GE=matrix$MTME$K_GE

      BMTME_ETA=list(Env=list(model="FIXED",X=Z_E),Lines=list(model="RKHS",K=K_expanded))

      nEnv=length(levels(as.factor(phenotype$ENV)))

      for(k in 1:nEnv){
        BMTME_ETA[[k+2]]<-list(K=matrix$BUN$KUN[[k]],model='RKHS')
      }

      BMTME_model_results= Multitrait(y = Y, ETA=BMTME_ETA,nIter =nIter, burnIn = burnIn)


      BMTME_predictions <- BMTME_model_results$ETAHat
      Predictions=list()
      BGLR_acc_results=list()
      for(j in 1:length(trait)){
        prediction=data.frame(Fold=rep(i,nrow(phenotype[-fold_indices,])),Trait=rep(trait[j],nrow(phenotype[-fold_indices,])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$ENV[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=BMTME_predictions[-fold_indices,j])
        results=prediction%>%
          group_by(Fold,Env)%>%
          summarise(ACC=cor(Obs, Pred, use = "pairwise.complete"),
                    SACC=cor(Obs, Pred, use = "pairwise.complete", method = c("spearman")),
                    Rsquare=R2(pred=Pred,obs=Obs,na.rm=TRUE),
                    RMSE=RMSE(pred=Pred,obs=Obs,na.rm=TRUE),
                    MAE=MAE(pred=Pred,obs=Obs,na.rm=TRUE))
        Predictions[[j]]=prediction
        BGLR_acc_results[[j]]=results
      }
      names(Predictions)=trait
      names(BGLR_acc_results)=trait
      Predictions_ALL[[i]]=Predictions
      BGLR_acc_results_ALL[[i]]=BGLR_acc_results
      #names(Predictions_ALL)=trait
      #names(BGLR_acc_results_ALL)=trait
    }

    BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 4))
    for (i in 1:length(BGLR_acc_results_ALL)){
      for(j in 1:length(trait)){
        try=BGLR_acc_results_ALL[[i]][[j]]
        results_long <- data.frame(try[,names(try) %in% c("Fold")],Trait=rep(trait[j], nrow(try)),try[,!names(try) %in% c("Fold")] )%>%pivot_longer(!c(Trait,Fold,Env), names_to = "Model", values_to = "Accuracy")
        BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
      }
    }

    names(BGLR_acc_table) <- c("Fold","Trait","Env","Model", "r")
    BGLR_acc_table=BGLR_acc_table %>% mutate_if(is.character,as.factor)
    results_ALL=BGLR_acc_table%>%
      group_by(Trait,Env,Model)%>%
      summarise(accuracy=mean(r,na.rm = TRUE))
    results_BMTME=list(Results=results_ALL,Predictions=Predictions_ALL)
    return(results_BMTME)
  }



}

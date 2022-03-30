MTME <- function(matrix,model="MT",trait=c("IT","SEV"),nIter = 80000, burnIn = 10000, folds = 5,UN=FALSE)
{


  if(UN==TRUE)
  {
    fold_list_mtme <- make_CV_sets(length(matrix$MTME$Y[,1]), k = folds)
    if(model=="ST"){
      Predictions_ALL=list()
      BGLR_acc_results_ALL=list()

      phenotype=matrix$MTME$YUN

      K=matrix$MTME$K
      Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE

      # Split into training and testing data
      for(j in 1:length(trait)){
        BGLR_acc_results <- c()
        Predictions<-c()
        for (i in 1:length(fold_list_mtme)){
          fold_indices <- which(fold_list_mtme[[i]])
          # Split into training and testing data

          Y <- as.matrix(phenotype[,trait])
          Y[-fold_indices,] <- NA
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
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,trait[j]])),Trait=rep(trait[j],length(phenotype[-fold_indices,trait[j]])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$Env[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=ST_IT_predictions[-fold_indices])
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

      phenotype=matrix$MTME$YUN
      #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
      #K_expanded=ZG%*%K%*%t(ZG)
      K=matrix$MTME$K
      Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE

      for (i in 1:length(fold_list_mtme)){
        fold_indices <- which(fold_list_mtme[[i]])
        # Split into training and testing data
        Y <- as.matrix(phenotype[,trait])
        Y[-fold_indices,] <- NA
        # Calculate the GS model using BGLR
        ##RKHS
        gc()
        MT_ETA<-list(list(K=K_expanded,model="RKHS"))
        MT_model_results <- Multitrait(y = Y, ETA = MT_ETA, nIter=nIter, burnIn=burnIn)
        MT_predictions <- MT_model_results$ETA[[1]]$u
        Predictions=list()
        BGLR_acc_results=list()
        for(j in 1:length(trait)){
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,trait[j]])),Trait=rep(trait[j],length(phenotype[-fold_indices,trait[j]])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$Env[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=MT_predictions[-fold_indices,j])

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
      phenotype=matrix$MTME$YUN
      #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
      #K_expanded=ZG%*%K%*%t(ZG)
      K=matrix$MTME$K
      Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE
      for(j in 1:length(trait)){
        BGLR_acc_results <- c()
        Predictions<-c()
        for (i in 1:length(fold_list_mtme)){
          fold_indices <- which(fold_list_mtme[[i]])
          # Split into training and testing data

          Y <- as.matrix(phenotype[,trait])
          Y[-fold_indices,] <- NA


          BME_IT_ETA=list(Env=list(model="FIXED",X=Z_E),Lines=list(model="RKHS",K=K_expanded),GE=list(model="RKHS",K=K_GE))

          BME_IT_ETA <- BGLR(y = Y[,j], ETA = BME_IT_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('st_IT_',nIter,'_',burnIn,'_'))
          BME_IT_predictions <- predict(BME_IT_ETA)

          #acc_IT <- cor(phenotype[-fold_indices,j], BME_IT_predictions[-fold_indices], use = "pairwise.complete")
          #sacc_IT <- cor(phenotype[-fold_indices,j], BME_IT_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
          #metrics_IT=postResample(pred=BME_IT_predictions[-fold_indices],obs=phenotype[-fold_indices,1])
          #results_IT=c(ACC=acc_IT,SACC=sacc_IT,metrics_IT)


          #model_results<-list(IT=ST_IT_model_results,SEV=ST_SEV_model_results)
          #results<-list(results_IT)
          #names(results)<-trait[j]
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,trait[j]])),Trait=rep(trait[j],length(phenotype[-fold_indices,trait[j]])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$Env[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=BME_IT_predictions[-fold_indices])

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

      phenotype=matrix$MTME$YUN
      #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
      #K_expanded=ZG%*%K%*%t(ZG)
      K=matrix$MTME$K
      Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE

      for (i in 1:length(fold_list_mtme)){
        fold_indices <- which(fold_list_mtme[[i]])

        Y <- as.matrix(phenotype[,trait])
        Y[-fold_indices,] <- NA

        #Z_E=matrix$MTME$Z_E
        #K_expanded=matrix$MTME$K_expanded
        #K_GE=matrix$MTME$K_GE

        BMTME_ETA=list(Env=list(model="FIXED",X=Z_E),Lines=list(model="RKHS",K=K_expanded),GE=list(model="RKHS",K=K_GE))
        BMTME_model_results= Multitrait(y = Y, ETA=BMTME_ETA,nIter =nIter, burnIn = burnIn)


        BMTME_predictions <- BMTME_model_results$ETAHat
        Predictions=list()
        BGLR_acc_results=list()
        for(j in 1:length(trait)){
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,trait[j]])),Trait=rep(trait[j],length(phenotype[-fold_indices,trait[j]])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$Env[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=BMTME_predictions[-fold_indices,j])
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

  }else{
    fold_list_bmtme <- make_CV_sets(length(matrix$BMTME$Y[,1]), k = folds)
    fold_list_bme <- make_CV_sets(length(data.frame(matrix$BME$Y[[1]])[,1]), k = folds)
    fold_list_mtme <- make_CV_sets(length(matrix$MTME$Y[,1]), k = folds)
    if(model=="ST"){
      Predictions_ALL=list()
      BGLR_acc_results_ALL=list()
      phenotype=matrix$MTME$YUN

      K=matrix$MTME$K
      Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE
      # Split into training and testing data
      for(j in 1:length(trait)){
        BGLR_acc_results <- c()
        Predictions<-c()
        for (i in 1:length(fold_list_mtme)){
          fold_indices <- which(fold_list_mtme[[i]])
          # Split into training and testing data

          Y <- as.matrix(phenotype[,trait])
          Y[-fold_indices,] <- NA
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
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,trait[j]])),Trait=rep(trait[j],length(phenotype[-fold_indices,trait[j]])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$Env[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=ST_IT_predictions[-fold_indices])

          #acc_IT <- cor(phenotype[-fold_indices,j], ST_IT_predictions[-fold_indices], use = "pairwise.complete")
          #sacc_IT <- cor(phenotype[-fold_indices,j], ST_IT_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
          #metrics_IT=postResample(pred=ST_IT_predictions[-fold_indices],obs=phenotype[-fold_indices,1])
          #results_IT=c(ACC=acc_IT,SACC=sacc_IT,metrics_IT)


          #model_results<-list(IT=ST_IT_model_results,SEV=ST_SEV_model_results)
          #results<-list(results_IT)
          #names(results)<-trait[j]
          #prediction=data.frame(Fold=rep(i,length(names[-fold_indices])),Trait=rep(trait[j],length(names[-fold_indices])),Genotype=names2[-fold_indices],Env=Env_names2[-fold_indices],Obs=phenotype[-fold_indices,j],Pred=ST_IT_predictions[-fold_indices])
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
      phenotype=matrix$MTME$YUN
      #ZG <- model.matrix(~0 + as.factor(rownames(Y)))
      #K_expanded=ZG%*%K%*%t(ZG)
      K=matrix$MTME$K
      Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE
      for (i in 1:length(fold_list_mtme)){
        fold_indices <- which(fold_list_mtme[[i]])
        # Split into training and testing data
        Y <- as.matrix(phenotype[,trait])
        Y[-fold_indices,] <- NA
        # Calculate the GS model using BGLR
        ##RKHS
        gc()
        MT_ETA<-list(list(K=K_expanded,model="RKHS"))
        MT_model_results <- Multitrait(y = Y, ETA = MT_ETA, nIter=nIter, burnIn=burnIn)
        MT_predictions <- MT_model_results$ETA[[1]]$u
        Predictions=list()
        BGLR_acc_results=list()
        for(j in 1:length(trait)){
          prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,trait[j]])),Trait=rep(trait[j],length(phenotype[-fold_indices,trait[j]])),Genotype=phenotype$Genotype[-fold_indices],Env=phenotype$Env[-fold_indices],Obs=phenotype[-fold_indices,trait[j]],Pred=MT_predictions[-fold_indices,j])
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
        BGLR_acc_results <- list()
        Predictions<-c()
        for (i in 1:length(fold_list_bme)){
          fold_indices <- which(fold_list_bme[[i]])
          phenotype_IT=matrix$BME$Y[[trait[j]]]
          pheno_IT=phenotype_IT[,-1]
          Y_IT=as.matrix(matrix$BME$Y[[trait[j]]])
          Y_IT[-fold_indices,] <- NA
          Y_IT=as.matrix(Y_IT[,-1])
          Y_IT=apply(Y_IT,2,as.numeric)
          ZI_IT=matrix$BME$Z1[[trait[j]]]

          # Calculate the GS model using BGLR
          ##RKHS
          gc()

          BME_IT_model_results <- BME(Y = Y_IT, Z1 = ZI_IT, nIter = nIter, burnIn = burnIn)

          acc_IT <- cor(pheno_IT[-fold_indices,], BME_IT_model_results$yHat[-fold_indices,], use = "pairwise.complete")
          sacc_IT <- cor(pheno_IT[-fold_indices,], BME_IT_model_results$yHat[-fold_indices,], use = "pairwise.complete", method = c("spearman"))

          R2_IT_C=c()
          RMSE_IT_C=c()
          MAE_IT_C=c()
          for(k in 1:ncol(BME_IT_model_results$yHat)){
            R2_IT=R2(pred=BME_IT_model_results$yHat[-fold_indices,k],obs=pheno_IT[-fold_indices,k][[1]],na.rm=TRUE)
            RMSE_IT=RMSE(pred=BME_IT_model_results$yHat[-fold_indices,k],obs=pheno_IT[-fold_indices,k][[1]],na.rm=TRUE)
            MAE_IT=MAE(pred=BME_IT_model_results$yHat[-fold_indices,k],obs=pheno_IT[-fold_indices,k][[1]],na.rm=TRUE)
            R2_IT_C=c(R2_IT_C,R2_IT)
            RMSE_IT_C=c(RMSE_IT_C,RMSE_IT)
            MAE_IT_C=c(MAE_IT_C,MAE_IT)
          }
          acc_IT_list=as.list(as.data.frame(acc_IT))
          sacc_IT_list=as.list(as.data.frame(sacc_IT))

          acc_IT_list=diag(acc_IT)
          acc_IT_list=unname(acc_IT_list)
          sacc_IT_list=diag(sacc_IT)
          sacc_IT_list=unname(sacc_IT_list)
          #metrics_IT=list(RMSE=RMSE_IT_C,Rsquare=R2_IT_C,MAE=MAE_IT_C)
          results_IT=list(ACC=acc_IT_list,SACC=sacc_IT_list,RMSE=RMSE_IT_C,Rsquare=R2_IT_C,MAE=MAE_IT_C)

          #model_results<-list(IT=BME_IT_model_results,SEV=BME_SEV_model_results)
          prediction_Obs=data.frame(Fold=rep(i,nrow(phenotype_IT[-fold_indices,])),Trait=rep(trait[j],nrow(phenotype_IT[-fold_indices,])),
                                    Genotype=phenotype_IT[-fold_indices,]$Genotype,Obs=phenotype_IT[-fold_indices,-1])
          prediction_Pred=data.frame(Fold=rep(i,nrow(phenotype_IT[-fold_indices,])),Trait=rep(trait[j],nrow(phenotype_IT[-fold_indices,])),
                                     Genotype=phenotype_IT[-fold_indices,]$Genotype,Pred=BME_IT_model_results$yHat[-fold_indices,])
          prediction_Obs_long=prediction_Obs%>%
            pivot_longer(
              cols = starts_with("Obs"),
              names_to = "Env",
              names_prefix = c("Obs."),
              values_to = c("Obs"),
              values_drop_na = TRUE
            )

          prediction_Pred_long=prediction_Pred%>%
            pivot_longer(
              cols = starts_with("Pred"),
              names_to = "Env",
              names_prefix = c("Pred."),
              values_to = c("Pred"),
              values_drop_na = TRUE
            )

          prediction <- merge(prediction_Obs_long,prediction_Pred_long,by=c("Fold","Trait","Genotype","Env"))
          prediction$Env<-gsub(".", " ", prediction$Env, fixed=TRUE)

          Predictions=rbind(Predictions,prediction)
          BGLR_acc_results[[i]] <-results_IT
        }

        Predictions_ALL[[j]]=Predictions
        BGLR_acc_results_ALL[[j]]=BGLR_acc_results
      }
      names(Predictions_ALL)=trait
      BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 4))
      for(j in 1:length(trait)){
        for (i in 1:length(BGLR_acc_results_ALL[[j]])){
          try=unlist(BGLR_acc_results_ALL[[j]][[i]])
          model_vect=names(try)
          env=colnames(matrix$BME$Y[[j]])[-1]
          env_vect=rep(env, (length(env)*(length(model_vect)/length(env)))/length(env))
          results_long <- data.frame(rep(i, length(model_vect)),rep(trait[j], length(model_vect)),env_vect, model_vect, unlist(BGLR_acc_results_ALL[[j]][[i]]))
          BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
        }
      }

      names(BGLR_acc_table) <- c("Fold","Trait","Env","Model", "r")
      BGLR_acc_table=BGLR_acc_table %>% mutate_if(is.character,as.factor)
      results_ALL=BGLR_acc_table%>%
        group_by(Trait,Env,Model)%>%
        summarise(accuracy=mean(r,na.rm = TRUE))
      results_ALL$Model <- (str_extract(results_ALL$Model, "[aA-zZ]+"))
      results_ALL=results_ALL %>% mutate_if(is.character,as.factor)
      results_BME=list(Results=results_ALL,Predictions=Predictions_ALL)
      return(results_BME)
    }

    if(model=="BMTME"){
      Predictions_ALL=list()
      BGLR_acc_results_ALL=list()
      for (i in 1:length(fold_list_bmtme)){
        fold_indices <- which(fold_list_bmtme[[i]])
        # Split into training and testing data
        phenotype=matrix$BMTME$Y[,trait]
        Y=as.matrix(matrix$BMTME$Y[,trait])
        Y[-fold_indices,] <- NA
        X=matrix$BMTME$X
        Z1=matrix$BMTME$Z1
        Z2=matrix$BMTME$Z2

        names=c()
        for(k in 1:ncol(matrix$BME$Y$IT[,-1])){
          names=rbind(names,as.character(data.frame(matrix$BME$Y$IT)$Genotype))
        }


        if(ncol(matrix$BME$Y$IT[,-1])==1){
          Env_names=rep( colnames(matrix$BME$Y$IT[,-1]),nrow(matrix$BME$Y$IT))
        }else{
          Z_E=matrix$BMTME$X
          Z_E_vector=as.factor(apply(Z_E,1,function(foo){return(colnames(Z_E)[which.max(foo)])}))
          Z_E_vector=as.character(Z_E_vector)
          Env_names <- sapply(strsplit(Z_E_vector, split=')', fixed=TRUE), function(x) (x[2]))
        }

        Env_names2=as.vector(Env_names)
        names2=as.vector(names)
        rownames(phenotype)=names2
        rownames(Y)=names2
        # Calculate the GS model using BGLR
        ##RKHS
        BMTME_model_results=BMTME(Y = Y, X = X, Z1 = Z1, Z2 = Z2, nIter = nIter, burnIn = burnIn)
        #str(BMTME_model_results)
        #acc_BMTME<- cor(phenotype[-fold_indices,], BMTME_model_results$yHat[-fold_indices,], use = "pairwise.complete")
        #sacc_BMTME <- cor(phenotype[-fold_indices,], BMTME_model_results$yHat[-fold_indices,], use = "pairwise.complete", method = c("spearman"))

        #R2_BMTME_C=c()
        #RMSE_BMTME_C=c()
        #MAE_BMTME_C=c()
        #for(i in 1:ncol(BMTME_model_results$yHat)){
        #R2_BMTME=R2(pred=BMTME_model_results$yHat[-fold_indices,i],obs=phenotype[-fold_indices,i],na.rm=TRUE)
        #RMSE_BMTME=RMSE(pred=BMTME_model_results$yHat[-fold_indices,i],obs=phenotype[-fold_indices,i],na.rm=TRUE)
        #MAE_BMTME=MAE(pred=BMTME_model_results$yHat[-fold_indices,i],obs=phenotype[-fold_indices,i],na.rm=TRUE)
        #R2_BMTME_C=c(R2_BMTME_C,R2_BMTME)
        #RMSE_BMTME_C=c(RMSE_BMTME_C,RMSE_BMTME)
        #MAE_BMTME_C=c(MAE_BMTME_C,MAE_BMTME)
        #}
        #metrics_IT=list(RMSE=RMSE_BMTME_C[1],Rsquare=R2_BMTME_C[1],MAE=MAE_BMTME_C[1])
        #results_IT=list(ACC=acc_BMTME[1,1],SACC=sacc_BMTME[1,1],metrics_IT)

        #metrics_BMTME=list(RMSE=RMSE_BMTME_C,Rsquare=R2_BMTME_C,MAE=MAE_BMTME_C)
        #results_BMTME=list(ACC=acc_BMTME,SACC=sacc_BMTME,metrics_BMTME)

        #acc_BMTME_list=as.list(as.data.frame(acc_BMTME))
        #sacc_BMTME_list=as.list(as.data.frame(sacc_BMTME))

        #acc_IT_list=diag(acc_BMTME)
        #acc_IT_list=unname(acc_BMTME_list)
        #sacc_IT_list=diag(sacc_BMTME)
        #sacc_IT_list=unname(acc_BMTME_list)

        #results_IT=list(ACC=acc_IT_list,SACC=sacc_IT_list,RMSE=RMSE_IT_C,Rsquare=R2_IT_C,MAE=MAE_IT_C)

        #model_results<-BMTME_model_results
        #results<-list(results_BMTME)
        #results<-list(IT=results_IT,SEV=results_SEV)
        #prediction=data.frame(Genotype=names[-fold_indices],Env=Env_names[-fold_indices],obs=phenotype[-fold_indices,],pred=BMTME_model_results$yHat[-fold_indices,])

        Predictions=list()
        BGLR_acc_results=list()
        for(j in 1:length(trait)){
          prediction=data.frame(Fold=rep(i,length(names[-fold_indices])),Trait=rep(trait[j], length(names[-fold_indices])),Genotype=names2[-fold_indices],Env=Env_names2[-fold_indices],Obs=phenotype[-fold_indices,j],Pred=BMTME_model_results$yHat[-fold_indices,j])
          results=prediction%>%
            group_by(Fold,Trait,Env)%>%
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
          results_long <- data.frame(try)%>%pivot_longer(!c(Trait,Fold,Env), names_to = "Model", values_to = "Accuracy")
          BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
        }
      }
      names(BGLR_acc_table) <- c("Fold","Trait","Env","Model", "r")
      BGLR_acc_table=BGLR_acc_table %>% mutate_if(is.character,as.factor)
      results_ALL=BGLR_acc_table%>%
        group_by(Trait,Env,Model)%>%
        summarise(Accuracy=mean(r,na.rm = TRUE))
      results_BMTME=list(Results=results_ALL,Predictions=Predictions_ALL)
      return(results_BMTME)
    }

  }

}

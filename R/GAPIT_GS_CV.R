GAPIT_GS_CV <- function(genotypes, phenotype,myGM,model="cBLUP",PCA.total=3,CV=NULL,kinship="VanRaden",markers=NULL, folds = 5,transformation=NULL)
{

  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])
    #If you want to sample markers
    if(!is.null(markers)){
      geno_mat=genotypes[,-1]
      samp=sample(1:ncol(geno_mat), markers)
      m_samp=geno_mat[,samp]
      myGD=cbind(genotypes[,1],m_samp)
    }else{
      myGD <- genotypes
    }

    # Split into training and testing data
    if(transformation=="sqrt"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] < 0, 0)
      phenotype[-fold_indices,2] <-sqrt(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-sqrt(phenotype[fold_indices,2])
      pheno_train <- phenotype
      pheno_train[-fold_indices,2] <- NA
      pheno_test=phenotype[-fold_indices,]
    }

    if(transformation=="log"){
      phenotype[,2]=replace(phenotype[,2], phenotype[,2] <= 0, 0.000001)
      phenotype[-fold_indices,2] <-log(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-log(phenotype[fold_indices,2])
      pheno_train <- phenotype
      pheno_train[-fold_indices,2] <- NA
      pheno_test=phenotype[-fold_indices,]

    }

    if(transformation=="boxcox"){
      phenotype[-fold_indices,2] <-boxcox_t(phenotype[-fold_indices,2])
      phenotype[fold_indices,2] <-boxcox_t(phenotype[fold_indices,2])
      pheno_train <- phenotype
      pheno_train[-fold_indices,2] <- NA
      pheno_test=phenotype[-fold_indices,]

    }

    if(transformation=="none"){
      pheno_train <- phenotype
      pheno_train[-fold_indices,2] <- NA
      pheno_test=phenotype[-fold_indices,]
    }
    #Make sure columns have same column names
    names(pheno_train)=c("Taxa", Trait)
    names(pheno_test)=c("Taxa", Trait)
    names(myGD)[1]=c("Taxa")
    if(!is.null(CV)){
      names(CV)[1]=c("Taxa")
    }

      myGAPIT=GAPIT(
        Y=pheno_train,
        GD=myGD,
        GM=myGM,
        CV=CV,
        kinship.algorithm=kinship,
        PCA.total=PCA.total,
        SNP.test=FALSE,
        model=model,
        file.output=FALSE)
      #Merge output
      gapit=merge(pheno_test,myGAPIT$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")

      acc <- cor(gapit[,2],gapit[,5], use = "pairwise.complete")
      sacc <- cor(gapit[,2],gapit[,5], use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=gapit[,5],obs=gapit[,2])
      results=c(ACC=acc,SACC=sacc,metrics)
      prediction=data.frame(Fold=rep(i,length(phenotype[-fold_indices,1])),phenotype[-fold_indices,1],Observed=gapit[,2],Predicted=gapit[,5])

        Predictions<-prediction
        BGLR_acc_results[[i]] <- list(results)

    Predictions_ALL=rbind(Predictions_ALL,Predictions)

    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
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

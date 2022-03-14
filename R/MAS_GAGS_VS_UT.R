
MAS_GAGS_VS_UT <- function(train_genotypes, train_phenotype,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05,threshold=NULL, QTN=10,markers=NULL,PCA.total=3,transformation=NULL)
{

  # Make the CV list

  genotypes<-rbind(train_genotypes,test_genotypes)
  if(!is.null(PCA)){
    PCA<-rbind(train_PCA,test_PCA)
  }

  # Split into training and testing data
  if(transformation=="sqrt"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)
    train_phenotype[,2] <-sqrt(train_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="log"){
    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)
    train_phenotype[,2] <-log(train_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="boxcox"){
    train_phenotype[,2] <-boxcox_t(train_phenotype[,2])
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])
  }

  if(transformation=="none"){
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- c(train=train_phenotype[,2],test=pheno_test[,2])

  }


  GWASR<- GAPIT(Y = train_phenotype,
                GD = train_GM,
                GM = train_GD,
                PCA.total=PCA.total,
                model = GWAS,
                file.output=F)
  if(threshold=="Bonferonni"){

    GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
    if(length(GWASSM)==0){

      prediction=data.frame(pheno_test,GEBV=NA)

    }else{
      sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

      CV <- genotypes[,GWASSM]


      MAS_train <- data.frame(CV)
      MAS_test  <- data.frame(CV)

      if(!is.null(train_PCA)){
        myPCA_train <-PCA
        myPCA_test <- PCA
        MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
        MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))

        GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

        names(GLM_data_PC)[1] <- "Y"
        #Linear model to calculate effects
        #You can run all signficant markers at once to see cumulative effect
        MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
        predictions <- predict(MAS_model_PC, MAS_test_PC)
        prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])


      }else{
        GLM_data <- data.frame(pheno_train, MAS_train)

        names(GLM_data)[1] <- "Y"
        #Linear model to calculate effects
        #You can run all signficant markers at once to see cumulative effect
        MAS_model <- lm(Y ~ ., data = GLM_data)
        predictions <- predict(MAS_model, MAS_test)
        prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
      }

    }




  }else{
    top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
    GWASSM=top10[1:QTN,]$SNP
    CV <- genotypes[,GWASSM]


    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)




    if(!is.null(train_PCA)){
      myPCA_train <- PCA
      myPCA_test <- PCA

      MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
      MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))

      GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

      names(GLM_data_PC)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
      predictions <- predict(MAS_model_PC, MAS_test_PC)
      prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
    }else{
      GLM_data <- data.frame(pheno_train, MAS_train)

      names(GLM_data)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model <- lm(Y ~ ., data = GLM_data)
      predictions <- predict(MAS_model, MAS_test)
      prediction=data.frame(pheno_test,GEBV=predictions[-c(1:length(train_phenotype[,2]))])
    }


  }

  results_ALL=list(Predictions=prediction)
  return(results_ALL)
}

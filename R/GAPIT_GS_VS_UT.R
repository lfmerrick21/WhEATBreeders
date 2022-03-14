

GAPIT_GS_VS_UT <- function(train_genotypes, train_phenotype,train_GM=NULL,test_genotypes, test_phenotype,test_GM=NULL,model="cBLUP",PCA.total=3,train_CV=NULL,test_CV=NULL,kinship="VanRaden",markers=NULL, transformation=NULL)
{

  #If you want to sample markers
  if(!is.null(markers)){
    genotypes=rbind(train=train_genotypes,test=test_genotypes)
    geno_mat=genotypes[,-1]
    samp=sample(1:ncol(geno_mat), markers)
    m_samp=geno_mat[,samp]
    myGD=cbind(genotypes[,1],m_samp)
    names(myGD)[1]=c("Taxa")
  }else{
    myGD_train <- train_genotypes
    myGD_test <- test_genotypes
    myGD <- rbind(train=train_genotypes,test=test_genotypes)
    names(myGD)[1]=c("Taxa")
  }
  # Split into training and testing data

  if(transformation=="sqrt"){
    names(train_phenotype)=c("Taxa", Trait)
    names(test_phenotype)=c("Taxa", Trait)

    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] < 0, 0)

    train_phenotype[,2] <-sqrt(train_phenotype[,2])

    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- rbind(train=train_phenotype,test=pheno_test)
  }

  if(transformation=="log"){
    names(train_phenotype)=c("Taxa",Trait)
    names(test_phenotype)=c("Taxa", Trait)

    train_phenotype[,2]=replace(train_phenotype[,2], train_phenotype[,2] <= 0, 0.000001)

    train_phenotype[,2] <-log(train_phenotype[,2])

    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- rbind(train=train_phenotype,test=pheno_test)

  }

  if(transformation=="boxcox"){
    names(train_phenotype)=c("Taxa", Trait)
    names(test_phenotype)=c("Taxa", Trait)

    train_phenotype[,2] <-boxcox_t(train_phenotype[,2])

    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- rbind(train=train_phenotype,test=pheno_test)

  }

  if(transformation=="none"){
    names(train_phenotype)=c("Taxa", Trait)
    names(test_phenotype)=c("Taxa", Trait)
    pheno_test=test_phenotype
    pheno_test[,2]<-NA
    pheno_train <- rbind(train=train_phenotype,test=pheno_test)
  }
  #Make sure columns have same column names


  if(!is.null(train_CV)){
    CV <- rbind(train=train_CV,test=test_CV)
    names(CV)[1]=c("Taxa")
  }

  dim(pheno_train)
  dim(myGD)

  myGAPIT=GAPIT(
    Y=pheno_train,
    GD=myGD,
    GM=train_GM,
    #CV=CV,
    #kinship.algorithm=kinship,
    PCA.total=PCA.total,
    model=model,
    file.output=FALSE)


  inclosure.from=10000
  bin.from=100000
  bin.to=100000
  bin.by=10000
  inclosure.to=10000
  inclosure.by=100
  sangwich.top=NULL
  sangwich.bottom=NULL
  SUPER_GS=T
  if(model=="gBLUP")
  { group.from=10000
  group.to=10000
  group.by=50}
  if(model=="cBLUP")
  { group.from=40
  group.to=10000
  group.by=40}
  if(model=="sBLUP")
  { group.from=10000
  group.to=10000
  inclosure.from=100
  bin.from=10000
  bin.to=100000
  bin.by=10000
  inclosure.to=400
  inclosure.by=100
  sangwich.top="MLM"
  sangwich.bottom="SUPER"
  SUPER_GS=T}

  myGAPIT=GAPIT(
    Y=pheno_train,
    GD=myGD,
    GM=train_GM,
    CV=CV,
    kinship.algorithm=kinship,
    group.from=group.from,group.to=group.to,group.by=group.by,
    PCA.total=PCA,
    model=model,
    SNP.test=FALSE,
    inclosure.from=inclosure.from,
    bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,
    inclosure.to=inclosure.to,
    inclosure.by=inclosure.by,
    sangwich.top=sangwich.top,
    sangwich.bottom=sangwich.bottom,
    SUPER_GS=SUPER_GS,
    file.output=FALSE)
  #Merge output
  gapit=merge(pheno_test,myGAPIT$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")

  prediction=data.frame(Genotype=pheno_test[,1],Observed=gapit[,2],Predicted=gapit[,5])

  results_ALL=list(Predictions=prediction)
  return(results_ALL)
}

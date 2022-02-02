GE_Matrix_IND <- function(genotypes, phenotype,trait,GE=FALSE,UN=TRUE,Kernel="Gaussian",Sparse=FALSE,m=NULL,degree=NULL, nL=NULL,model="BMTME"){
  library(BGLR)
  library(BMTME)
  library(rrBLUP)
  library(tidyr)

  calc_maf_apply <- function(gt_mat, encoding = c(-1, 0, 1))
  {
    col_func1 <- function(gt_col)
    {
      allele1_ct <- (sum(gt_col == -1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
      allele2_ct <- (sum(gt_col == 1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)

      maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
    }

    col_func2 <- function(gt_col)
    {
      allele1_ct <- (sum(gt_col == 0, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
      allele2_ct <- (sum(gt_col == 2, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)

      maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
    }

    if (all(encoding == c(-1, 0, 1)))
    {
      maf_vect <- apply(gt_mat, 2, col_func1)
    } else if (all(encoding == c(0, 1, 2)))
    {
      maf_vect <- apply(gt_mat, 2, col_func2)
    } else{
      print('Encoding not recognized, returning NULL')
      maf_vect <- NULL
    }

    return(maf_vect)
  }

  Pheno=phenotype
  Pheno=droplevels(Pheno)
  Pheno=complete(Pheno, Genotype, ENV)
  maf <- calc_maf_apply(genotypes, encoding = c(0, 1, 2))
  mono_indices <- which(maf ==0)
  if(length(mono_indices)!=0){
    genotypes = genotypes[,-mono_indices]
  }
  if(GE==TRUE){

    if(Kernel=="Endelman"){
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
      K=A.mat(X)
    }else{
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
    if(model=="BMTME"){
      #G=K
      #diag(G)=diag(G)+1/1e4
      #L=t(chol(G))
      LG <- cholesky(K)
      ZG <- model.matrix(~0 + as.factor(Pheno$Genotype))
      Z.G <- ZG %*% LG
      Z.E <- model.matrix(~0 + as.factor(Pheno$ENV))
      ZEG <- model.matrix(~0 + as.factor(Pheno$Genotype):as.factor(Pheno$ENV))
      G2 <- kronecker(diag(length(unique(Pheno$ENV))), data.matrix(K))
      LG2 <- cholesky(G2)
      Z.EG <- ZEG %*% LG2
      Y <- as.matrix(Pheno[,trait])
      YUN=phenotype[,c("Genotype","ENV",trait)]
      YUN=droplevels(YUN)
      BMTME_Matrix=list(Y=Y,YUN=YUN,K=K,LG=LG,X=Z.E,Z1=Z.G,Z2=Z.EG)
    }
    if(model=="MTME"){
      Y <- as.matrix(Pheno[,trait])
      YUN=phenotype[,c("Genotype","ENV",trait)]
      YUN=droplevels(YUN)

      Z_L=model.matrix(~0+Genotype,data=YUN)
      dim(Z_L)
      Z_E=model.matrix(~0+ENV,data=YUN)
      dim(Z_E)
      K_E=Z_E%*%t(Z_E)

      K_expanded=Z_L%*%K%*%t(Z_L)
      K_GE=K_expanded*K_E


      MTME_Matrix=list(Y=Y,YUN=YUN,K=K,Z_L=Z_L,Z_E=Z_E,K_E=K_E,K_expanded=K_expanded,K_GE=K_GE)
    }
    if(UN==TRUE){
      if(Kernel=="Endelman"){
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
        K=A.mat(X)
      }else{
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
      if(model=="MEI"){
        YUN=phenotype[,c("Genotype","ENV",trait)]
        YUN=droplevels(YUN)
        env=YUN$ENV
        env=droplevels(env)
        #MEI Kernel
        Z_L=model.matrix(~0+Genotype,data=YUN)
        dim(Z_L)
        Z_E=model.matrix(~0+ENV,data=YUN)
        dim(Z_E)
        K_E=Z_E%*%t(Z_E)
        K_expanded=Z_L%*%K%*%t(Z_L)
        K_GE=K_expanded*K_E
        #MDe
        KUN=list()
        for(k in 1:length(levels(env))){
          ZEE<-matrix(0,nrow=nrow(Z_E),ncol=ncol(Z_E))
          ZEE[,k]<-Z_E[,k]
          ZEEZ<-(ZEE%*%t(Z_E))
          K3<-K_expanded*ZEEZ
          KUN[[paste0(levels(env)[k])]]=K3
        }
        #MEI SNP Explicitly
        X=as.matrix(genotypes)
        X=apply(genotypes,2,as.numeric)
        #sum(rowSums(is.na(X)))
        XS<-scale(X,center=TRUE,scale=TRUE)
        #sum(rowSums(is.na(XS)))
        library(caret)
        nzv <- nearZeroVar(XS)
        XSZV <- XS[, -nzv]
        #sum(rowSums(is.na(XSZV)))
        #K=A.mat(XSZV)
        X=XSZV/sqrt(ncol(XSZV))
        rownames(X)<-rownames(genotypes)
        y=YUN[,3]
        names(y)=YUN$Genotype
        env=YUN$ENV
        env=droplevels(env)
        X0=X[names(y),] # Matrix for main effects
        stopifnot(all(rownames(X0)==names(y)))

        XUN=list()
        # now interactions
        for(i in 1:length(levels(env))){
          X_trait=X0
          for(j in 1:nrow(X0)){
            X_trait[j,]<-(env[j]==levels(env)[i])*X0[j,]
          }
          XUN[[paste0(levels(env)[i])]]=X_trait
        }
        BUN_Matrix=list(Y=YUN,X0=X0,XUN=XUN,KUN=KUN,K=K,Z_L=Z_L,Z_E=Z_E,K_E=K_E,K_expanded=K_expanded,K_GE=K_GE,KUN)
      }
    }
  }else{
    if(Kernel=="Endelman"){
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
      K=A.mat(X)
    }else{
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
    if(model=="BMTME"){
      #G=K
      #diag(G)=diag(G)+1/1e4
      #L=t(chol(G))
      LG <- cholesky(K)
      ZG <- model.matrix(~0 + as.factor(Pheno$Genotype))
      Z.G <- ZG %*% LG
      Y <- as.matrix(Pheno[,trait])
      YUN=phenotype[,c("Genotype","ENV",trait)]
      YUN=droplevels(YUN)
      BMTME_Matrix=list(Y=Y,YUN=YUN,K=K,LG=LG,Z1=Z.G)
    }
    if(model=="MTME"){
      Y <- as.matrix(Pheno[,trait])
      YUN=phenotype[,c("Genotype","ENV",trait)]
      YUN=droplevels(YUN)
      Z_L=model.matrix(~0+Genotype,data=YUN)
      #dim(Z_L)
      #Z_E=model.matrix(~0+ENV,data=Pheno)
      #dim(Z_E)
      #K_E=Z_E%*%t(Z_E)

      K_expanded=Z_L%*%K%*%t(Z_L)
      #K_GE=K_expanded*K_E


      MTME_Matrix=list(Y=Y,YUN=YUN,K=K,Z_L=Z_L,K_expanded=K_expanded)

    }


    if(UN==TRUE){
      if(Kernel=="Endelman"){
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
        K=A.mat(X)
      }else{
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
      if(model=="MEI"){
        YUN=phenotype[,c("Genotype","ENV",trait)]
        YUN=droplevels(YUN)
        X=apply(genotypes,2,as.numeric)
        #sum(rowSums(is.na(X)))
        XS<-scale(X,center=TRUE,scale=TRUE)
        #sum(rowSums(is.na(XS)))
        library(caret)
        nzv <- nearZeroVar(XS)
        XSZV <- XS[, -nzv]
        #sum(rowSums(is.na(XSZV)))
        #K=A.mat(XSZV)
        X=XSZV/sqrt(ncol(XSZV))
        rownames(X)<-rownames(genotypes)

        BUN_Matrix=list(Y=YUN,K=K,X=X)
      }

    }
  }
  if(model=="BME"){
    Y=list()
    Z1=list()
    for(i in 1:length(trait)){
      Trait_Pheno=Pheno[,c("Genotype","ENV",paste0(trait[i]))] %>% spread(ENV,c(paste0(trait[i])))
      if(Kernel=="Endelman"){
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
        K=A.mat(X)
      }else{
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
      #G=K
      #diag(G)=diag(G)+1/1e4
      #L=t(chol(G))
      LG <- cholesky(K)
      Trait_ZG <- model.matrix(~0 + as.factor(Trait_Pheno$Genotype))
      Trait_ZG <- Trait_ZG %*% LG
      Y[[paste0(trait[i])]]=Trait_Pheno
      Z1[[paste0(trait[i])]]=Trait_ZG
      YUN=phenotype[,c("Genotype","ENV",trait)]
      YUN=droplevels(YUN)
    }
    BME_Matrix=list(Y=Y,YUN,K=K,LG=LG,Z1=Z1)

  }


  if(model=="BMTME")
  {
    results=list(BMTME=BMTME_Matrix)
  }
  if(model=="MTME")
  {
    results=list(MTME=MTME_Matrix)
  }
  if(model=="BME")
  {
    results=list(BME=BME_Matrix)
  }
  if(model=="MEI")
  {
    results=list(BUN=BUN_Matrix)
  }
  return(results)
}

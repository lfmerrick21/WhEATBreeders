PGEandCV <- function(Pheno_Em_bl=NULL,myGD_EM=NULL,mytaxa_EM=NULL,myGM_EM=NULL,CV=NULL){
  colnames(Pheno_Em_bl)[1]<-c("Genotype")
  #rownames(Pheno_Em_bl) <- Pheno_Em_bl$Genotype
  Pheno_Em_bl=droplevels(Pheno_Em_bl)
  Markers_Em_bl<-myGD_EM

  colnames(mytaxa_EM)[1]<-c("Genotype")
  rownames(Markers_Em_bl) <- mytaxa_EM$Genotype
  colnames(Markers_Em_bl) <- myGM_EM$SNP

  colnames(CV)[1]<-c("Genotype")
  rownames(CV) <- CV$Genotype

  for(i in 2:ncol(CV)){
    CV[,i]=impute(CV[,i])
    CV[,i]=as.numeric(CV[,i])
  }

  Pheno_Em_bl <- Pheno_Em_bl[Pheno_Em_bl$Genotype %in% rownames(Markers_Em_bl),]
  Markers_Em_bl <- Markers_Em_bl[rownames(Markers_Em_bl) %in% Pheno_Em_bl$Genotype,]

  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]
  Pheno_Em_bl=droplevels(Pheno_Em_bl)

  CV <- CV[rownames(CV) %in% rownames(Pheno_Em_bl),]
  CV <- CV[order(rownames(CV)),]

  myGM<-myGM_EM[myGM_EM$SNP %in% colnames(Markers_Em_bl),]
  myGD<-data.frame(rownames(Markers_Em_bl),Markers_Em_bl)
  colnames(myGD)[1]<-c("taxa")

  PCA_bl_em=prcomp(Markers_Em_bl)
  myPCA_bl_em=PCA_bl_em$x
  PGPCA=list(pheno=Pheno_Em_bl,geno=Markers_Em_bl,map=myGM,numeric=myGD,PC=myPCA_bl_em,CV=CV)
  return(PGPCA)
}

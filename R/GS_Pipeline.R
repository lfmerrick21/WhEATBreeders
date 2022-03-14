Genotype<-fread("F:\\OneDrive\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Genomic Selection\\Genomic Selection Pipeline\\Jason_GBS\\WAC_2016-2020_production_filt.hmp.txt",fill=TRUE)
Phenotype<-read.csv("F:\\OneDrive\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Genomic Selection\\Genomic Selection Pipeline\\GS-Complex-Traits\\BL_EM_Pheno.csv",header=T)
Geno_Type="Hapmap"
Imputation="Beagle" #Beagle, KNN, or Middle
Missing_Rate=0.20
MAF=0.05
Scheme="K-Fold"
Method="Two-Step"
Type="Regression"
Model="rrBLUP"
Replications=1
Training="F5_2015"
Prediction="DH_2020"
CV=NULL
Trait="EM"
Study="Tutorial"
Outcome="Tested" #Tested or Untested
Trial=c("F5_2015","DH_2020","BL_2015_2020")
i=1
###############################################################################
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,bigmemory,biganalytics,ggplot2,tidyverse,dplyr,rrBLUP,knitr,cvTools,vcfR,compiler,BGLR,caret,gdata,WeightIt)


###############################################################################
Phenotype=Phenotype[,-1]

###############################################################################
Phenotype1=Phenotype %>%
  pivot_longer(!Taxa, names_to = "Env", values_to = "EM")
Phenotype1
Phenotype=data.frame(Name=l20d_IT_ap[,1],
                      ENV=rep(l20d[1,c(8)],length(l20d_IT_ap[,c(1)])),
                      Trial=rep(l20d[1,c(9)],length(l20d_IT_ap[,c(1)])),
                      Loc=rep(l20d[1,c(10)],length(l20d_IT_ap[,c(1)])),
                      Year=rep(l20d[1,c(11)],length(l20d_IT_ap[,c(1)])),
                      IT=l20d_IT_ap[,2],
                      SEV=l20d_SEV_ap[,2],
                      HD=l20d_hd_ap[,2],
                      PH=l20d_ph_ap[,2],
                      GY=l20d_kgha_ap[,2]
)
###############################################################################
WHEAT<-function(Phenotype,
                Genotype,
                QC=TRUE,
                GS=TRUE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Missing_Rate=0.20,
                MAF=0.05,
                #If QC Info is FALSE
                GDre=NULL,
                GT=NULL,
                GIre=NULL,
                #GS Info
                Type="Regression",
                Replications=1,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Tested", #Tested or Untested
                Trial=c("F5_2015","DH_2020","BL_2015_2020"),
                Scheme="K-Fold",
                Method="Two-Step",
                Package="rrBLUP",
                Model="rrBLUP",
                Kernel="Markers",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=NULL,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation=NULL,
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold=NULL,
                GE=TRUE,
                UN=FALSE,
                model="MTME")){
if(QC==TRUE){
#############PandG#########################
###############################################################################
#Get Taxa
names(Phenotype)[1]<-"Taxa"
Phenotype$Taxa=as.character(Phenotype$Taxa)
gname=as.vector(Phenotype$Taxa)
gname<-clean_names(gname)
#Filter Hapmap for Taxa
if(Geno_Type=="VCF"){
  genotypes_num<- VCF_2_numeric(Xbgl)
  GD=genotypes_num$genotypes
  names.use <- names(GD)[(names(GD) %in% gname)]
  GD <- GD[, names.use, with = FALSE]
  GI=genotypes_num$marker_map
  GT
  #Save Raw Files
  fwrite(GT,paste0(Study,"_GT_taxa.csv"))
  fwrite(GD,paste0(Study,"_GD_numeric_Raw.csv"))
  fwrite(GI,paste0(Study,"_GI_map_Raw.csv"))
  save(GD,GI,GT,file=paste0(Study,"_Raw.RData"))
}

if(Geno_Type=="Hapmap"){
  names.use <- names(Genotype)[(names(Genotype) %in% gname)]
  hapmap <- Genotype[, names.use, with = FALSE]
  hapmap=cbind(Genotype[,1:11],hapmap)
  #str(output)
  hapmap[hapmap=="NA"]<-"NA"
  hapmap$`assembly#`=as.character(hapmap$`assembly#`)
  hapmap$center=as.character(hapmap$center)
  hapmap$protLSID=as.character(hapmap$protLSID)
  hapmap$assayLSID=as.character(hapmap$assayLSID)
  hapmap$panelLSID=as.character(hapmap$panelLSID)
  hapmap$QCcode=as.character(hapmap$QCcode)
  dim(hapmap)
  #Save filtered hapmap
  save(hapmap,file=paste0(Study,"_Filtered_Hapmap.RData"))
  fwrite(hapmap,paste0(Study,"_Filtered_Hapmap.hmp.txt"),sep="\t",row.names = FALSE,col.names = TRUE)
  #Convert to Numeric
  if(Impuation=="Middle"){
    outG=GAPIT.HapMap(hapmap,SNP.impute="Middle")
  }else{
    outG=GAPIT.HapMap(hapmap,SNP.impute="None")
  }
  GT=outG$GT
  GT=data.frame(V1=rownames(GT),GT)
  GD=outG$GD
  GI=outG$GI
  rownames(GD)=rownames(GT)
  colnames(GD)=GI$SNP
  #Save Raw Files
  fwrite(GT,paste0(Study,"_GT_taxa.csv"))
  fwrite(GD,paste0(Study,"_GD_numeric_Raw.csv"))
  fwrite(GI,paste0(Study,"_GI_map_Raw.csv"))
  save(GD,GI,GT,file=paste0(Study,"_Raw.RData"))
}

#### Filter
##### Remove makers based on missing data
if(Filter=="None"){
  GDmf = GD
  GImf = GI
}else{
mr=calc_missrate(GD)
mr_indices <- which(mr > Missing_Rate)
if(length(mr_indices)!=0){
  GDmr = GD[,-mr_indices]
  GImr = GI[-mr_indices,]
}
dim(GDmr)
dim(GImr)

##### Remove Monomorphic Markers
maf <- calc_maf_apply(GDmr, encoding = c(0, 1, 2))
mono_indices <- which(maf ==0)
if(length(mono_indices)!=0){
  GDmo = GDmr[,-mono_indices]
  GImo = GImr[-mono_indices,]
}
dim(GDmo)
dim(GImo)

##### Remove MAF
maf <- calc_maf_apply(GDmo, encoding = c(0, 1, 2))
mono_indices <- which(maf < MAF)
if(length(mono_indices)!=0){
  GDmf = GDmo[,-mono_indices]
  GImf = GImo[-mono_indices,]
}
dim(GDmf)
dim(GImf)
}
#### Imputation
##### Imputation Using BEAGLE
if(Impuation=="None"){
  GDEre = GDmf
  GIEre = GImf
  dim(GDEre)
  dim(GIEre)

  fwrite(GDre,paste0(Study,"_GD_Filt_Imputed.csv"))
  fwrite(GIre[,1:3],paste0(Study,"_GI_Filt_Imputed.csv"))
  save(GDre,GIre,GT,file=paste0(Study,"_Filt_Imputed.RData"))
}else{

if(Impuation=="Beagle"){
  GImf$chr <- GImf$Chromosome
  GImf$rs <- GImf$SNP
  GImf$pos <- GImf$Position

  LD_file <- paste0(Study)
  vcf_LD <- numeric_2_VCF(GDmf, GImf)
  write_vcf(vcf_LD, outfile = LD_file)#Exports vcf
  # Assign parameters
  genotype_file = paste0(Study,".vcf")
  outfile = paste0(Study,"_imp")

  # Define a system command
  command1_prefix <- "java -Xmx1g -jar beagle.25Nov19.28d.jar"
  command_args <- paste(" gt=", genotype_file, " out=", outfile, sep = "")
  command1 <- paste(command1_prefix, command_args)
  test <- system(command1, intern = T)
  # Run BEAGLE using the system function, this will produce a gzip .vcf file
  Xvcf=read.vcfR(paste0(Study,"_imp",".vcf.gz"))
  Xbgl=data.frame(Xvcf@fix,Xvcf@gt)
  genotypes_imp <- VCF_2_numeric(Xbgl)[[1]]
  #Pre-imputed genotype matrix
  dim(GDmf)
  sum(is.na(GDmf))
  #Imputed genotype matrix
  dim(genotypes_imp)
  sum(is.na(genotypes_imp))

  #Remove MAF <5%
  maf <- calc_maf_apply(genotypes_imp, encoding = c(0, 1, 2))
  mono_indices <- which(maf < MAF)
  GDre = genotypes_imp[,-mono_indices]
  GIre = GImf[-mono_indices,]
  GIre=GIre[,-c("chr","rs","pos")]
  dim(GDre)
  dim(GIre)

  fwrite(GDre,paste0(Study,"_GD_Filt_Imputed.csv"))
  fwrite(GIre[,1:3],paste0(Study,"_GI_Filt_Imputed.csv"))
  save(GDre,GIre,GT,file=paste0(Study,"_Filt_Imputed.RData"))
  dim(GDre)
  dim(GIre)
  dim(GT)

}

if(Impuation=="KNN"){
  #Section for downloading package
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("impute")

  x=impute::impute.knn(as.matrix(t(GDEmf)))
  myGD_imp=t(x$data)
  myGD_imp<-round(myGD_imp,0)
  sum(is.na(myGD_imp))
  #Filter again
  #Remove MAF <5%
  maf <- calc_maf_apply(myGD_imp, encoding = c(0, 1, 2))
  mono_indices <- which(maf < MAF)
  GDEre = myGD_imp[,-mono_indices]
  GIEre = GIEmf[-mono_indices,]
  dim(GDEre)
  dim(GIEre)

  fwrite(GDre,paste0(Study,"_GD_Filt_Imputed.csv"))
  fwrite(GIre[,1:3],paste0(Study,"_GI_Filt_Imputed.csv"))
  save(GDre,GIre,GT,file=paste0(Study,"_Filt_Imputed.RData"))
}
}
##########PG################################
if(Method=="Two-Step"){
  if(!is.null(CV)){
    GBS.list=c()
    for(i in 1:length(Trial)){
      Pheno=Phenotype1 %>% filter(Env %in% c(Trial[i]))
      Pheno=Pheno[,c("Taxa","Env",Trait)]
      for(j in 1:length(Trait)){
        Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
        GBS<<-PGandCV(Pheno,GDre,GT,GIre,CV)
        mv(from = "GBS", to = paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),envir = globalenv())

        GBS.list=c(GBS.list,paste0("GBS_2_CV_",Trial[i],"_",Trait[j]))
        #save(list=paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_CV_",Trial[i],"_",Trait[j],".RData"))
      }

    }
    save(list = GBS.list, file=paste0("GBS_2_CV_",Study,".RData"))
  }else{
    GBS.list=c()
    for(i in 1:length(Trial)){
      Pheno=Phenotype1 %>% filter(Env %in% c(Trial[i]))
      Pheno=Pheno[,c("Taxa","Env",Trait)]
      for(j in 1:length(Trait)){
        Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
        GBS<<-PandG(Pheno,GDre,GT,GIre)
        mv(from = "GBS", to = paste0("GBS_2_",Trial[i],"_",Trait[j]),envir = globalenv())

        GBS.list=c(GBS.list,paste0("GBS_2_",Trial[i],"_",Trait[j]))
        #save(list=paste0("GBS_2_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_",Trial[i],"_",Trait[j],".RData"))
      }

    }
    save(list = try, file=paste0("GBS_2_",Study,".RData"))
  }
}

if(Method=="One-Step"){
  if(!is.null(CV)){
    Pheno=Phenotype1 %>% filter(Env %in% c(Trial))
    Pheno=Phenotype1[,c("Taxa","Env",Trait)]
    GBS<<-PGEandCV(Pheno,GDre,GT,GIre,CV)
    mv(from = "GBS", to = paste0("GBS_1_CV_",Study),envir = globalenv())
    save(list=paste0("GBS_1_CV_",Study),file=paste0("GBS_1_CV_",Study,".RData"))
  }else{
    Pheno=Phenotype1 %>% filter(Env %in% c(Trial))
    Pheno=Phenotype1[,c("Taxa","Env",Trait)]
    GBS<<-PGandE(Pheno,GDre,GT,GIre)
    mv(from = "GBS", to = paste0("GBS_1_",Study),envir = globalenv())
    save(list=paste0("GBS_1_",Study),file=paste0("GBS_1_",Study,".RData"))
  }
}


#############GS#########################
#Load files
if(Method=="Two-Step"){
  if(!is.null(CV)){
    load(file=paste0("GBS_2_CV_",Study,".RData"))
  }else{
    load(file=paste0("GBS_2_",Study,".RData"))
  }
}

if(Method=="One-Step"){
  if(!is.null(CV)){
    load(file=paste0("GBS_1_CV_",Study,".RData"))
  }else{
    load(file=paste0("GBS_1_",Study,".RData"))
  }
}

if(Method=="One-Step"){
  Matrix<<-GE_Matrix_IND(genotypes=get(paste0("GBS_1_",Study,"$geno")), phenotype=get(paste0("GBS_1_",Study,"$pheno")),trait=Trait,GE=GE,UN=UN,model=model)
  mv(from = "Matrix", to = paste0("Matrix_",Study,model),envir = globalenv())
  save(list=paste0("Matrix_",Study),file=paste0("GBS_1_",Study,".RData"))
}
}
if(GS==TRUE){
#####################GS######################################################
if(QC==FALSE){

  if(Method=="Two-Step"){
    if(!is.null(CV)){
      GBS.list=c()
      for(i in 1:length(Trial)){
        Pheno=Phenotype1 %>% filter(Env %in% c(Trial[i]))
        Pheno=Pheno[,c("Taxa","Env",Trait)]
        for(j in 1:length(Trait)){
          Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
          GBS<<-PGandCV(Pheno,GDre,GT,GIre,CV)
          mv(from = "GBS", to = paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),envir = globalenv())

          GBS.list=c(GBS.list,paste0("GBS_2_CV_",Trial[i],"_",Trait[j]))
          #save(list=paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_CV_",Trial[i],"_",Trait[j],".RData"))
        }

      }
      save(list = GBS.list, file=paste0("GBS_2_CV_",Study,".RData"))
    }else{
      GBS.list=c()
      for(i in 1:length(Trial)){
        Pheno=Phenotype1 %>% filter(Env %in% c(Trial[i]))
        Pheno=Pheno[,c("Taxa","Env",Trait)]
        for(j in 1:length(Trait)){
          Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
          GBS<<-PandG(Pheno,GDre,GT,GIre)
          mv(from = "GBS", to = paste0("GBS_2_",Trial[i],"_",Trait[j]),envir = globalenv())

          GBS.list=c(GBS.list,paste0("GBS_2_",Trial[i],"_",Trait[j]))
          #save(list=paste0("GBS_2_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_",Trial[i],"_",Trait[j],".RData"))
        }

      }
      save(list = try, file=paste0("GBS_2_",Study,".RData"))
    }
  }

  if(Method=="One-Step"){
    if(!is.null(CV)){
      Pheno=Phenotype1 %>% filter(Env %in% c(Trial))
      Pheno=Phenotype1[,c("Taxa","Env",Trait)]
      GBS<<-PGEandCV(Pheno,GDre,GT,GIre,CV)
      mv(from = "GBS", to = paste0("GBS_1_CV_",Study),envir = globalenv())
      save(list=paste0("GBS_1_CV_",Study),file=paste0("GBS_1_CV_",Study,".RData"))
    }else{
      Pheno=Phenotype1 %>% filter(Env %in% c(Trial))
      Pheno=Phenotype1[,c("Taxa","Env",Trait)]
      GBS<<-PGandE(Pheno,GDre,GT,GIre)
      mv(from = "GBS", to = paste0("GBS_1_",Study),envir = globalenv())
      save(list=paste0("GBS_1_",Study),file=paste0("GBS_1_",Study,".RData"))
    }
  }
}


if(Method=="Two-Step"){
  if(Outcome=="Tested"){
    if(Scheme=="K-Fold"){

      if(Package=="rrBLUP"){
        if(GAGS==TRUE){
          if(!is.null(PC)){
            #No CV
            Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  Kernel=Kernel,
                                                                                  PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                  CV=CV,
                                                                                  markers=markers,
                                                                                  folds = folds,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  transformation=transformation,
                                                                                  Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total
            )})
          }else{
            #No CV, No PC
            Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  Kernel=Kernel,
                                                                                  PCA=PC,
                                                                                  CV=CV,
                                                                                  markers=markers,
                                                                                  folds = folds,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  transformation=transformation,
                                                                                  Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total
            )})
          }
        }else{
          if(!is.null(CV)){
            if(!is.null(PC)){
              #No GAGS
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               Kernel=Kernel,
                                                                               PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }else{
              #No GAGS, No PC
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               Kernel=Kernel,
                                                                               PCA=PC,
                                                                               CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }
          }else{
            if(!is.null(PC)){
              #No GAGS, No CV
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               Kernel=Kernel,
                                                                               PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               CV=CV,
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }else{
              #No GAGS, No CV, No PC
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               Kernel=Kernel,
                                                                               PCA=PC,
                                                                               CV=CV,
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }
          }

        }


      }

      if(Package=="MAS"){
        if(GAGS==TRUE){
          if(!is.null(PC)){
            #No CV
            Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               CV=CV,
                                                                               folds = folds,
                                                                               transformation=transformation,
                                                                               Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                               GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                               GWAS=GWAS,
                                                                               alpha=alpha,
                                                                               threshold=threshold,
                                                                               QTN=QTN,
                                                                               PCA.total=PCA.total
            )})
          }else{
            #No CV, No PC
            Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               PCA=PC,
                                                                               CV=CV,
                                                                               folds = folds,
                                                                               transformation=transformation,
                                                                               Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                               GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                               GWAS=GWAS,
                                                                               alpha=alpha,
                                                                               threshold=threshold,
                                                                               QTN=QTN,
                                                                               PCA.total=PCA.total
            )})
          }
        }else{
          if(!is.null(CV)){
            if(!is.null(PC)){
              #No GAGS
              Results=sapply(1:Replications, function(i,...){Results=MAS_CV(phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                            PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                            CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                            folds = folds,
                                                                            transformation=transformation
              )})
            }else{
              #No GAGS, No PC
              Results=sapply(1:Replications, function(i,...){Results=MAS_CV(phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                            PCA=PC,
                                                                            CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                            folds = folds,
                                                                            transformation=transformation
              )})
            }
          }

        }


      }

      if(Package=="BGLR"){
        if(Type=="Ordinal"){
          if(GAGS==TRUE){
            if(!is.null(PC)){
              #No CV
              Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                          phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                          model=model,
                                                                                          Kernel=Kernel,
                                                                                          PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                          CV=CV,
                                                                                          markers=markers,
                                                                                          folds = folds,
                                                                                          nIter = nIter,
                                                                                          burnIn = burnIn,
                                                                                          Sparse=Sparse,
                                                                                          m=m,
                                                                                          degree=degree,
                                                                                          nL=nL,
                                                                                          Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                          GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                          GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                          GWAS=GWAS,
                                                                                          alpha=alpha,
                                                                                          threshold=threshold,
                                                                                          QTN=QTN,
                                                                                          PCA.total=PCA.total
              )})
            }else{
              #No CV, No PC
              Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                          phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                          model=model,
                                                                                          Kernel=Kernel,
                                                                                          PCA=PC,
                                                                                          CV=CV,
                                                                                          markers=markers,
                                                                                          folds = folds,
                                                                                          nIter = nIter,
                                                                                          burnIn = burnIn,
                                                                                          Sparse=Sparse,
                                                                                          m=m,
                                                                                          degree=degree,
                                                                                          nL=nL,
                                                                                          Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                          GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                          GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                          GWAS=GWAS,
                                                                                          alpha=alpha,
                                                                                          threshold=threshold,
                                                                                          QTN=QTN,
                                                                                          PCA.total=PCA.total
              )})
            }
          }else{
            if(!is.null(CV)){
              if(!is.null(PC)){
                #No GAGS
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                       CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                                       markers=markers,
                                                                                       folds = folds,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }else{
                #No GAGS, No PC
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       PCA=PC,
                                                                                       CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                                       markers=markers,
                                                                                       folds = folds,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }
            }else{
              if(!is.null(PC)){
                #No GAGS, No CV
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                       CV=CV,
                                                                                       markers=markers,
                                                                                       folds = folds,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }else{
                #No GAGS, No CV, No PC
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       PCA=PC,
                                                                                       CV=CV,
                                                                                       markers=markers,
                                                                                       folds = folds,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }
            }

          }
        }else{

          if(GAGS==TRUE){
            if(!is.null(PC)){
              #No CV
              Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  model=model,
                                                                                  Kernel=Kernel,
                                                                                  PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                  CV=CV,
                                                                                  markers=markers,
                                                                                  folds = folds,
                                                                                  nIter = nIter,
                                                                                  burnIn = burnIn,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  transformation=transformation,
                                                                                  Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total
              )})
            }else{
              #No CV, No PC
              Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  model=model,
                                                                                  Kernel=Kernel,
                                                                                  PCA=PC,
                                                                                  CV=CV,
                                                                                  markers=markers,
                                                                                  folds = folds,
                                                                                  nIter = nIter,
                                                                                  burnIn = burnIn,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  transformation=transformation,
                                                                                  Y=get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total
              )})
            }
          }else{
            if(!is.null(CV)){
              if(!is.null(PC)){
                #No GAGS
                Results=sapply(1:Replications, function(i,...){Results=BLGR_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }else{
                #No GAGS, No PC
                Results=sapply(1:Replications, function(i,...){Results=BLGR_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               PCA=PC,
                                                                               CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }
            }else{
              if(!is.null(PC)){
                #No GAGS, No CV
                Results=sapply(1:Replications, function(i,...){Results=BLGR_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               CV=CV,
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }else{
                #No GAGS, No CV, No PC
                Results=sapply(1:Replications, function(i,...){Results=BLGR_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               PCA=PC,
                                                                               CV=CV,
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }
            }

          }

        }
      }

      if(Package=="GAPIT"){
        if(!is.null(CV)){
          #No GAGS
          Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                             phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                             myGD=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                             model=model,
                                                                             PCA.total=PC,
                                                                             CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                             kinship=kinship,
                                                                             markers=markers,
                                                                             folds = folds,
                                                                             transformation=transformation
          )})

        }else{
          #No GAGS, No CV
          Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                             phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                             myGD=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                             model=model,
                                                                             PCA.total=PC,
                                                                             CV=CV,
                                                                             kinship=kinship,
                                                                             markers=markers,
                                                                             folds = folds,
                                                                             transformation=transformation
          )})

        }
      }

      if(Package=="caret"){
        Results=sapply(1:Replications, function(i,...){Results=Caret_Models_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               type=type,
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               folds = folds,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation,
                                                                               sampling=sampling,
                                                                               repeats=repeats,
                                                                               method=method
        )})
      }

      if(Package=="GLM"){
        Results=sapply(1:Replications, function(i,...){Results=GLM_CV(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                      phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                      fam=fam,
                                                                      Kernel=Kernel,
                                                                      markers=markers,
                                                                      folds = folds,
                                                                      Sparse=Sparse,
                                                                      m=m,
                                                                      degree=degree,
                                                                      nL=nL
        )})
      }

    }
    if(Scheme=="VS"){
      if(Package=="rrBLUP"){
        if(GAGS==TRUE){
          if(!is.null(PC)){
            #No CV
            Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                  test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                  test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                  test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                                  Kernel=Kernel,
                                                                                  markers=markers,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total,
                                                                                  transformation=transformation
            )})

          }else{
            #No CV, No PC
            Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  train_PCA=PC,
                                                                                  test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                  test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                  test_PCA=PC,
                                                                                  Kernel=Kernel,
                                                                                  markers=markers,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total,
                                                                                  transformation=transformation
            )})
          }
        }else{
          if(!is.null(CV)){
            if(!is.null(PC)){
              #No GAGS
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                               test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }else{
              #No GAGS, No PC
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=PC,
                                                                               train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=PC,
                                                                               test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }
          }else{
            if(!is.null(PC)){
              #No GAGS, No CV
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               train_CV=CV,
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                               test_CV=CV,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }else{
              #No GAGS, No CV, No PC
              Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=PC,
                                                                               train_CV=CV,
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=PC,
                                                                               test_CV=CV,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
              )})
            }
          }

        }


      }

      if(Package=="MAS"){
        if(GAGS==TRUE){
          if(!is.null(PC)){
            #No CV
            Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                               train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                               train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                               markers=markers,
                                                                               GWAS=GWAS,
                                                                               alpha=alpha,
                                                                               threshold=threshold,
                                                                               QTN=QTN,
                                                                               PCA.total=PCA.total,
                                                                               transformation=transformation
            )})
          }else{
            #No CV, No PC
            Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                               train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                               train_PCA=PC,
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=PC,
                                                                               markers=markers,
                                                                               GWAS=GWAS,
                                                                               alpha=alpha,
                                                                               threshold=threshold,
                                                                               QTN=QTN,
                                                                               PCA.total=PCA.total,
                                                                               transformation=transformation
            )})
          }
        }else{
          if(!is.null(CV)){
            if(!is.null(PC)){
              #No GAGS
              Results=sapply(1:Replications, function(i,...){Results=MAS_VS(train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                            train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                            train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                            test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                            test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                            test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                            transformation=transformation
              )})
            }else{
              #No GAGS, No PC
              Results=sapply(1:Replications, function(i,...){Results=MAS_VS(train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                            train_PCA=PC,
                                                                            train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                            test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                            test_PCA=PC,
                                                                            test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                            transformation=transformation
              )})
            }
          }

        }


      }
      if(Package=="BGLR"){
        if(Type=="Ordinal"){
          if(GAGS==TRUE){
            if(!is.null(PC)){
              #No CV
              Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                          train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                          train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                          train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                          train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                          test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                          test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                          test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                                          model=model,
                                                                                          Kernel=Kernel,
                                                                                          markers=markers,
                                                                                          nIter = nIter,
                                                                                          burnIn = burnIn,
                                                                                          Sparse=Sparse,
                                                                                          m=m,
                                                                                          degree=degree,
                                                                                          nL=nL,
                                                                                          GWAS=GWAS,
                                                                                          alpha=alpha,
                                                                                          threshold=threshold,
                                                                                          QTN=QTN,
                                                                                          PCA.total=PCA.total
              )})
            }else{
              #No CV, No PC
              Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                          train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                          train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                          train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                          train_PCA=PC,
                                                                                          test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                          test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                          test_PCA=PC,
                                                                                          model=model,
                                                                                          Kernel=Kernel,
                                                                                          markers=markers,
                                                                                          nIter = nIter,
                                                                                          burnIn = burnIn,
                                                                                          Sparse=Sparse,
                                                                                          m=m,
                                                                                          degree=degree,
                                                                                          nL=nL,
                                                                                          GWAS=GWAS,
                                                                                          alpha=alpha,
                                                                                          threshold=threshold,
                                                                                          QTN=QTN,
                                                                                          PCA.total=PCA.total
              )})
            }
          }else{
            if(!is.null(CV)){
              if(!is.null(PC)){
                #No GAGS
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                       train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                                       test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                       test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                       test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                                       test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       markers=markers,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }else{
                #No GAGS, No PC
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       train_PCA=PC,
                                                                                       train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                                       test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                       test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                       test_PCA=PC,
                                                                                       test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       markers=markers,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }
            }else{
              if(!is.null(PC)){
                #No GAGS, No CV
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                       train_CV=CV,
                                                                                       test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                       test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                       test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                                       test_CV=CV,
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       markers=markers,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }else{
                #No GAGS, No CV, No PC
                Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                       train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                       train_PCA=PC,
                                                                                       train_CV=CV,
                                                                                       test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                       test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                       test_PCA=PC,
                                                                                       test_CV=CV,
                                                                                       model=model,
                                                                                       Kernel=Kernel,
                                                                                       markers=markers,
                                                                                       nIter = nIter,
                                                                                       burnIn = burnIn,
                                                                                       Sparse=Sparse,
                                                                                       m=m,
                                                                                       degree=degree,
                                                                                       nL=nL
                )})
              }
            }

          }
        }else{

          if(GAGS==TRUE){
            if(!is.null(PC)){
              #No CV
              Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                                  test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                  test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                  test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                                  model=model,
                                                                                  Kernel=Kernel,
                                                                                  markers=markers,
                                                                                  nIter = nIter,
                                                                                  burnIn = burnIn,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  transformation=transformation,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total
              )})
            }else{
              #No CV, No PC
              Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                                  train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                                  train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                                  train_GD=get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                                  train_PCA=PC,
                                                                                  test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                                  test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                                  test_PCA=PC,
                                                                                  model=model,
                                                                                  Kernel=Kernel,
                                                                                  markers=markers,
                                                                                  nIter = nIter,
                                                                                  burnIn = burnIn,
                                                                                  Sparse=Sparse,
                                                                                  m=m,
                                                                                  degree=degree,
                                                                                  nL=nL,
                                                                                  transformation=transformation,
                                                                                  GWAS=GWAS,
                                                                                  alpha=alpha,
                                                                                  threshold=threshold,
                                                                                  QTN=QTN,
                                                                                  PCA.total=PCA.total
              )})
            }
          }else{
            if(!is.null(CV)){
              if(!is.null(PC)){
                #No GAGS
                Results=sapply(1:Replications, function(i,...){Results=BLGR_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                               test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }else{
                #No GAGS, No PC
                Results=sapply(1:Replications, function(i,...){Results=BLGR_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=PC,
                                                                               train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=PC,
                                                                               test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }
            }else{
              if(!is.null(PC)){
                #No GAGS, No CV
                Results=sapply(1:Replications, function(i,...){Results=BLGR_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=get(paste0("GBS_2_",Training,"_",Trait))$PC[,1:PC],
                                                                               train_CV=CV,
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=get(paste0("GBS_2_",Prediction,"_",Trait))$PC[,1:PC],
                                                                               test_CV=CV,
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }else{
                #No GAGS, No CV, No PC
                Results=sapply(1:Replications, function(i,...){Results=BLGR_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               train_PCA=PC,
                                                                               train_CV=CV,
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               test_PCA=PC,
                                                                               test_CV=CV,
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               nIter = nIter,
                                                                               burnIn = burnIn,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation
                )})
              }
            }

          }

        }
      }
      if(Package=="caret"){
        Results=sapply(1:Replications, function(i,...){Results=Caret_Models_CV(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                               train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                               test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                               test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                               type=type,
                                                                               model=model,
                                                                               Kernel=Kernel,
                                                                               markers=markers,
                                                                               Sparse=Sparse,
                                                                               m=m,
                                                                               degree=degree,
                                                                               nL=nL,
                                                                               transformation=transformation,
                                                                               sampling="up",
                                                                               repeats=5,
                                                                               method="repeatedcv"
        )})
      }


      if(Package=="GAPIT"){
        if(!is.null(CV)){
          #No GAGS
          Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_VS(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                             phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                             train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                             train_CV=get(paste0("GBS_2_",Training,"_",Trait))$CV[,-1],
                                                                             test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$numeric,
                                                                             test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                             test_GM=get(paste0("GBS_2_",Prediction,"_",Trait))$map,
                                                                             test_CV=get(paste0("GBS_2_",Prediction,"_",Trait))$CV[,-1],
                                                                             model=model,
                                                                             PCA.total=PC,
                                                                             kinship=kinship,
                                                                             markers=markers,
                                                                             transformation=transformation
          )})

        }else{
          #No GAGS, No CV
          Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_VS(genotypes = get(paste0("GBS_2_",Training,"_",Trait))$numeric,
                                                                             phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                             train_GM=get(paste0("GBS_2_",Training,"_",Trait))$map,
                                                                             train_CV=CV,
                                                                             test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$numeric,
                                                                             test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                             test_GM=get(paste0("GBS_2_",Prediction,"_",Trait))$map,
                                                                             test_CV=CV,
                                                                             model=model,
                                                                             PCA.total=PC,
                                                                             kinship=kinship,
                                                                             markers=markers,
                                                                             transformation=transformation
          )})

        }
      }
      if(Package=="GLM"){
        Results=sapply(1:Replications, function(i,...){Results=GLM_VS(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                      train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                      test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                      test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                      fam=fam,
                                                                      Kernel=Kernel,
                                                                      markers=markers,
                                                                      Sparse=Sparse,
                                                                      m=m,
                                                                      degree=degree,
                                                                      nL=nL
        )})
      }

    }
  }
  if(Outcome=="Untested"){
    if(Scheme=="K-Fold"){
      #Is there even a point for K-fold for untested?
#We can do this if we include blank values, but not in rrBLUP or DL
    if(package=="rrBLUP"){
      if(Kernel=="Markers"){
        NULL
      }else{}
    }
    }
    if(Scheme=="VS"){

    }
  }
}


if(Method=="One-Step"){
  #We can do this with rrBLUP if we use Kernel/kin.blup
  if(Outcome=="Tested"){
    if(Scheme=="K-Fold"){

    }
    if(Scheme=="VS"){

    }
  }
  if(Outcome=="Untested"){
    if(Scheme=="K-Fold"){
      #We can do this if we include blank values, but not in rrBLUP or DL
      if(package=="rrBLUP"){
        if(Kernel=="Markers"){
          NULL
        }else{}
    }
    if(Scheme=="VS"){

    }
  }
}



if(!is.null(CV)){
  CV_Message=TRUE
}else{
  CV_Message=NULL
}

if(!is.null(PCA)){
  PCA_Message=TRUE
}else{
  PCA_Message=NULL
}

if(GAGS==TRUE){
  GAGS_Message=TRUE
}else{
  GAGS_Message=NULL
}

if(Outcome=="Tested"){
Results_Accuracy=Extract_ACC(Results,Replications,Training,Model,Kernel,CV_Message,Trait)

Results_Predictions=Extract_Pred(Results,Replications,Training,Model,Kernel,CV_Message,Trait)

Results_All=list(Accuracy=Results_Accuracy,Predictions=Results_Predictions)
}

return(Results_All)
}
}

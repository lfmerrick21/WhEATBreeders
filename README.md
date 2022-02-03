User Manual and Tutorial for Genomic Selection in WhEATBreeders
================
Lance Merrick
February 2, 2022

<style type="text/css">

h1.title {
  font-size: 38px;
  color: black;
  text-align: center;
}
h4.author { /* Header 4 - and the author and data headers use this too  */
    font-size: 25px;
  font-family: "Times New Roman", Times, serif;
  color: black;
  text-align: center;
}
h4.date { /* Header 4 - and the author and data headers use this too  */
  font-size: 18px;
  font-family: "Times New Roman", Times, serif;
  color: black;
  text-align: center;
}
</style>

<p>

 

</p>

<p>

 

</p>

<!-- test inserting image-->

<center>

![](wheatclipart.png)

</center>

<p>

 

</p>

<center>

<font size="6"> **Wrappers for Easy Association, Tools, and Selection
for Breeders (WhEATBreeders)** </font>

</center>

<p>

 

</p>

<p>

 

</p>

### **Table of contents**

#### Introduction

#### Description of package functions

#### Getting started

#### Required Inputes

#### Optional Inputes

#### Outputs and Examples

#### Examples with know QTNs

#### Rapid run

#### Rapid run options

#### Frequently asked questions

#### Further Information

<p>

 

</p>

<p>

 

</p>

<p>

 

</p>

### **Introduction**

A common problem associated with performing Genome Wide Association
Studies (GWAS) is the abundance of false positive and false negative
associated with population structure. One way of dealing with this issue
is to first calculate Principal Components (PC) and then use those PCs
as to account for population structure when running the GWAS. This
method has shown to not only reduce false positives but also increase
the power of the test. Here we present an R based package that can
preform GWAS using PCA and user imputed covariate to calculate SNPs
associated with a phenotype called “GWhEAT”. It will also check to see
if the PC are in linear dependence to the covariates and remove the PCs
if they are. Using the freely availably R studio software this packaged
is designed to quickly and efficiently calculate a P-value for every SNP
phenotype association. The results from this package are including but
not limited to a Manhattan plot, QQ-plot and table with every recorded
p-value.

### **Description of package functions**

For a full list of functions within GWhEAT see “Reference\_Manual.pdf”
this file contains not only the full list of functions but also a
description of each. The pdf also has each functions arguments listed.
And like with all R packages once GWhEAT is installed and loaded you can
type ?function\_name and that specific function’s full descriptions will
appear in the help tab on RStudio.

### **Getting started**

First if you do not already have R and R studio installed on your
computer head over to <https://www.r-project.org/> and install the
version appropriate for you machine. Once R and R studio are installed
you will need to install the GWhEAT package since this is a working
package in it’s early stages of development it’s only available through
Github. To download files off Github first download and load the library
of the packaged “devtools” using the code below.

\#Read in GBS data

``` r
#file.choose()
rm(list = ls(all.names = TRUE))
gc()
library(data.table)
library(dplyr)
#library(GAPIT3)
library(bigmemory)
library(biganalytics)
library(ggplot2)
library(dplyr)
library(rrBLUP)
library(knitr)
library(cvTools)
library(vcfR)
require(compiler) #for cmpfun
#source("http://zzlab.net/GAPIT/GAPIT.library.R")
#source("http://zzlab.net/GAPIT/gapit_functions.txt") #make sure compiler is running
#source("http://zzlab.net/GAPIT/emma.txt")

#Better for FDR function
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

source('IMPUTATION_functions.R')
source('FUNCTION_numeric_2_VCF.R')
source('FUNCTION_VCF_2_numeric.R')
source('Impute_functions.R')
source("FUNCTIONS_GS_COMP.R")
source("FUNCTIONS_GS_RVC.R")
#memory.size()
#gc()
gen1 <- fread("WAC_2016-2020_production_filt.hmp.txt",fill=TRUE)
Pheno
```

\#Read in 2020 Lind Taxa

``` r
#file.choose()
#N x 1 with just taxa list
gname=read.csv("Lind_2020_taxa.csv",header=T)
#str(gname)
names(gname)<-"taxa"
gname$taxa=as.character(gname$taxa)
gname=as.vector(gname$taxa)
str(gname)
#gname=gname[-1]
#str(gname)
```

\#Filter GBS for Lind 2020

``` r
names.use <- names(gen1)[(names(gen1) %in% gname)]
output <- gen1[, names.use, with = FALSE]
output=cbind(gen1[,1:11],output)
#str(output)
output[output=="NA"]<-"NA"
output$`assembly#`=as.character(output$`assembly#`)
output$center=as.character(output$center)
output$protLSID=as.character(output$protLSID)
output$assayLSID=as.character(output$assayLSID)
output$panelLSID=as.character(output$panelLSID)
output$QCcode=as.character(output$QCcode)
str(output)
fwrite(output,"WA16_20_LIND_20_filt.hmp.txt",sep="\t",row.names = FALSE,col.names = TRUE)
save(output,file="WA16_20_LIND_20_filt.hmp.txt.RData")
#readLines(output)
```

\#Convert to Hapmap

``` r
getwd()
#Only needed if you're reloading
output=read_hapmap("WA16_20_LIND_20_filt.hmp.txt")

save(output,file="LIND_20_Hapmap_Raw.RData")
load(file="WA16_20_LIND_20_filt.hmp.txt.RData")
#Convert to Numeric Output
#Impute is none because we will impute using BEAGLE
outG=GAPIT.HapMap(output,SNP.impute="None")#None could be "Middle"
#outG$GD[1:10,1:10]
GTL=outG$GT
GTL=cbind(rownames(GTL),GTL)
fwrite(GTL,"GT_L20_taxa.csv")
GDL=outG$GD
GIL=outG$GI

rownames(GDL)=rownames(GTL)
colnames(GDL)=GIL$SNP
fwrite(GDL,"GD_L20_numeric_Raw.csv")
fwrite(GIL,"GI_L20_map_Raw.csv")
save(GDL,GIL,GTL,file="LIND_20_GD_GI_Raw.RData")
load(file="LIND_20_GD_GI_Raw.RData")
GTL=read.csv(file="GT_L20_taxa.csv")
GTL
```

\#Remove Missing Data \>80%, MAF \<5% and Monomorphic Markers

``` r
GDL[1:10,1:10]
dim(GDL)
#Remove missing data with more than 20%
mr=calc_missrate(GDL)
mr_indices <- which(mr > 0.20)
GDLmr = GDL[,-mr_indices]
GILmr = GIL[-mr_indices,]

dim(GDLmr)
dim(GILmr)

#Remove Monomorphic Markers
maf <- calc_maf_apply(GDLmr, encoding = c(0, 1, 2))
mono_indices <- which(maf ==0)
GDLmo = GDLmr[,-mono_indices]
GILmo = GILmr[-mono_indices,]

dim(GDLmo)
dim(GILmo)

#Remove Markers with MAF less than 5%
maf <- calc_maf_apply(GDLmo, encoding = c(0, 1, 2))
mono_indices <- which(maf <0.05)
GDLmf = GDLmo[,-mono_indices]
GILmf = GILmo[-mono_indices,]

dim(GDLmf)
dim(GILmf)
GDLmf[1:10,1:10]
fwrite(GDLmf,"GD_L20_numeric_Filtered.csv")
fwrite(GILmf,"GI_L20_map_Filtered.csv")
save(GDLmf,GILmf,GTL,file="LIND_20_GD_GI_Filtered.RData")
```

\#Imputation Using BEAGLE

``` r
GDLmf[1:10,1:10]
  GILmf$chr <- GILmf$Chromosome
  GILmf$rs <- GILmf$SNP
  GILmf$pos <- GILmf$Position

LD_file <- "LD_Numeric"
vcf_LD <- numeric_2_VCF(GDLmf, GILmf)
write_vcf(vcf_LD, outfile = LD_file)#Exports vcf 
# Assign parameters
  genotype_file = "LD_Numeric.vcf" 
  outfile = "LD_Numeric_imp"
  
  # Define a system command
  command1_prefix <- "java -Xmx1g -jar beagle.25Nov19.28d.jar"
  command_args <- paste(" gt=", genotype_file, " out=", outfile, sep = "")
  command1 <- paste(command1_prefix, command_args)
  test <- system(command1, intern = T)
  # Run BEAGLE using the system function, this will produce a gzip .vcf file
  Xvcf=read.vcfR("LD_Numeric_imp.vcf.gz")
  Xbgl=data.frame(Xvcf@fix,Xvcf@gt)
  genotypes_imp <- VCF_2_numeric(Xbgl)[[1]]
  #Pre-imputed genotype matrix
  dim(GDLmf)
  sum(is.na(GDLmf))
  #Imputed genotype matrix
  dim(genotypes_imp)
  sum(is.na(genotypes_imp))
  
genotypes_imp[1:10,1:10]
#Remove MAF <5%
maf <- calc_maf_apply(genotypes_imp, encoding = c(0, 1, 2))
mono_indices <- which(maf <0.05)
GDLre = genotypes_imp[,-mono_indices]
GILre = GILmf[-mono_indices,]
dim(GDLre)
#dim(GIEmo)
dim(GDLmf)
dim(genotypes_imp)
GILre
#If Beagle didn't work could us KNN
#x=impute::impute.knn(as.matrix(t(GDEmf)))
#myGD_imp=t(x$data)
#myGD_imp<-round(myGD_imp,0)
#sum(is.na(myGD_imp))
#View(myGD_imp)

fwrite(GDLre,"GD_L20_numeric_Filt_Imputed.csv")
fwrite(GILre[,1:3],"GI_L20_map_Filt_Imputed.csv")
save(GDLre,GILre,GTL,file="LIND_20_GD_GI_Filt_Imputed.RData")
dim(GDLre)
dim(GILre)
dim(GTL)
```

\#Function to get data in order as above chunck but integrated into a
function to also combine all files for a nice list with PCs
included.

``` r
PandG <- function(Pheno_Em_bl=NULL,myGD_EM=NULL,mytaxa_EM=NULL,myGM_EM=NULL){
  colnames(Pheno_Em_bl)[1]<-c("Genotype")
  Pheno_Em_bl=data.frame(Pheno_Em_bl)
  rownames(Pheno_Em_bl) <- Pheno_Em_bl$Genotype
  
  Markers_Em_bl<-myGD_EM
  
  rownames(Markers_Em_bl) <- mytaxa_EM$V1
  colnames(Markers_Em_bl) <- myGM_EM$SNP
  
  
  Pheno_Em_bl <- Pheno_Em_bl[rownames(Pheno_Em_bl) %in% rownames(Markers_Em_bl),]
  Markers_Em_bl <- Markers_Em_bl[rownames(Markers_Em_bl) %in% rownames(Pheno_Em_bl),]
  
  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]
  
  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]
  
  myGM<-myGM_EM[myGM_EM$SNP %in% colnames(Markers_Em_bl),]
  myGD<-data.frame(rownames(Markers_Em_bl),Markers_Em_bl)
  colnames(myGD)[1]<-c("taxa")
  
  PCA_bl_em=prcomp(Markers_Em_bl)
  myPCA_bl_em=PCA_bl_em$x
  PGPCA=list(pheno=Pheno_Em_bl,geno=Markers_Em_bl,map=myGM,numeric=myGD,PC=myPCA_bl_em)
  return(PGPCA)
}

#If you also want to integrate covariates into the list
PGandCV <- function(Pheno_Em_bl=NULL,myGD_EM=NULL,mytaxa_EM=NULL,myGM_EM=NULL,CV=NULL){
  colnames(Pheno_Em_bl)[1]<-c("Genotype")
  rownames(Pheno_Em_bl) <- Pheno_Em_bl$Genotype
  
  Markers_Em_bl<-myGD_EM
  
  rownames(Markers_Em_bl) <- mytaxa_EM$V1
  colnames(Markers_Em_bl) <- myGM_EM$SNP
  
  colnames(CV)[1]<-c("Genotype")
  rownames(CV) <- CV$Genotype
  library(tidyr)
  library(Hmisc)
  for(i in 2:ncol(CV)){
    CV[,i]=impute(CV[,i])
    CV[,i]=as.numeric(CV[,i])
  }
  
  
  Pheno_Em_bl <- Pheno_Em_bl[rownames(Pheno_Em_bl) %in% rownames(Markers_Em_bl),]
  Markers_Em_bl <- Markers_Em_bl[rownames(Markers_Em_bl) %in% rownames(Pheno_Em_bl),]
  
  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]
  
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
```

\#\#\#\#\#\#START HERE IF GENOTYPE DATA IS
DONE\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
\#This gets into my own files This section gets your data in order but
using a function

``` r
#load phenotypic and genotypic data in.
#setwd("F:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/Documents/Genomic Selection/Genomic Selection Pipeline/GWAS Pipeline")
load(file="LIND_20_GD_GI_Filt_Imputed.RData")
#Without Covariates
#Remove missing data
Pheno_LND_20_1<-Pheno[complete.cases(Pheno[,2]),]
GBS_LND_20_1=PandG(Pheno_LND_20_1,GDLre,GTL,GILre)

Pheno_LND_20_2<-Pheno[complete.cases(Pheno[,3]),]
GBS_LND_20_2=PandG(Pheno_LND_20_2,GDLre,GTL,GILre)

#If you have covariates
#Remove missing data
Pheno_LND_20<-Pheno[complete.cases(Pheno[,2]),]
GBS_LND_20_CV=PGandCV(Pheno_LND_20,GDLre,GTL,GILre,myCV_EM)

save(GBS_LND_20,GBS_LND_20_CV,file="LND_20_GBS.RData")
```

\#GS without CV

``` r
load(file="LND_20_GBS.RData")
LND_20=test_all_models_BLUP_pc_mean_recode(genotypes = GBS_LND_20$geno, 
                                    phenotype = GBS_LND_20$pheno[,c(1,2)],
                                    PCA=GBS_LND_20$PC[,1:3])
#Replications
LND_20_50=Sapply(1:50, function(i,...){LND_20_50=test_all_models_BLUP_pc_mean_recode(genotypes = GBS_LND_20$geno, 
                                    phenotype = GBS_LND_20$pheno[,c(1,2)],
                                    PCA=GBS_LND_20$PC[,1:3])})



Extract_ACC=function(results,nrep,pop,year,loc,model,trait,pheno){
model_vect <- c("MAE","MAE_PC","Pearson","Pearson_PC","R2","R2_PC","RMSE","RMSE_PC","Spearman" ,"Spearman_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(results[1,])){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(results[1,i]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
colnames(BGLR_acc_table)<-c("Rep","Metric","Accuracy")
po=rep(pop,nrep)
yr=rep(year,nrep)
loc=rep(loc,nrep)
mo=rep(model,nrep)
tr=rep(trait,nrep)
ph=rep(pheno,nrep)
ac=data.frame(Population=po,Year=yr,Location=loc,Prediction=mo,Trait=tr,Adjustment=ph,BGLR_acc_table)
return(ac)
}

SVM_BL17_28=Extract_ACC(LND_20_50,50,"Breeding Lines","2017","Lind","SVM","IT","Factor")
RVC_ALL_FIL=rbind(SVM_BL17_28)
```

\#REG BL
17

``` r
BL_REG_17_acc=RVC_ALL_FIL %>% filter(Population %in% c("Breeding Lines") & Location %in% c("Lind") & Year %in% c("2017") & Metric %in% c("Pearson")& Prediction %in% c("rrBLUP","GLM","SVMR"))

Mean_BL_REG_17_acc=BL_REG_17_acc %>%
  group_by(Specific,Metric) %>%
  summarise(mean = mean(Accuracy))


BL_REG_17_acc_IT=BL_REG_17_acc %>% filter(Trait %in% c("IT"))


Tukey_test_BL_REG_17_acc_IT <- aov(Accuracy~Specific, data=BL_REG_17_acc_IT) %>%
HSD.test("Specific", group=TRUE) %>%
.$groups %>%
as_tibble(rownames="Specific") %>%
rename("Letters"="groups")

Tukey_test_BL_REG_17_std_IT<-aov(Accuracy~Specific, data=BL_REG_17_acc_IT) %>%
HSD.test("Specific", group=TRUE)%>%.$means%>%
as_tibble(rownames="Specific")%>%.[,c("Specific","std")]


ggviolin(BL_REG_17_acc_IT, x = "Specific", y = "Accuracy", fill = "Specific",
         add = "boxplot", add.params = list(fill = "white"),legend="none",title="Infection Type")+ geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(colour="black",face="bold" ,size=8, angle = 45, vjust=0.5))+
    geom_label(data=Tukey_test_BL_REG_17_acc_IT, aes(label=Letters))
```

\#GS with CV

``` r
load(file="LND_20_GBS.RData")
LND_20_CV=test_all_models_BLUP_pc_mean_recode(genotypes = GBS_LND_20_CV$geno, 
                                    phenotype = GBS_LND_20_CV$pheno[,c(1,2)],
                                    PCA=GBS_LND_20_CV$PC[,1:3],
                                    CV=GBS_LND_20_CV$CV[,2])
#Replicatoins
LND_20_CV_50=Sapply(1:50, function(i,...){LND_20_CV=test_all_models_BLUP_pc_mean_recode(genotypes = GBS_LND_20_CV$geno, 
                                    phenotype = GBS_LND_20_CV$pheno[,c(1,2)],
                                    PCA=GBS_LND_20_CV$PC[,1:3],
                                    CV=GBS_LND_20_CV$CV[,2])})
```

\#Without CV

``` r
LND_20=Caret_Models_Mean_Matrix(genotypes=GBS_LND_20$geno,
                                phenotype= GBS_LND_20$pheno[,c(1,2)],
                                model="svmRadial",
                                Matrix="VanRaden",
                                sampling="up")#Unbalanced
```

\#With
CV

``` r
LND_20_CV=Caret_Models_Mean_Matrix(genotypes=cbind(GBS_LND_20_CV$CV[,c(2:5)],GBS_LND_20_CV$geno),
                                   phenotype= GBS_LND_20_CV$pheno[,c(1,2)],
                                   model="svmRadial",
                                   Matrix="VanRaden",
                                   sampling="up")
```

\#Validation GS \#Need more Info to show you \#Without
CV

``` r
test_all_models_BLUP_vs_pc_mean_recode(train_genotypes = GBS_BQALL_comp2$geno, 
                                       train_phenotype = GBS_BQALL_comp2$pheno[,c(1,2)],
                                       train_PCA=GBS_BQALL_comp2$PC[,1:3],
                                       test_genotypes = GBS_BQALL_comp4$geno,
                                       test_phenotype = GBS_BQALL_comp4$pheno[,c(1,4)],
                                       test_PCA=GBS_BQALL_comp4$PC[,1:3])
```

\#With
CV

``` r
test_all_models_BLUP_vs_pc_mean_recode(train_genotypes = GBS_BQALL_comp2$geno, 
                                       train_phenotype = GBS_BQALL_comp2$pheno[,c(1,2)],
                                       train_PCA=GBS_BQALL_comp2$PC[,1:3],
                                       train_CV=GBS_BQALL_comp2$CV[,2],
                                       test_genotypes = GBS_BQALL_comp4$geno,
                                       test_phenotype = GBS_BQALL_comp4$pheno[,c(1,4)],
                                       test_PCA=GBS_BQALL_comp4$PC[,1:3],
                                       test_CV=GBS_BQALL_comp4$CV[,2])
```

\#R Script for Kamiak

``` r
  source("FUNCTIONS_GS_COMP.R")
  load(file = "LND_20_GBS.RData")
  library(dplyr)
  library(rrBLUP)
  library(BGLR)
  library(tidyr)
  library(caret)
  library(Metrics)
  library(mpath)

library(parallel)
cores=10
cl <- makeForkCluster(cores)
clusterSetRNGStream(cl)

LND_20=parSapply(cl, 1:50, function(i,...){LND_20=test_all_models_BLUP_pc_mean_recode(genotypes = GBS_LND_20$geno, 
                                    phenotype = GBS_LND_20$pheno[,c(1,2)],
                                    PCA=GBS_LND_20$PC[,1:3])
save(LND_20,file=paste0("LND_20_",i,".RData"))
LND_20})

stopCluster(cl)
save(LND_20,file = "LND_20_Jason.RData")
```

\#R Script for Kamiak with CV

``` r
  source("FUNCTIONS_GS_COMP.R")
  load(file = "LND_20_GBS.RData")
  library(dplyr)
  library(rrBLUP)
  library(BGLR)
  library(tidyr)
  library(caret)
  library(Metrics)
  library(mpath)

library(parallel)
cores=10
cl <- makeForkCluster(cores)
clusterSetRNGStream(cl)

LND_20_CV=parSapply(cl, 1:50, function(i,...){LND_20_CV=test_all_models_BLUP_pc_mean_recode(genotypes = GBS_LND_20_CV$geno, 
                                    phenotype = GBS_LND_20_CV$pheno[,c(1,2)],
                                    PCA=GBS_LND_20_CV$PC[,1:3],
                                    CV=GBS_LND_20_CV$CV[,2])
save(LND_20_CV,file=paste0("LND_20_CV_",i,".RData"))
LND_20_CV})

stopCluster(cl)
save(LND_20_CV,file = "LND_20_CV_Jason.RData")
```

### **Frequently asked questions**

1.  Where can I get GWhEAT? GWhEAT is currently only available on github
    via Samuel Prather’s public repository found at
    <https://github.com/stp4freedom/GWhEAT>.  
      
2.  How to download GWhEAT? Please refer the “Getting started” section
    above where I show how to download and install GWhEAT using
    devtools.  
      
3.  What SNP format does GWhEAT take? Currently GWhEAT only accepts SNP
    data in numerical form.  
      
4.  What are the required inputs to run GWhEAT? You will need SNP data
    in numerical form, phenotypic data and a genetic map that has SNP
    name, Chromosome and position on chromosome. Without all three
    components GWhEAT will not be able to run successfully.  
      
5.  How is GWASapply\_rapid different than GWASapply? Both GWASapply and
    GWASapply\_rapid run the calculations the same way. The only
    difference is GWAsapply\_rapid will not give you any of the output
    automatically, which makes it run faster and allows the user the
    develop their own output.  
      
6.  Where can I get help with issues regarding GWhEAT? All questions
    should be directed to Samuel Prather at <stp4freedom@gmail.com>.  
      
7.  Are there other R packages that can perform GWAS? Yes, I would
    strongly recommend you check out the packaged GAPIT
    <http://www.zzlab.net/GAPIT/> as it has a much wider range of
    options and functions.

### **Further Information**

For more information on individual functions please see the
“Reference\_Manual.pdf” or type ?FUNCITON\_NAME into the R console,
this will pull up specific information of each function inside GWhEAT.
For example typing ?manhattan\_plot will pull of the help page with
details about the function that creates the Manhattan plots.

PopVar\_Exmaple
================

### The majority of this code is my pipeline I use to take a hapmap file, convert it to numeric, filter and impute it.

If you just want the code for what I use for kamiak, you can just skip
to the bottom or look at the R Script I provided.

# Packages needed

``` r
#Better for FDR function
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("impute")
#From the source
require(compiler) #for cmpfun
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt") #make sure compiler is running
source("http://zzlab.net/GAPIT/emma.txt")

#Better for FDR function
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  library(dplyr)
  library(Hmisc)
  library(WeightIt)
  library(mpath)
  library(glmnetUtils)
  library(glmnet)
  library(MASS)
  library(Metrics)
  library(stringr)
  library(lsa)
  library(keras)
  library(tensorflow)
  library(BMTME)
  library(plyr)
  library(data.table)
  library(bigmemory)
  library(biganalytics)
  library(ggplot2)
  library(tidyverse)
  library(knitr)
  library(cvTools)
  library(vcfR)
  library(compiler)
  library(gdata)
  library(BiocManager)
  library(impute)
  library(PopVar)
  library(BLR)
  library(sommer)
  library(emmeans)
  library(heritability)
  library(arm)
  library(optimx)
  library(purrr)
  library(psych)
  library(lme4)
  library(lmerTest)
  library(gridExtra)
  library(grid)
  library(readxl)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,bigmemory,biganalytics,dplyr,compiler)
```

### **Load Package**

Next using the code below download and install the package GWhEAT from
Github. The bottom two line of code in the chunk below make sure the
dependencies GWhEAT relies on are also downloaded and installed.

``` r
#install package
install_github("lfmerrick21/WhEATBreeders")
library(WhEATBreeders)#package name
```

### Read in phenotypic file and rename columns

``` r
Phenotype <- read.csv(file="F:\\OneDrive\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\MN Grant\\TCAP\\TCAP_EDX_GP.csv", header=T, sep=",", stringsAsFactors=F) #Fill-in phenotype file 
colnames(Phenotype)<-c("Genotype","Mg","P","K","Ca","Mn","Fe","Cu","Zn")
length(unique(Phenotype$Genotype))
Phenotype$Env="TCAP_MN"
Phenotype=Phenotype[,c(1,10,2:9)]
Genotype <- read_excel("F:\\OneDrive\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\MN Grant\\TCAP\\242_19192_nuc_hapmap.xlsx")
sum(is.na(Genotype))
```

### Get data in order

#### In this example the phenotype file consisted of the first colum of genotype and the other 8 as nutrients

``` r
TCAP_QC=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="KNN",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("Mg","P","K","Ca","Mn","Fe","Cu","Zn"),
                Study="TCAP",
                Outcome="Tested", #Tested or Untested
                #Trial=c("F5_2015","DH_2020","BL_2015_2020"),
                Trial=c("TCAP_MN"),
                Scheme="K-Fold",#K-Fold or VS
                Method="Two-Step", #Two-Step or #One-Step
                Messages=TRUE)
```

### PopVar need to be in a certain format

#### Here I focused on T9 which was zinc

You can do it for all columns but only needs to be done for the trait
you intend to use in PopVar

``` r
load(file="GBS_2_TCAP.RData")
GBS_2_TCAP_MN_Zn$pheno
View(GBS_2_TCAP_MN_Zn$numeric)
#Remove taxa column
num=GBS_2_TCAP_MN_Zn$numeric[,-1]
#Convert numeric alleles to -1,0,1
num=apply(num,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
#Makes sure all columns are numeric
num=apply(num,2,as.numeric)
View(num)
num=data.frame(GBS_2_TCAP_MN_Zn$numeric$taxa,num)
num=rbind(colnames(num),num)
num=num[-1,]
colnames(num)[1]<-"taxa"
num$taxa=as.character(num$taxa)
num[1,1]="taxa"
View(num)
GBS_2_TCAP_MN_Zn$PopVar=as.matrix(num)
View(GBS_2_TCAP_MN_Zn$PopVar)

GBS_2_TCAP_MN_Zn$pheno=cbind(
GBS_2_TCAP_MN_Zn$pheno,
Ca=GBS_2_TCAP_MN_Ca$pheno[,2],
Cu=GBS_2_TCAP_MN_Cu$pheno[,2],
Fe=GBS_2_TCAP_MN_Fe$pheno[,2],
K=GBS_2_TCAP_MN_K$pheno[,2],
Mg=GBS_2_TCAP_MN_Mg$pheno[,2],
Mn=GBS_2_TCAP_MN_Mn$pheno[,2],
P=GBS_2_TCAP_MN_P$pheno[,2])

colnames(GBS_2_TCAP_MN_Zn$pheno)[1]="taxa"
save(GBS_2_TCAP_MN_Zn,file="TCAP_MN_Zn.RData")

dim(GBS_2_TCAP_MN_Zn$PopVar)
dim(GBS_2_TCAP_MN_Zn$pheno)
dim(GBS_2_TCAP_MN_Zn$map)
```

### This chunk is what I put in as my R Script for kamiak

This is where you can specify different models and such and can be
customized following the PopVar documentation.

``` r
load(file="TCAP_MN_Zn.RData")
library(PopVar)
TCAP_PopVar=pop.predict(G.in = as.matrix(GBS_2_TCAP_MN_Zn$PopVar), y.in =as.matrix(GBS_2_TCAP_MN_Zn$pheno), map.in = as.matrix(GBS_2_TCAP_MN_Zn$map),
                                                                  crossing.table = NULL, parents = "TP", tail.p = 0.1, nInd = 200,
                                                                  map.plot = T, min.maf = 0.01, mkr.cutoff = 0.50, entry.cutoff = 0.50,
                                                                  remove.dups = F, impute = "pass", nSim = 25, frac.train = 0.8,
                                                                  nCV.iter = 10, nFold = 5, nFold.reps = 10, nIter = 12000,
                                                                  burnIn = 3000, models = c("rrBLUP"), return.raw = T)
save(TCAP_PopVar,file = "TCAP_PopVar_NC.RData")
```

The PopVar output is a ton of lists which is hard to visualize and use.
It consisted of 31 columns due to the secondary selection for the other
7 traits than Zn

``` r
load(file="TCAP_PopVar_NC.RData")
TCAP_PopVar$CVs$Zn
TCAP_PopVar$predictions$Zn_param.df
TCAP_PopVar$preds.per.sim$Zn_param.df
TCAP_PopVar$models.chosen
TCAP_PopVar$markers.removed
TCAP_PopVar$entries.removed
try=data.frame(TCAP_PopVar$predictions$Zn_param.df)
try1=unlist(try$Par1)
try2=unlist(try$Par2)
try3=unlist(try$midPar.Pheno)
try4=unlist(try$midPar.GEBV)
try5=unlist(try$pred.mu)
try6=unlist(try$pred.mu_sd)
try7=unlist(try$pred.varG)
try8=unlist(try$pred.varG_sd)
try9=unlist(try$mu.sp_low)
try10=unlist(try$mu.sp_high)
try11=unlist(try$low.resp_Mg)
try12=unlist(try$low.resp_P)
try13=unlist(try$low.resp_Ca)
try14=unlist(try$low.resp_Mn)
try15=unlist(try$low.resp_Fe)
try16=unlist(try$low.resp_Cu)
try17=unlist(try$low.resp_Zn)
try18=unlist(try$high.resp_Mg)
try19=unlist(try$high.resp_P)
try20=unlist(try$high.resp_Ca)
try21=unlist(try$high.resp_Mn)
try22=unlist(try$high.resp_Fe)
try23=unlist(try$high.resp_Cu)
try24=unlist(try$high.resp_Zn)
try25=unlist(try$cor_w._Mg)
try26=unlist(try$cor_w._P)
try27=unlist(try$cor_w._Ca)
try28=unlist(try$cor_w._Mn)
try29=unlist(try$cor_w._Fe)
try30=unlist(try$cor_w._Cu)
try31=unlist(try$cor_w._Zn)
PopVar_Names=colnames(TCAP_PopVar$predictions$Zn_param.df)
Popvar_Zn=data.frame(try1,try2,try3,try4,try5,try6,try7,try8,try9,try10,try11,try12,try13,try14,try15,try16,try18,try19,try20,try21,try22,try23,try25,try26,try27,try28,try29,try30)
PopVar_Zn_Names=PopVar_Names[-c(17,24,31)]
colnames(Popvar_Zn)<-PopVar_Zn_Names
Popvar_Zn[1:10,1:28]

View(Popvar_Zn[1:10,1:28])
```

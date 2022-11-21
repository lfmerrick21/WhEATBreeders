User Manual and Tutorial for WhEATBreeders
================
Lance F. Merrick
November 22, 2022

<p>

 

</p>

<p>

 

</p>

<!-- test inserting image-->

<center>

![](WinterWheat.png)

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

#### Load Package

#### Phenotypic Data

#### Genotypic Data

#### Genomic Selection Tutorial

#### Genome-Wide Association Tutorial using GAPIT

#### Cross Prediction

#### Frequently asked questions

#### Further Information

<p>

 

</p>

<p>

 

</p>

<p>

 

</p>

### **Introduction**

WhEATBreeders was created to lower the bar for implementing genomic
selection models for plant breeders to utilize within their own breeding
programs. Not only does include functions for genotype quality control
and filtering, but it includes easy to use wrappers for the most
commonly use models in many scenarios with K-Fold cross-validation or
validation sets. You can also implement GWAS assisted genomic selection.
We created a full wrapper for quality control and genomci selection in
our function “WHEAT”. Additionaly we walk through the set up of
unrpelicated data using adjuste means and calculate cullis heritability.
We also go through multi-output and multi-trait wrappers for GWAS in
GAPIt. Finally we walk through cross-prediction using PopVar, rrBLUP,
and sommer.

### **Description of package functions**

For a full list of functions within WhEATBreeders see
“Reference\_Manual.pdf” this file contains not only the full list of
functions but also a description of each. The pdf also has each
functions arguments listed. And like with all R packages once
WhEATBreeders is installed and loaded you can type ?function\_name and
that specific function’s full descriptions will appear in the help tab
on RStudio.

### **Getting started**

First if you do not already have R and R studio installed on your
computer head over to <https://www.r-project.org/> and install the
version appropriate for you machine. Once R and R studio are installed
you will need to install the WhEATBreeders package since this is a
working package in it’s early stages of development it’s only available
through Github. To download files off Github first download and load the
library of the package “devtools” using the code below.

### **Deep Dive Into Code**

For a deep dive into all the code in this package and the inner working
of the functions, please review the file
WhEATBreeders\_DeepDive\_Into\_Code.Rmd. There is a lot including the
adjusted means for single plot trials. All the quality control and GS
models, GWAS, Manhattan plots, and even popvar tutorials. There are also
GBS and heterozygote calling pipelines available in the relevant folder.

# Packages needed

``` r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools)
library(devtools)
#Better for FDR function
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("impute")
#From the source
#require(compiler) #for cmpfun
#Only if you want the source code
#source("http://zzlab.net/GAPIT/GAPIT.library.R")
#source("http://zzlab.net/GAPIT/gapit_functions.txt") #make sure compiler is running
#source("http://zzlab.net/GAPIT/emma.txt")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(BGLR,
  rrBLUP,
  caret,
  tidyr,
  dplyr,
  Hmisc,
  WeightIt,
  mpath,
  glmnetUtils,
  glmnet,
  MASS,
  Metrics,
  stringr,
  lsa,
  keras,
  tensorflow,
  BMTME,
  plyr,
  data.table,
  bigmemory,
  biganalytics,
  ggplot2,
  tidyverse,
  knitr,
  cvTools,
  vcfR,
  compiler,
  gdata,
  PopVar,
  BLR,
  sommer,
  heritability,
  arm,
  optimx,
  purrr,
  psych,
  lme4,
  lmerTest,
  gridExtra,
  grid,
  readxl,
  devtools)
```

### **Load Package**

Next using the code below download and install the package WhEATBreeders
from Github. The bottom two line of code in the chunk below make sure
the dependencies WhEATBreeders relies on are also downloaded and
installed.

``` r
#install package
install_github("lfmerrick21/WhEATBreeders")
library(WhEATBreeders)#package name
```

#### **Genotypic Data**

##### Read in Genotype Data and/or Phenotype Data

``` r
library(data.table)
Genotype<-fread("F:\\OneDrive\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Genomic Selection\\Genomic Selection Pipeline\\Jason_GBS\\WAC_2016-2020_production_filt.hmp.txt",fill=TRUE)
Phenotype<-fread("F:\\OneDrive\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Genomic Selection\\Genomic Selection Pipeline\\GS-Complex-Traits\\BL_EM_Pheno.csv",header=T)
```

#### Phenotypic data wide to long format

``` r
Phenotype=Phenotype[,-1]
Phenotype1=Phenotype %>%
  pivot_longer(!Genotype, names_to = "Env", values_to = "EM")
Phenotype=Phenotype1
```

#### You can conduct just quality control on phenotype and genotype data using the wrapper WHEAT function without running any GS models

##### Make QC=TRUE and GS=FALSE

``` r
LIND_QC=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="Two-Step", #Two-Step or #One-Step
                Messages=TRUE)
```

### **START HERE IF GENOTYPE DATA IS DONE**

### **Genomic Selection Tutorial**

### Optional Input

#### if you want to use the numeric files

##### GDre (Numeric Matrix from output)

##### GT (Taxa file from output)

##### GIre (Map file from output)

##### GBS\_Train (GBS RData output you want to use for training population)

##### GBS\_Predict (GBS RData output you want to use for validation set if you’re using VS)

### Multipe Trait Output

#### Method: “Two-Step”

#### Scheme: Cross-Validation (“K-Fold”) or Validation Sets (“VS”)

#### fold (Number of folds in K-Fold)

#### Training (Select trial for training population)

#### Prediction (Select trial if using validation sets)

#### Outcome: “Tested” or “Untested”

##### Tested allows predictions and accuracy and error comparisons

##### Untested for new predictions without comparisons

#### Packages: “rrBLUP”, “MAS”, “BGLR”, “caret”, “GLM”,

#### Type: “Regression” or “Classification” with caret package

#### Models: rrBLUP, MAS, BGLR(BayesA,BayesB,BayesC,BayesL,BayesRR,RKHS, Ordinal), caret(any model available in the caret package), GAPIT(gBLUP,cBLUP, or sBLUP)

#### nIter (Number of iterations for BGLR models)

#### burnIn (Number of burn-ins for BGLR models)

#### With/without Principal Components (“PC”: Number included or NULL)

#### With/without Covariates (“CV”: CV set (Genotypes and CVs) or NULL)

#### Phenotype transformation (“none”,“sqrt”,“log”,“boxcox”)

#### Subsampling of markers (markers: Number of markers to sample or NULL)

#### Kernel (“Markers”,“Linear”,Gaussian“,”Polynomial“,”Sigmoid“,”Arc-cosine with 1 layer or multiple layers:nl (“AK1”,“AKL”),“Exponential”)

#### Kernel (caret package allows above kernels with addition of “VanRaden”,“PC”,“Zhang”,“Endelman”(A-mat),)

#### nl (Number of layers for arc-cosine kernel)

#### degree (degree for polynomial kernel)

#### Sparse Matrix (Sparese=TRUE or FALSE)

#### m (Number of lines to keep in sparse matrix)

#### Replications (Number of Replications)

#### GWAS-Assisted GS (GAGS:TRUE or FALSE)

#### PCA.total (Number of PCs for GAGS and GAPIT models)

#### QTN (Number of top QTN or markers to include for GAGS)

#### GWAS (GWAS model for GAGS c("BLINK) or any GAPIT model)

#### threshold (“Bonferonni” or “none” for GAGS)

#### alpha (0.5, or any other threshold for Bonferonni for GAGS)

#### Messages (TRUE or FALSE)

#### All of the above can be done with the wrapper function WHEAT

#### Option shown below displays the function after QC and just GS with the GBS RData as input for validation sets using Two-step rrBLUP

#### Input using numeric data.

``` r
load(file="Tutorial_Filt_Imputed.RData")
#Input
F515_DH20=WHEAT(Phenotype=Phenotype,
                GDre=GDre,
                GT=GT,
                GIre=GIre,
                #GS Info
                Type="Regression",
                Replications=2,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Untested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="VS",
                Method="Two-Step",
                Package="rrBLUP",
                model="rrBLUP",
                Kernel="Markers",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation="none",
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold="none",
                GE=TRUE,
                UN=FALSE,
                GE_model="MTME",
                sampling="up",
                repeats=5,
                method="repeatedcv",
                digits=4,
                nCVI=5,
                Messages=TRUE)
```

#### Input using GBS RData.

``` r
load(file ='GBS_2_Tutorial.RData')
#Input
F515_DH20=WHEAT(Phenotype=Phenotype,
                GBS_Train=GBS_2_F5_2015_EM,
                GBS_Predict=GBS_2_Untested_DH_2020_EM,
                #GS Info
                Type="Regression",
                Replications=2,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Untested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="VS",
                Method="Two-Step",
                Package="rrBLUP",
                model="rrBLUP",
                Kernel="Markers",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation="none",
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold="none",
                GE=TRUE,
                UN=FALSE,
                GE_model="MTME",
                sampling="up",
                repeats=5,
                method="repeatedcv",
                digits=4,
                nCVI=5,
                Messages=TRUE)
```

### Full Wrapper for Quality Control and Genomic Selection

#### Input using Phenotype with Genotype, Environment, and Trait(s)

#### Use hapmap as Genotype file

#### Set QC and GS to TRUE

#### Option shown below displays the function after QC and just GS with the GBS RData as input for validation sets using Two-step rrBLUP

``` r
F515_DH20=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=TRUE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                Filter_Ind=TRUE,
                Missing_Rate_Ind=0.80,
                #If QC Info is FALSE
                GDre=NULL,
                GT=NULL,
                GIre=NULL,
                GBS_Train=NULL,
                GBS_Predict=NULL,
                Matrix=NULL,
                #GS Info
                Type="Regression",
                Replications=2,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Untested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="VS",
                Method="Two-Step",
                Package="rrBLUP",
                model="rrBLUP",
                Kernel="Markers",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation="none",
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold="none",
                GE=TRUE,
                UN=FALSE,
                GE_model="MTME",
                Messages=TRUE)
```

### One-Step Genomic Selection

#### Input using numeric data.

##### Gaussian Kernel output: Kernel=“Gaussian”

``` r
load(file="Tutorial_Filt_Imputed.RData")
#Input
F515_DH20=WHEAT(Phenotype=Phenotype,
                GDre=GDre,
                GT=GT,
                GIre=GIre,
                #GS Info
                Type="Regression",
                Replications=2,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Untested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",
                Method="One-Step",
                Package="BGLR",
                Kernel="Gaussian",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation="none",
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold="none",
                GE=TRUE,
                UN=FALSE,
                GE_model="MTME",
                sampling="up",
                repeats=5,
                method="repeatedcv",
                digits=4,
                nCVI=5,
                Messages=TRUE)
```

#### Input using GBS RData using Matrix Input.

#### BGLR using MTME for

##### Gaussian Kernel output: Kernel=“Gaussian”

``` r
load(file ='GBS_2_Tutorial.RData')
#Input
F515_DH20=WHEAT(Phenotype=Phenotype,
                Matrix=Matrix,
                #GS Info
                Type="Regression",
                Replications=2,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Untested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",
                Method="One-Step",
                Package="BGLR",
                Kernel="Gaussian",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation="none",
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold="none",
                GE=TRUE,
                UN=FALSE,
                GE_model="MTME",
                sampling="up",
                repeats=5,
                method="repeatedcv",
                digits=4,
                nCVI=5,
                Messages=TRUE)
```

### Individual Functions

#### We can get the proper matrices with the wrapper WHEAT function

##### Make QC=TRUE and GS=FALSE

##### Gaussian Kernel output: Kernel=“Gaussian”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=FALSE
                GE_model="MTME"
                Messages=TRUE)
```

#### Bayesian without GE: model=“ST”

#### Output with 3-Way Covariance Matrices

#### GE\_model=“BMTME”

#### UN=“FALSE”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=FALSE
                GE_model="BMTME"
                Messages=TRUE)
```

#### Bayesian with 3-Way Covariance: model=“BME” and UN=FALSE

#### Output with Factor Analytic Matrix

#### GE\_model=“MTME”

#### UN=“FALSE”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=FALSE
                GE_model="MTME"
                Messages=TRUE)
```

#### Bayesian with Factor Analytic Model: model=“BME” and UN=FALSE

#### Output with Marker-Environment Interaction Matrix

#### GE\_model=“MEI”

#### UN=“TRUE”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=TRUE
                GE_model="MEI"
                Messages=TRUE)
```

#### Bayesian with Marker-Environment Interaction Model: model=“BME” and function is MTME\_UN

#### One-Step With Machine Linear MLP Model

#### We can get the proper matrices with the wrapper WHEAT function

##### Make QC=TRUE and GS=FALSE

##### Gaussian Kernel output: Kernel=“Gaussian”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=FALSE
                GE_model="MTME"
                Messages=TRUE)
```

#### MLP without GE: model=“ST”

#### Output with 3-Way Covariance Matrices

#### GE\_model=“BMTME”

#### UN=“FALSE”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=FALSE
                GE_model="BMTME"
                Messages=TRUE)
```

#### MLP with 3-Way Covariance: model=“BME” and UN=FALSE

#### Output with Factor Analytic Matrix

#### GE\_model=“MTME”

#### UN=“FALSE”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE
                UN=FALSE
                GE_model="MTME"
                Messages=TRUE)
```

#### MLP with Factor Analytic Model: model=“BME” and UN=FALSE

#### Output with Marker-Environment Interaction Matrix

#### GE\_model=“MEI”

#### UN=“TRUE”

``` r
LIND_QC_One_Step=WHEAT(Phenotype=Phenotype,
                Genotype=Genotype,
                QC=TRUE,
                GS=FALSE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                #Do not remove individuals
                Filter_Ind=FALSE,
                Missing_Rate_Ind=0.80,
                Trait=c("EM"),
                Study="Tutorial",
                Outcome="Tested",
                Trial=c("F5_2015","DH_2020"),
                Scheme="K-Fold",#K-Fold or VS
                Method="One-Step", #Two-Step or #One-Step
                Kernel="Gaussian",
                folds = 5,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                GE=TRUE,
                UN=TRUE,
                GE_model="MEI"
                Messages=TRUE)
```

#### MLP with Marker-Environment Interaction Model: model=“BME” and function is MLP\_CV\_UN

### **Frequently asked questions**

1.  Where can I get WhEATBreeders? WhEATBreeders is currently only
    available on github via Lance F. Merrick’s public repository found
    at <https://github.com/lfmerrick21/WhEATBreeders>.
2.  How to download WhEATBreeders? Please refer the “Getting started”
    section above where I show how to download and install WhEATBreeders
    using devtools.
3.  Where can I get help with issues regarding WhEATBreeders? All
    questions should be directed to Lance F. Merrick at
    <lance.merrick21@gmail.com>.

### **Further Information**

For more information on individual functions please see the
“Reference\_Manual.pdf” or type ?FUNCITON\_NAME into the R console,
this will pull up specific information of each function inside
WhEATBreeders. For example typing ?manhattan\_plot will pull of the help
page with details about the function that creates the Manhattan plots.

User Manual and Tutorial for WhEATBreeders
================
Lance F. Merrick
February 2, 2022

<p>

 

</p>

<p>

 

</p>

<!-- test inserting image-->

<center>

![](Wheat.jpg)

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

``` r
#file.choose()
rm(list = ls(all.names = TRUE))
gc()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,bigmemory,biganalytics,dplyr,compiler)
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
library(devtools)
library(emmeans)
library(heritability)
library(arm)
library(optimx)
library(purrr)
library(psych)
library(lme4)
library(lmerTest)
library(gridExtra)
require(compiler) #for cmpfun
#source("http://zzlab.net/GAPIT/GAPIT.library.R")
#source("http://zzlab.net/GAPIT/gapit_functions.txt") #make sure compiler is running
#source("http://zzlab.net/GAPIT/emma.txt")

#Better for FDR function
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)
```

### **Load Package**

Next using the code below download and install the package GWhEAT from
Github. The bottom two line of code in the chunk below make sure the
dependencies GWhEAT relies on are also downloaded and installed.

``` r
#install package
install_github("lfmerrick21/WhEATBreeders")
library(GWhEAT)#package name
#Load dependencies 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,knitr,gridExtra)
```

In order to have an effective tutorial you’ll need some data to play
with, this code below downloads and loads into your environment data
used for this tutorial. \#\#\#\# Read in GBS data

``` r
if (!file.exists("GAPIT_Tutorial_Data.zip"))
{
  download.file("http://zzlab.net/GAPIT/GAPIT_Tutorial_Data.zip", destfile = "GAPIT_Tutorial_Data.zip")
  unzip("GAPIT_Tutorial_Data.zip")
}
download.file("http://zzlab.net/GAPIT/data/CROP545_Covariates.txt", destfile = "CROPS545_Covariates.txt")
download.file("http://zzlab.net/GAPIT/data/CROP545_Phenotype.txt", destfile = "CROPS545_Phenotype.txt")
# Import the GAPIT demo data genotypes
gt_scan <- data.frame(read.table("GAPIT_Tutorial_Data/mdp_numeric.txt", header = T, stringsAsFactors = F, sep = "\t", nrows = 1))
classes <- sapply(gt_scan, class)
genotypes <- data.frame(read.table("GAPIT_Tutorial_Data/mdp_numeric.txt", header = T, row.names = 1, colClasses = classes, stringsAsFactors = F, sep = "\t"))

GM <- read.table("GAPIT_Tutorial_Data/mdp_SNP_information.txt", header = T, stringsAsFactors = F, sep = "\t")
CV <- read.table("CROPS545_Covariates.txt", header = T, stringsAsFactors = F, sep = "\t")
phenotypes <- read.table("CROPS545_Phenotype.txt", header = T, stringsAsFactors = F, sep = "\t")


hapmap <- fread("https://zzlab.net/GAPIT/data/mdp_genotype_test.hmp.txt",fill=TRUE)
```

#### Read in Phenotypic Data

``` r
#Read in phenotypic files
em_trials=read.csv("C:\\Users\\lance\\OneDrive - Washington State University (email.wsu.edu)\\Documents\\Genomic Selection\\Genomic Selection Pipeline\\Selected Trials_C\\Selected Trials_Emergence\\Emergence_Trials.csv",header=TRUE)
str(em_trials)
em_trials$bloc=as.factor(em_trials$bloc)
em_trials$checks=as.factor(em_trials$checks)
em_trials$prow=as.factor(em_trials$prow)
em_trials$pcol=as.factor(em_trials$pcol)
em_trials$ibloc=as.factor(em_trials$ibloc)
em_trials$year1=as.factor(em_trials$year1)
em_trials$r_expt=as.factor(em_trials$r_expt)
colnames(em_trials)[8]<-c("Name")

#We have to create new identifiers according to the model indicated in the review paper I sent you.
em_trials$new.ind=em_trials$checks
em_trials=em_trials %>% mutate(check.ind = recode(checks,"1"="0", "0"="1"))
em_trials$check.ind=as.factor(em_trials$check.ind)
em_trials$new.ind=as.factor(em_trials$new.ind)
str(em_trials)
View(em_trials)

#I just subset the different trials I needed
lq15=subset(em_trials,em_trials$r_expt==levels(em_trials$r_expt)[2])
lq17=subset(em_trials,em_trials$r_expt==levels(em_trials$r_expt)[3])
lq18=subset(em_trials,em_trials$r_expt==levels(em_trials$r_expt)[4])
lq19=subset(em_trials,em_trials$r_expt==levels(em_trials$r_expt)[5])
lq15_17=em_trials[em_trials$r_expt == "2015 QAM Plots Lind" | em_trials$r_expt == "2017 QAM Plots Lind", ]

lq15_18=em_trials[em_trials$r_expt == "2015 QAM Plots Lind" | em_trials$r_expt == "2017 QAM Plots Lind" |em_trials$r_expt == "2018 QAM Rows Lind", ]
lq15_19=em_trials[em_trials$r_expt == "2015 QAM Plots Lind" | em_trials$r_expt == "2017 QAM Plots Lind" |em_trials$r_expt == "2018 QAM Rows Lind" |em_trials$r_expt == "2019 QAM Rows Lind", ]
```

### **Phenotypic Data**

#### Cullis Heritability (From Shantel Martinez Github)

##### Single Environment Including genotype

``` r
lq15_bp <- lmer(emer~check.ind + (1|Name:new.ind) + (1|ibloc), data=lq15)
Cullis_H2(lq15_bp)
summary(lq15_bp)
anova(lq15_bp)
ranova(lq15_bp)
```

##### Multi Environment including genotype

``` r
lq15_17_bp <- lmer(emer~check.ind + (1|Name:new.ind) + (1|ibloc)+(1|r_expt)+check.ind:r_expt + (1|Name:new.ind:r_expt) + (1|ibloc:r_expt), data=lq15_17)
Cullis_H2(lq15_17_bp)
summary(lq15_17_bp)
anova(lq15_17_bp)
ranova(lq15_17_bp)
```

\#\#\#Singularity
Problems

##### If you run into singularity or convergence problems you can always add in an lmerControl and then mess around with the maxfun number to get rid of the problem.

``` r
lb1_2_bp <- lmer(emer~check.ind + (1|Name:new.ind) + (1|ibloc)+(1|r_expt)+check.ind:r_expt + (1|Name:new.ind:r_expt) + (1|ibloc:r_expt), data=lbl1_2,
            control=lmerControl(optimizer="Nelder_Mead",
                                 optCtrl=list(maxfun=1e3)))
```

#### Adjusted Means

##### Single Environment

``` r
#Adjustments
lq15_r <- lm(emer~ibloc+checks, data=lq15)
#residuals(lq15_r)
lq15=lq15 %>% mutate(ap_adj = emer+lq15_r$residuals)
lq15_ap=aggregate(lq15[,c(21)],list(lq15$Name),mean)
colnames(lq15_ap)<-c("Name","emer_f15_ap")
```

##### Multi Environment

``` r
lq15_17_r <- lm(emer~ibloc+checks+r_expt+ibloc*r_expt+checks*r_expt, data=lq15_17)
#lq15_17_r$residuals
lq15_17=lq15_17 %>% mutate(ap_adj = emer+lq15_17_r$residuals)
lq15_17_ap=aggregate(lq15_17[,c(21)],list(lq15_17$Name),mean)
colnames(lq15_17_ap)<-c("Name","emer_15_17_ap")
```

#### True version of Cullis heritability code that I’m working on which has been adapted from the Comparison of heritability papers and the github “<https://github.com/PaulSchmidtGit/Heritability>”. I’m still working on this code, but I figured I’d give you an overview.

##### Single Environment

``` r
g.ran <- lmer(data    = lq15,
              formula = emer~check.ind + (1|Name:new.ind) + (1|ibloc))

### handle model estimates
# to my knowledge, lme4 does not offer a function to
# extract variance-covariance-matrices for BLUPs (a.k.a. prediction error variance [PEV] matrix).
# therefore, I here manually reconstruct mixed model equation for this specific example.
# notice that this solution therefore only works for this specific model!

vc <- g.ran %>% VarCorr %>% as_tibble # extract estimated variance components (vc)

# R = varcov-matrix for error term
n <- g.ran %>% summary %>% pluck(residuals) %>% length # numer of observations
vc_e <- vc %>% filter(grp=="Residual") %>% pull(vcov)  # error vc
R    <- diag(n)*vc_e                                   # R matrix = I_n * vc_e

# G = varcov-matrx for all random effects
# subset of G regarding genotypic effects
n_g  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("Name:new.ind") # number of genotypes
vc_g <- vc %>% filter(grp=="Name:new.ind") %>% pull(vcov)              # genotypic vc
G_g  <- diag(n_g)*vc_g                                        # gen part of G matrix = I * vc.g

# subset of G regarding incomplete block effects
n_b  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("ibloc") # number of incomplete blocks
vc_b <- vc %>% filter(grp=="ibloc") %>% pull(vcov)              # incomplete block vc
G_b  <- diag(n_b)*vc_b                                              # incomplete block part of G matrix = I * vc.b

G <- bdiag(G_g, G_b) # G is blockdiagonal with G.g and G.b in this example
G <- G_g
# Design Matrices
X <- g.ran %>% getME("X") %>% as.matrix # Design matrix fixed effects
Z <- g.ran %>% getME("Z") %>% as.matrix # Design matrix random effects

# Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
C11 <- t(X) %*% solve(R) %*% X
C12 <- t(X) %*% solve(R) %*% Z
C21 <- t(Z) %*% solve(R) %*% X
C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 

C <- rbind(cbind(C11, C12),  
           cbind(C21, C22)) %>% as.matrix # Combine components into one matrix C

# Mixed Model Equation Solutions 
C_inv <- C %>% solve# Inverse of C
dim(C_inv)
View(C_inv)
colnames(C_inv)
unique(length(levels(lq15$Name)))
lq15$Name=droplevels(lq15$Name)
C22_g <- C_inv[3:483, 3:483] # subset of C.inv that refers to genotypic BLUPs
head(C_inv)
# Mean variance of BLUP-difference from C22 matrix of genotypic BLUPs
one        <- matrix(1, nrow=n_g, ncol=1)      # vector of 1s
P_mu       <- diag(n_g, n_g) - one %*% t(one)  # P_mu = matrix that centers for overall-mean
vdBLUP_sum <- psych::tr(P_mu %*% C22_g)        # sum of all variance of differences = trace of P_mu*C22_g
vdBLUP_avg <- vdBLUP_sum * (2/(n_g*(n_g-1)))   # mean variance of BLUP-difference = divide sum by number of genotype pairs

### H2 Cullis
H2Cullis <- 1 - (vdBLUP_avg / 2 / vc_g)
H2Cullis #0.4911897
```

##### Multi Environment

``` r
g.ran <- lmer(data    = lq15_17,
              formula = emer~check.ind + (1|Name:new.ind) + (1|ibloc)+(1|r_expt)+check.ind:r_expt + (1|Name:new.ind:r_expt) + (1|ibloc:r_expt))

### handle model estimates
# to my knowledge, lme4 does not offer a function to
# extract variance-covariance-matrices for BLUPs (a.k.a. prediction error variance [PEV] matrix).
# therefore, I here manually reconstruct mixed model equation for this specific example.
# notice that this solution therefore only works for this specific model!

vc <- g.ran %>% VarCorr %>% as_tibble # extract estimated variance components (vc)

# R = varcov-matrix for error term
n <- g.ran %>% summary %>% pluck(residuals) %>% length # numer of observations
vc_e <- vc %>% filter(grp=="Residual") %>% pull(vcov)  # error vc
R    <- diag(n)*vc_e                                   # R matrix = I_n * vc_e

# G = varcov-matrx for all random effects
# subset of G regarding genotypic effects
n_g  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("Name:new.ind") # number of genotypes
vc_g <- vc %>% filter(grp=="Name:new.ind") %>% pull(vcov)              # genotypic vc
G_g  <- diag(n_g)*vc_g                                        # gen part of G matrix = I * vc.g

# subset of G regarding incomplete block effects
n_nr  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("Name:new.ind:r_expt") # number of incomplete blocks
vc_nr <- vc %>% filter(grp=="Name:new.ind:r_expt") %>% pull(vcov)              # incomplete block vc
G_nr <- diag(n_nr)*vc_nr   


n_br  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("ibloc:r_expt") # number of incomplete blocks
vc_br <- vc %>% filter(grp=="ibloc:r_expt") %>% pull(vcov)              # incomplete block vc
G_br  <- diag(n_br)*vc_br                                              # incomplete block part of G matrix = I * vc.b

n_b  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("ibloc") # number of incomplete blocks
vc_b <- vc %>% filter(grp=="ibloc") %>% pull(vcov)              # incomplete block vc
G_b <- diag(n_b)*vc_b  

n_r  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("r_expt") # number of incomplete blocks
vc_r <- vc %>% filter(grp=="r_expt") %>% pull(vcov)              # incomplete block vc
G_r <- diag(n_r)*vc_r   

G <- bdiag(G_g,G_b,G_r,G_nr,G_br) # G is blockdiagonal with G.g and G.b in this example
# Design Matrices
X <- g.ran %>% getME("X") %>% as.matrix # Design matrix fixed effects
Z <- g.ran %>% getME("Z") %>% as.matrix # Design matrix random effects
dim(X)
dim(Z)
# Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
C11 <- t(X) %*% solve(R) %*% X
C12 <- t(X) %*% solve(R) %*% Z
C21 <- t(Z) %*% solve(R) %*% X
C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 
dim(C11)
dim(C12)
dim(C21)
dim(C22)
dim(G)
C <- rbind(cbind(C11, C12),  
           cbind(C21, C22)) %>% as.matrix # Combine components into one matrix C

# Mixed Model Equation Solutions 
C_inv <- C %>% solve# Inverse of C
dim(C_inv)
colnames(C_inv)[1447]
lq15_17$Name=droplevels(lq15_17$Name)
unique(length(levels(lq15_17$Name)))
unique(length(levels(lq15_17$ibloc)))
unique(length(levels(lq15_17$checks)))

C22_g <- C_inv[967:1447, 967:1447] # subset of C.inv that refers to genotypic BLUPs
# Mean variance of BLUP-difference from C22 matrix of genotypic BLUPs
one        <- matrix(1, nrow=n_g, ncol=1)      # vector of 1s
P_mu       <- diag(n_g, n_g) - one %*% t(one)  # P_mu = matrix that centers for overall-mean
vdBLUP_sum <- psych::tr(P_mu %*% C22_g)        # sum of all variance of differences = trace of P_mu*C22_g
vdBLUP_avg <- vdBLUP_sum * (2/(n_g*(n_g-1)))   # mean variance of BLUP-difference = divide sum by number of genotype pairs

### H2 Cullis
H2Cullis <- 1 - (vdBLUP_avg / 2 / vc_g)
H2Cullis #0.4911897
```

#### **Genotypic Data**

##### Read in 2020 Lind Taxa

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

##### IF VCF

##### Read in VCF

##### Currently do not have a VCF file in the folder

``` r
#setwd("F:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/Documents/Genomic Selection/Genomic Selection Pipeline/GWAS Pipeline")
Xvcf=read.vcfR("YR_Numeric_COMP1.vcf")
Xvcf@fix
Xvcf@gt
```

#### Convert to Numeric

``` r
Xbgl=data.frame(Xvcf@fix,Xvcf@gt)
genotypes_num<- VCF_2_numeric(Xbgl)
GDE=genotypes_num$genotypes
GIE=genotypes_num$marker_map
View(GDE[1:10,1:10])
View(GIE[1:10,1:3])
```

The majority of this code is my pipeline I use to take a hapmap file,
convert it to numeric, filter and impute it. If you just want the code
for what I use for kamiak, you can just skip to the bottom or look at
the R Script I provided.

#### Filter GBS for Lind 2020

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

#### Convert to Numeric

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

#### Remove Missing Data \>80%, MAF \<5% and Monomorphic Markers

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

#### Filter

##### Remove makers with more than 20% missing

``` r
mr=calc_missrate(GDE)
mr_indices <- which(mr > 0.20)
GDEmr = GDE[,-mr_indices]
GIEmr = GIE[-mr_indices,]
dim(GDEmr)
dim(GIEmr)
```

##### Remove Monomorphic Markers

``` r
maf <- calc_maf_apply(GDEmr, encoding = c(0, 1, 2))
mono_indices <- which(maf ==0)
GDEmo = GDEmr[,-mono_indices]
GIEmo = GIEmr[-mono_indices,]
dim(GDEmo)
dim(GIEmo)
```

##### Remove MAF \<5%

``` r
maf <- calc_maf_apply(GDEmo, encoding = c(0, 1, 2))
mono_indices <- which(maf <0.05)
GDEmf = GDEmo[,-mono_indices]
GIEmf = GIEmo[-mono_indices,]
dim(GDEmf)
dim(GIEmf)
```

#### Imputation

##### Imputation Using BEAGLE

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

##### Or Impute with KNN

``` r
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
mono_indices <- which(maf <0.05)
GDEre = myGD_imp[,-mono_indices]
GIEre = GIEmf[-mono_indices,]
dim(GDEre)
dim(GIEre)
```

##### Save Files for easy reproducibility

``` r
fwrite(GDEre,"GDE_numeric.csv")
colnames(GIEre)
fwrite(GIEre[,1:3],"GIE_map.csv")
```

### **START HERE IF GENOTYPE DATA IS DONE**

#### This gets into my own files This section gets your data in order but using a function

#### Function to get data in order as above chunck but integrated into a function to also combine all files for a nice list with PCs included.

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

### **Genomic Selection Tutorial**

### Cross-Validation

#### GS without CV

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

#### HSD and ggplot comparison

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

#### GS with CV

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

#### Without CV

``` r
LND_20=Caret_Models_Mean_Matrix(genotypes=GBS_LND_20$geno,
                                phenotype= GBS_LND_20$pheno[,c(1,2)],
                                model="svmRadial",
                                Matrix="VanRaden",
                                sampling="up")#Unbalanced
```

#### With CV

``` r
LND_20_CV=Caret_Models_Mean_Matrix(genotypes=cbind(GBS_LND_20_CV$CV[,c(2:5)],GBS_LND_20_CV$geno),
                                   phenotype= GBS_LND_20_CV$pheno[,c(1,2)],
                                   model="svmRadial",
                                   Matrix="VanRaden",
                                   sampling="up")
```

### Validation GS

#### Without CV

``` r
test_all_models_BLUP_vs_pc_mean_recode(train_genotypes = GBS_BQALL_comp2$geno, 
                                       train_phenotype = GBS_BQALL_comp2$pheno[,c(1,2)],
                                       train_PCA=GBS_BQALL_comp2$PC[,1:3],
                                       test_genotypes = GBS_BQALL_comp4$geno,
                                       test_phenotype = GBS_BQALL_comp4$pheno[,c(1,4)],
                                       test_PCA=GBS_BQALL_comp4$PC[,1:3])
```

#### With CV

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

#### R Script for Kamiak

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

#### R Script for HPC with CV

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

### **Genome-Wide Association Tutorial using GAPIT**

#### Model Selection in GAPIT (aka PC number selection)

``` r
#QAM
myGAPIT_MS<- GAPIT(Y = GBS_qam_adj18$pheno[,c(1,18)],
                         GD = GBS_qam_adj18$numeric,
                         GM = GBS_qam_adj18$map,
                         model = "MLM",
                         PCA.total=3,
                         Model.selection = TRUE,
                         file.output=T)
```

#### GAPIT

``` r
myGAPIT_BLINK_qam_adj_23<- GAPIT(Y = GBS_qam_adj23$pheno[,c(1,23)],
                         GD = GBS_qam_adj23$numeric,
                         GM = GBS_qam_adj23$map,
                         PCA.total=3,
                         model = "BLINK",
                         file.output=T)
#I like saving the results in a R data file.
save(myGAPIT_BLINK_qam_adj_23,file="GWAS_EM_QAM_ADJ_BLINK_ZZ_Redo_5_11.RData")

myGAPIT_BLINK_adj_13<- GAPIT(Y = GBS_adj13$pheno[,c(1,13)],
                         GD = GBS_adj13$numeric,
                         GM = GBS_adj13$map,
                         PCA.total=3,
                         model = "BLINK",
                         file.output=T)
save(myGAPIT_BLINK_adj_13,file="GWAS_EM_adj_BLINK_ZZ_Redo_5_11.RData")

FarmCPU_qam_adj_23<- GAPIT(Y = GBS_qam_adj23$pheno[,c(1,23)],
                         GD = GBS_qam_adj23$numeric,
                         GM = GBS_qam_adj23$map,
                         PCA.total=3,
                         model = "FarmCPU",
                         file.output=T)
save(FarmCPU_qam_adj_23,file="GWAS_EM_QAM_ADJ_ZZ_FarmCPU_Redo_5_11.RData")

FarmCPU_adj_13<- GAPIT(Y = GBS_adj13$pheno[,c(1,13)],
                         GD = GBS_adj13$numeric,
                         GM = GBS_adj13$map,
                         PCA.total=3,
                         model = "FarmCPU",
                         file.output=T)
save(FarmCPU_adj_13,file="GWAS_EM_adj_ZZ_FarmCPU_Redo_5_11.RData")

myGAPIT_MLM_qam_adj_23<- GAPIT(Y = GBS_qam_adj23$pheno[,c(1,23)],
                         GD = GBS_qam_adj23$numeric,
                         GM = GBS_qam_adj23$map,
                         PCA.total=3,
                         model = c("MLM"),
                         file.output=T)
save(myGAPIT_MLM_qam_adj_23,file="GWAS_EM_QAM_ADJ_MLM_ZZ_Redo_5_11.RData")
myGAPIT_MLM_adj_13<- GAPIT(Y = GBS_adj13$pheno[,c(1,13)],
                         GD = GBS_adj13$numeric,
                         GM = GBS_adj13$map,
                         PCA.total=3,
                         model = c("MLM"),
                         file.output=T)
save(myGAPIT_MLM_adj_13,file="GWAS_EM_adj_MLM_ZZ_Redo_5_11.RData")
```

#### Marker effects and variation

##### This can be either FDR or Bonferonni

``` r
BLINK_Effects=Marker_Effects(Pheno=myY[,2],GWAS=myGAPIT_BLINK_qam_adj_23,alpha=0.05,correction="Bonferonni",messages=TRUE,model="BLINK")
BLINK_Effects
```

##### Read in Hapmap

``` r
load(file="EM_Hapmap_Redo.RData")
output[1:10,1:14]
#Extract Map and Allele Information
EM_GBS_SNPs=output[,1:5]
```

This version is much better but requires a hapmap. It creates a nice
table with alleles in which allows you to identify the favorable allele.
The numericalization makes the 2 the allele second in alphabetical
order. If it is A/T then A=0 and T=2. So if the marker effect is
negative such as -2, then A is the favorable allele because -2 is the
effect of the T allele.

#### Calculate Effects

#### This can be either FDR or Bonferonni

IF you prefer the online functions for gapit, you will need to remove
all functions from the website since they don’t all work. I use gapit
functions in the below function. So to make it work, download the Github
version if you did not do that above.

``` r
#This removes just functions.
#rm(list=lsf.str())
#install.packages("devtools")
#devtools::install_github("jiabowang/GAPIT3",force=TRUE)
#library(GAPIT3)
FDR_sig_B6 <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_BLINK_qam_adj_23,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs)
FDR_sig_B6
FDR_sig_F6 <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=FarmCPU_qam_adj_23,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs)
FDR_sig_F6
FDR_sig_M6 <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_MLM_qam_adj_23,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs,model="MLM")
FDR_sig_M6
FDR_sig_B1_bl <- Marker_Effects_Allele(Pheno=GBS_adj13$pheno[,c(13)],GWAS=myGAPIT_BLINK_adj_13,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs)
FDR_sig_B1_bl
FDR_sig_F1_bl <- Marker_Effects_Allele(Pheno=GBS_adj13$pheno[,c(13)],GWAS=FarmCPU_adj_13,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs)
FDR_sig_F1_bl
FDR_sig_M1_bl <- Marker_Effects_Allele(Pheno=GBS_adj13$pheno[,c(13)],GWAS=myGAPIT_MLM_adj_13,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs,model="MLM")
FDR_sig_M1_bl
```

#### Secondary Trait

This is for another trait that I don’t want to create a manhattan plot
but overlay significant markers within the other manhattan
pltos

``` r
FDR_sig_B1_Col_Len_blup=Marker_Effects_Allele(Pheno=GBS_qam_adj18$CV[,c(6)],GWAS=BLINK_Col_Length_qam_blup_18,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs,model="BLINK")
FDR_sig_F1_Col_Len_blup=Marker_Effects_Allele(Pheno=GBS_qam_adj18$CV[,c(6)],GWAS=FarmCPU_Col_Length_qam_blup_18,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs,model="FarmCPU")
FDR_sig_M1_Col_Len_blup=Marker_Effects_Allele(Pheno=GBS_qam_adj18$CV[,c(6)],GWAS=MLM_Col_Length_qam_blup_18,alpha=0.05,correction="FDR",messages=TRUE,Markers=EM_GBS_SNPs,model="MLM")
Col_lenth=union(FDR_sig_B1_Col_Len_blup$SNP,FDR_sig_F1_Col_Len_blup$SNP)
```

I reacreated the manahttan plots from the gapit code using ggplot which
allows for far more customization, and the ability to edit the plot
after it is initially made.

I did this to allow easy labeling and allow multiple traits/models from
two different GWAS results to be
overlaid.

#### First you need to change Numeric to Number and Letter Chromosomes

``` r
FarmCPU_qam_adj_23$GWAS=FarmCPU_qam_adj_23$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))  
  
FarmCPU_adj_13$GWAS=FarmCPU_adj_13$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))

myGAPIT_BLINK_qam_adj_23$GWAS=myGAPIT_BLINK_qam_adj_23$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))  
  
myGAPIT_BLINK_adj_13$GWAS=myGAPIT_BLINK_adj_13$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))

myGAPIT_MLM_qam_adj_23$GWAS=myGAPIT_MLM_qam_adj_23$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))  
  
myGAPIT_MLM_adj_13$GWAS=myGAPIT_MLM_adj_13$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))

myGAPIT_MLM_qam_adj_23$GWAS=cbind(FarmCPU_Col_qam_adj_23$GWAS[,1:3],myGAPIT_MLM_qam_adj_23$GWAS[,-c(1:3)])
myGAPIT_MLM_adj_13$GWAS=cbind(FarmCPU_Col_qam_adj_23$GWAS[,1:3],myGAPIT_MLM_adj_13$GWAS[,-c(1:3)])
```

\#\#\#\#Single Manhattan Plot

If you don’t want to overlay GWAS results do not use Multi options.
Third labels is for a third trait where I just want to show the
significant markers without plotting the GWAS results.

QTN is if you want to highlight certain markers with a vertical line.

``` r
manhattan_plot_Multi(FarmCPU_qam_adj_23$GWAS,labels =NULL,
                     model="FarmCPU",QTN_index=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP),
                     Multi=NULL,Multi_labels=NULL,
                     Third_Labels=Col_lenth,model_cor ="BLINK")

manhattan_plot_Multi(myGAPIT_BLINK_qam_adj_23$GWAS,labels =FDR_sig_B6,
                     model="BLINK",QTN_index=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Multi=NULL,Multi_labels=NULL,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "BLINK")
manhattan_plot_Multi(myGAPIT_MLM_qam_adj_23$GWAS,labels =FDR_sig_M6,
                     model="MLM",QTN_index=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Multi=NULL,Multi_labels=NULL,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "MLM")
```

#### Overall Mulit Plot

This is to overlay multiple GWAS resutls. Multi is the second GWAS
results and Multi\_labels if you want to label those differently.

``` r
intersect(FDR_sig_F6$SNP,FDR_sig_F1_bl$SNP)

fmmp23=manhattan_plot_Multi(FarmCPU_qam_adj_23$GWAS,labels =FDR_sig_F6,
                     model="FarmCPU",QTN_index=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Multi=FarmCPU_adj_13$GWAS,Multi_labels=FDR_sig_F1_bl,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "BLINK")
bmmp23=manhattan_plot_Multi(myGAPIT_BLINK_qam_adj_23$GWAS,labels =FDR_sig_B6,
                     model="BLINK",QTN_index=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Multi=myGAPIT_BLINK_adj_13$GWAS,Multi_labels=FDR_sig_B1_bl,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "BLINK")
mmp23=manhattan_plot_Multi(myGAPIT_MLM_qam_adj_23$GWAS,labels =FDR_sig_M6,
                     model="MLM",QTN_index=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Multi=myGAPIT_MLM_adj_13$GWAS,Multi_labels=FDR_sig_M1_bl,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "MLM")


mmp23=mmp23+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x= element_blank())
fmmp23=fmmp23+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x= element_blank())
bmmp23=bmmp23
```

#### Zhiwu’s version

This allows for labeling significant markers that are found in all
models with a single large label at the top manahattan plot.

QTN\_index2 is differnt then above. It’s a different color that I used
if a marker was found in 2 models.

QTN\_index3 is like QTN\_index2 but for 3 or more models.

Sig\_L and Sig\_multi allow the ability to create different threshold
lines for the two different GWAS results since they may not be the same
especially since FDR is based on p-values. Whereas for bonferonni, if
the number of markers are the same so will the threshold line plotted.

``` r
intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)
intersect(FDR_sig_F6$SNP,FDR_sig_M6$SNP)
intersect(FDR_sig_B6$SNP,FDR_sig_M6$SNP)

fmmp23=manhattan_plot_Multi(FarmCPU_qam_adj_23$GWAS,labels =NULL,
                     model="FarmCPU",
                     QTN_index2=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[1],
                     
                     QTN_index3=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Sig_L=FDR_sig_F6,
                     Sig_Multi=FDR_sig_F1_bl,
                     
                     ZZ_label =intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP),
                     
                     Multi=FarmCPU_adj_13$GWAS,
                     Multi_labels=NULL,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "BLINK")

bmmp23=manhattan_plot_Multi(myGAPIT_BLINK_qam_adj_23$GWAS,labels =NULL,
                     model="BLINK",
                     QTN_index2=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[1],
                     
                     QTN_index3=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Sig_L=FDR_sig_B6,
                     Sig_Multi=FDR_sig_B1_bl,
                     
                     Multi=myGAPIT_BLINK_adj_13$GWAS,
                     Multi_labels=NULL,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "BLINK")

mmp23=manhattan_plot_Multi(myGAPIT_MLM_qam_adj_23$GWAS,labels =NULL,
                     model="MLM",
                     QTN_index2=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[1],
                     
                     QTN_index3=intersect(FDR_sig_F6$SNP,FDR_sig_B6$SNP)[2],
                     Sig_L=FDR_sig_M6,
                     Sig_Multi=FDR_sig_M1_bl,
                     
                     Multi=myGAPIT_MLM_adj_13$GWAS,
                     Multi_labels=NULL,
                     Third_Labels=Col_lenth,Trait_one="Diversity Panel",Trait_two="Breeding Lines",model_cor = "MLM")


fmmp23=fmmp23+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x= element_blank())
bmmp23=bmmp23+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x= element_blank())
mmp23=mmp23
```

#### Stacking of multiple manhattan plots.

``` r
library(grid)
grid.newpage()
grid.draw(rbind(
                ggplotGrob(fmmp23),
                ggplotGrob(bmmp23),
                ggplotGrob(mmp23),
                size = "last"))
```

#### Multi Trait GWAS based on sommer multivariate gwas but includes results for three different p-values for common, interaction, and full effects.

``` r
Gen_Table_JMP_DP23=left_join(GBS_qam_adj23$pheno[,c(1,23)],GBS_qam_adj23$CV[,c(1,6)],by="Genotype")
myGM23=GBS_qam_adj23$geno
myGM23=apply(myGM23,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
myGM23=apply(myGM23,2,as.numeric)
rownames(myGM23)=rownames(GBS_qam_adj23$geno)
A23 <- A.mat(myGM23)
colnames(Gen_Table_JMP_DP23)[2]<-"Y"

ansx23 <- GWAS(cbind(Y,DP_BLUP)~1,
random=~ vs(Genotype, Gu=A23, Gtc=unsm(2)),
rcov=~ vs(units, Gtc=unsm(2)),
data=Gen_Table_JMP_DP23,
M=myGM23,n.PC=3,
gTerm = "u:Genotype", verbose = TRUE)

MTMM23=sommer_MTMM(Y=Gen_Table_JMP_DP23,SNP_INFO=GBS_qam_adj23$map,model=ansx23,X=myGM23,A=A23)
```

#### Extract Values

``` r
myGAPIT_BLINK_qam_adj_23$GWAS=myGAPIT_BLINK_qam_adj_23$GWAS%>%mutate(Chromosome=recode(Chromosome,"1"="1A","2"="1B","3"="1D","4"="2A","5"="2B","6"="2D","7"="3A","8"="3B","9"="3D","10"="4A","11"="4B","12"="4D","13"="5A","14"="5B","15"="5D","16"="6A","17"="6B","18"="6D","19"="7A","20"="7B","21"="7D","22"="UN"))

#This is simply to creat teh same data frame to input our results so they can work with the other gapit functions we use.
myGAPIT_MTMM_qam_adj23_EM=myGAPIT_BLINK_qam_adj_23
myGAPIT_MTMM_qam_adj23_CL=myGAPIT_BLINK_qam_adj_23
myGAPIT_MTMM_qam_adj23_FULL=myGAPIT_BLINK_qam_adj_23
myGAPIT_MTMM_qam_adj23_IE=myGAPIT_BLINK_qam_adj_23
myGAPIT_MTMM_qam_adj23_COM=myGAPIT_BLINK_qam_adj_23

#INput our results into the data frames.
myGAPIT_MTMM_qam_adj23_EM$GWAS$P.value=MTMM23$pvals$Y.score
myGAPIT_MTMM_qam_adj23_CL$GWAS$P.value=MTMM23$pvals$DP_BLUP.score
myGAPIT_MTMM_qam_adj23_FULL$GWAS$P.value=MTMM23$pvals$pval_full
myGAPIT_MTMM_qam_adj23_IE$GWAS$P.value=MTMM23$pvals$pval_trait_specific
myGAPIT_MTMM_qam_adj23_COM$GWAS$P.value=MTMM23$pvals$pval_trait_common
```

#### Marker Effects

``` r
Bon_sig_MTMM23_EM <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_MTMM_qam_adj23_EM,alpha=0.05,correction="Bon",messages=TRUE,Markers=EM_GBS_SNPs)
Bon_sig_MTMM23_EM
Bon_sig_MTMM23_CL <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_MTMM_qam_adj23_CL,alpha=0.05,correction="Bon",messages=TRUE,Markers=EM_GBS_SNPs)
Bon_sig_MTMM23_CL
Bon_sig_MTMM23_FULL <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_MTMM_qam_adj23_FULL,alpha=0.05,correction="Bon",messages=TRUE,Markers=EM_GBS_SNPs)
Bon_sig_MTMM23_FULL
#Bon_sig_MTMM23_IE <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_MTMM_qam_adj23_IE,alpha=0.05,correction="Bon",messages=TRUE,Markers=EM_GBS_SNPs)
#Bon_sig_MTMM23_IE
Bon_sig_MTMM23_COM <- Marker_Effects_Allele(Pheno=GBS_qam_adj23$pheno[,c(23)],GWAS=myGAPIT_MTMM_qam_adj23_COM,alpha=0.05,correction="Bon",messages=TRUE,Markers=EM_GBS_SNPs)
Bon_sig_MTMM23_COM
```

#### Prepare for Venn Diagram

You can use these with the above manhattan plot functions, but the below
code is to create a venn diagram.

``` r
Extract_Sig=function(results,pop,year,model,cov){
nrep=nrow(results)
po=rep(pop,nrep)
yr=rep(year,nrep)
mo=rep(model,nrep)
cv=rep(cov,nrep)

aca=data.frame(Population=po,Year=yr,Model=mo,Covariate=cv,results)
#rownames(aca)<-names(ac)
return(aca)
}
```

#### Combine results

``` r
Bon_M23_EM=Extract_Sig(Bon_sig_MTMM23_EM,"DP","2015-2018","MTMM","EM")
Bon_M22_CL=Extract_Sig(Bon_sig_MTMM22_CL,"DP","2015-2017","MTMM","CL")
Bon_M23_FULL=Extract_Sig(Bon_sig_MTMM23_FULL,"DP","2015-2018","MTMM","FULL")
#Bon_M23_IE=Extract_Sig(Bon_sig_MTMM23_IE,"DP","2015-2018","MTMM","IE")
Bon_M23_COM=Extract_Sig(Bon_sig_MTMM23_COM,"DP","2015-2018","MTMM","COM")

MTMMB_Sig=rbind(
Bon_M23_EM,
Bon_M23_CL,
Bon_M23_FULL,
#Bon_M23_IE,
Bon_M23_COM)
install.packages("xlsx")
library(xlsx)
write.xlsx(MTMMB_Sig,"Significant Markers_MTMM_Bon.xlsx", sheetName = "Sheet1", 
  col.names = TRUE, row.names = FALSE, append = FALSE)
```

#### Venn Diagram

``` r
library(VennDiagram)
venn.diagram(
  x = list(
    MTMMB_Sig %>% filter(Covariate=="EM") %>% select(SNP) %>% unlist() , 
    MTMMB_Sig %>% filter(Covariate=="CL") %>% select(SNP) %>% unlist() , 
    MTMMB_Sig %>% filter(Covariate=="FULL") %>% select(SNP) %>% unlist(),
    MTMMB_Sig %>% filter(Covariate=="IE") %>% select(SNP) %>% unlist(),
    MTMMB_Sig %>% filter(Covariate=="COM") %>% select(SNP) %>% unlist()
    ),
  category.names = c("EM" , "CL" , "FULL","IE","COM"),
  filename = 'vennb.png',
  output = TRUE ,
          imagetype="png" ,
          col = "black",
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    alpha = 0.50,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
     1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.05
        )
```

### **Frequently asked questions**

1.  Where can I get WhEATBreeders? WhEATBreeders is currently only
    available on github via Lance F. Merrick’s public repository found
    at <https://github.com/lfmerrick21/WhEATBreeders>.  
      
2.  How to download WhEATBreeders? Please refer the “Getting started”
    section above where I show how to download and install WhEATBreeders
    using devtools.  
      
3.  Where can I get help with issues regarding GWhEAT? All questions
    should be directed to Lance F. Merrick at
    <lance.merrick21@gmail.com>.

### **Further Information**

For more information on individual functions please see the
“Reference\_Manual.pdf” or type ?FUNCITON\_NAME into the R console,
this will pull up specific information of each function inside
WhEATBreeders. For example typing ?manhattan\_plot will pull of the help
page with details about the function that creates the Manhattan plots.

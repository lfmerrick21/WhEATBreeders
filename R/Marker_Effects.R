#just need phenotypic vector and GWAS results

Marker_Effects<-function(Pheno=NULL,GWAS=NULL,alpha=0.05,correction="Bonferonni",messages=TRUE,model="BLINK"){
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(dplyr, compiler)
  #The downloaded link at: http://cran.r-project.org/package=scatterplot3d
  #source("http://zzlab.net/GAPIT/emma.txt")
  #source("http://zzlab.net/GAPIT/gapit_functions.txt")
  #use GAPIT3 development package instead of source code
  library(GAPIT3)
  #Input your phenotypic vector
  myY_train <- Pheno
  if(messages==TRUE){print(paste0("Marker effects are being calculated and using the ",correction," correction."))}
  #Correction Method
  #Bonferonni Correction
  #FDR Correction
  if(model=="BLINK"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1$Rs=GWAS1$effect
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}

  if(model=="FarmCPU"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1$Rs=GWAS1$effect
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}

  if(model=="SUPER"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1=GWAS1[,!names(GWAS$GWAS) %in% c("effect")]
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}

  if(model=="MLM"){
    if(correction=="FDR"){
      GWAS1=GWAS$GWAS
      GWAS1=GWAS1[,!names(GWAS$GWAS) %in% c("effect")]
      GWASFDR=GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(GWAS1)
      sig_mar <- which(GWASFDR$PWIP$`FDR_Adjusted_P-values` <= alpha)
      sig_markers <-as.numeric(c(rownames(GWASFDR$PWIP)[sig_mar]))
    }else{
      sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))
    }}



  #Regular P-value
  if(correction=="P-value"){sig_markers <- which(GWAS$GWAS$P.value <= alpha)}
  PC_mat <-GWAS$PCA[,-1]
  MAS_mat<-GWAS$GD[,-1]
  MAS_mat<-data.frame(MAS_mat[,sig_markers])
  if(length(sig_markers)==1){names(MAS_mat)<-GWAS$GWAS[sig_markers,]$SNP[[1]]}
  if(length(sig_markers)==0){print("Your GWAS is terrible and Halle thinks you suck.")}
  if(length(sig_markers)==0){stop("Your GWAS and correction method found zero significant markers.")}
  if(messages==TRUE){print(paste0("There are ",length(sig_markers)," significant marker(s)."))}
  GLM_data <- data.frame(myY_train,PC_mat,MAS_mat)
  names(GLM_data)[1] <- "Y"
  MAS_model <- lm(Y ~ ., data = GLM_data)
  Marker_Results=data.frame()
  cov_length=ncol(data.frame(myY_train))+ncol(PC_mat)
  for(i in cov_length+1:ncol(as.data.frame(MAS_mat)) ){
    #Null Model
    GLM_data_null <- data.frame(myY_train,PC_mat)
    names(GLM_data_null)[1] <- "Y"
    MAS_model_null <- lm(Y ~ ., data = GLM_data_null)
    R2null=summary(MAS_model_null)$r.squared

    #Model with 1 marker
    SNP_Name=names(GLM_data)[i]
    GLM_data_w1 <- data.frame(myY_train,PC_mat,GLM_data[,i])
    names(GLM_data_w1)[1] <- "Y"
    names(GLM_data_w1)[cov_length+1] <- SNP_Name
    MAS_model_w1 <- lm(Y ~ ., data = GLM_data_w1)
    R2w1=summary(MAS_model_w1)$r.squared
    R2w1m=R2w1-R2null #r2 for significant marker

    #Model without marker
    GLM_data_sm1 <- GLM_data[,-i]
    MAS_model_sm1 <- lm(Y ~ ., data = GLM_data_sm1)
    R2full=summary(MAS_model)$r.squared
    R2wo=summary(MAS_model_sm1)$r.squared
    R2wo1m=summary(MAS_model)$r.squared-summary(MAS_model_sm1)$r.squared #r2 for significant marker
    Full_Effects=MAS_model$coefficients[i]
    SNP_Effects=MAS_model_w1$coefficients[cov_length+1]

    results=data.frame(SNP_Name,R2null,R2w1,R2w1m,R2full,R2wo,R2wo1m,Full_Effects,SNP_Effects)
    names(results)=c("SNP","R2.null.model","R2.add.single.marker","R2.of.single.marker","R2.full.model","R2.minus.single.marker","R2.diff","Effects.marker.full.model","Effects.single.marker")
    Marker_Results=rbind(Marker_Results,results)

    if(correction=="FDR"){
      GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect","Rs")],Marker_Results,by="SNP")
      GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
    }else{
      GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
    }

    if(model=="BLINK"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect","Rs")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}

    if(model=="FarmCPU"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect","Rs")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}

    if(model=="SUPER"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}

    if(model=="MLM"){
      if(correction=="FDR"){
        GWAS_Table=left_join(GWASFDR$PWIP[sig_mar,!names(GWASFDR$PWIP) %in% c("effect")],Marker_Results,by="SNP")
        GWAS_Table <- GWAS_Table[order(GWAS_Table$SNP),]
      }else{
        GWAS_Table=left_join(GWAS$GWAS[sig_markers,!names(GWAS$GWAS) %in% c("effect")],Marker_Results,by="SNP")
      }}




  }
  if(messages==TRUE){print(paste0("Wow! Your signficant markers and covariates accounted for ",round(R2full*100,2),"% of the variation."))}
  return(GWAS_Table)
}

Marker_Effects_Allele<-function(Pheno=NULL,GWAS=NULL,alpha=0.05,correction="Bonferonni",messages=TRUE,Markers=NULL,model="BLINK"){

  #Pheno=GBS_qam_adj22$pheno[,c(22)];GWAS=myGAPIT_MLM_qam_adj_22;alpha=0.05;correction="FDR";messages=TRUE;Markers=EM_GBS_SNPs;model="MLM"
  #Pheno=GBS_qam_adj24$pheno[,c(24)];GWAS=FarmCPU_qam_adj_24;alpha=0.05;correction="Bonferonni";messages=TRUE;Markers=EM_GBS_SNPs



  corr_sig=Marker_Effects(Pheno=Pheno,GWAS=GWAS,alpha=alpha,correction=correction,messages=messages,model=model)
  names(Markers)[1]<-"SNP"
  GWAS_Results=left_join(corr_sig,Markers[,1:3],by="SNP")

  if(model=="BLINK"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,16,17,4:15)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,15,16,4:14)]
    }}

  if(model=="FarmCPU"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,16,17,4:15)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,15,16,4:14)]
    }}

  if(model=="SUPER"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,17,18,4:16)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,17,18,4:16)]
    }}

  if(model=="MLM"){
    if(correction=="FDR"){
      GWAS_Results=GWAS_Results[,c(1:3,18,19,4:6,9:17)]
    }else{
      GWAS_Results=GWAS_Results[,c(1:3,17,18,4:6,9:16)]
    }}
  return(GWAS_Results)
}

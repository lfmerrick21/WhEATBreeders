#Sommer_MM_GWAS
sommer_MTMM=function(Y,SNP_INFO,model,X,A,MAFX=NULL){
  colnames(SNP_INFO)<-c('SNP','Chr','Pos')
  ## preparing genotype and Kinship data
  Y_<-Y
  X<-X[which(rownames(X)%in%Y_[,1]),]
  Y<-Y_[which(Y_[,1]%in%rownames(X)),]
  ### enter here filtering step for 2 different kinships for different subsets .

  Y1<-(Y[,2])
  Y2<-(Y[,3])
  names(Y1)<-Y[,1]
  names(Y2)<-Y[,1]
  Y1_<-na.omit(Y1)
  Y2_<-na.omit(Y2)

  ecot_id1<-as.integer(names(Y1_))
  ecot_id2<-as.integer(names(Y2_))

  K<-as.matrix(A[which(rownames(A)%in%names(Y1)),which(colnames(A)%in%names(Y1))])
  K1<-as.matrix(A[which(rownames(A)%in%ecot_id1),which(colnames(A)%in%ecot_id1)])
  K2<-as.matrix(A[which(rownames(A)%in%ecot_id2),which(colnames(A)%in%ecot_id2)])

  n<-nrow(Y)
  n1<-length(ecot_id1)
  n2<-length(ecot_id2)


  # combining the traits

  Y_ok<-c(Y1,Y2)

  #environment

  Env<-c(rep(0,n),rep(1,n))

  #standardize kinship

  K_stand<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
  K_inv<<-solve(K_stand)
  #K_stand1<-(n1-1)/sum((diag(n1)-matrix(1,n1,n1)/n1)*K1)*K1
  #K_stand2<-(n2-1)/sum((diag(n2)-matrix(1,n2,n2)/n2)*K2)*K2

  count2<-function(x) {length(which(x==2))}
  count1<-function(x) {length(which(x==1))}

  AC_2 <- data.frame(colnames(X),apply(X,2,count2))
  colnames(AC_2)<-c('SNP','AC_2')
  AC_1 <- data.frame(colnames(X),apply(X,2,count1))
  colnames(AC_1)<-c('SNP','AC_1')

  MAF_1<-data.frame(AC_2,AC_1$AC_1,AC_0=nrow(X)-AC_1$AC_1-AC_2$AC_2)

  MAF_2<-data.frame(MAF_1,MAC=apply(MAF_1[,c(2,4)],1,min))

  MAF_3<-data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X)))

  MAF_ok<<-merge(SNP_INFO,MAF_3,by='SNP')

  rm(AC_1,MAF_1,MAF_2,MAF_3)

  #Filter for MAF




  if(!is.null(MAFX)){
    MAF<-subset(MAF_ok,MAF==0)[,1]
    X_ok<-X[,!colnames(X) %in% MAF]
  }else{
    X_ok<-X
  }

  Xo<-rep(1,2*n)
  Xo1<-rep(1,n1)
  ex1<-as.matrix(Xo1)
  Xo2<-rep(1,n2)
  ex2<-as.matrix(Xo2)


  varcov<-model$sigma[,,1]
  ve<-model$sigma[,,2]

  rho_g<-varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2])
  rho_e<-ve[1,2]/sqrt(ve[1,1]*ve[2,2])
  rho_p<-(varcov[1,2]+ve[1,2])/sqrt((varcov[1,1]+ve[1,1])*(varcov[2,2]+ve[2,2]))
  pearson<-cor(Y1,Y2)
  heriti1<-varcov[1,1]/(varcov[1,1]+ve[1,1])
  heriti2<-varcov[2,2]/(varcov[2,2]+ve[2,2])
  correlation<-list(pearson=pearson,gen_cor=rho_g,env_cor=rho_e,phen_cor=rho_p,h1_joint=heriti1,h2_joint=heriti2,converge=model$convergence)

  K_comb<-kronecker(varcov,K_stand)
  rownames(K_comb)<-c(ecot_id1,ecot_id2)
  colnames(K_comb)<-c(ecot_id1,ecot_id2)

  I_comb<-kronecker(ve,diag(n))
  rownames(I_comb)<-rownames(K_comb)
  colnames(I_comb)<-colnames(K_comb)

  bigK<-K_comb+I_comb


  M<-solve(chol(bigK))


  # scaling of the SNPs to interpret GxE results
  X_ok1<-X_ok*sqrt(varcov[1,1])
  X_ok2<-X_ok*sqrt(varcov[2,2])

  Y_t<-crossprod(M,Y_ok)
  cof_t<-crossprod(M,cbind(Xo,Env))

  RSS_env<-sum(lsfit(cof_t,Y_t,intercept = FALSE)$residuals^2)

  nbchunks<-5
  m<-ncol(X_ok)

  RSS_full<-list()
  RSS_glob<-list()

  for (j in 1:(nbchunks-1)){
    X_<-rbind(X_ok1[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))],X_ok2[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    X_t<-array(dim=c(nrow(X_),ncol(X_),2))
    X_t[,,1]<-crossprod(M,X_)
    X_t[,,2]<-crossprod(M,(Env*X_))
    RSS_full[[j]]<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
    RSS_glob[[j]]<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
    rm(X_,X_t)}
  X_<-rbind(X_ok1[,((j)*round(m/nbchunks)+1):m],X_ok2[,((j)*round(m/nbchunks)+1):m])
  X_t<-array(dim=c(nrow(X_),ncol(X_),2))
  X_t[,,1]<-crossprod(M,X_)
  X_t[,,2]<-crossprod(M,(Env*X_))
  RSS_full[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
  RSS_glob[[nbchunks]]<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
  rm(X_,X_t,j)

  RSS_full<-unlist(RSS_full)
  RSS_glob<-unlist(RSS_glob)

  #nb parameters in each models
  #env
  par_env<-ncol(cof_t)
  par_glob<-par_env+1
  par_full<-par_glob+1

  ##FTESTS
  #FULL vs NULL
  F_full<-(rep(RSS_env,m)/RSS_full-1)*(2*n-par_full)/(par_full-par_env)
  F_ge<-(RSS_glob/RSS_full-1)*(2*n-par_full)/(par_full-par_glob)
  F_glob<-(rep(RSS_env,m)/RSS_glob-1)*(2*n-par_glob)/(par_glob-par_env)

  pval_full<-pf(F_full,par_full-par_env,2*n-par_full,lower.tail=FALSE)
  pval_ge<-pf(F_ge,par_full-par_glob,2*n-par_full,lower.tail=FALSE)
  pval_glob<-pf(F_glob,par_glob-par_env,2*n-par_glob,lower.tail=FALSE)

  #outputs
  outfull<<-data.frame('SNP'=colnames(X_ok),pval=pval_full)
  outge<<-data.frame('SNP'=colnames(X_ok),pval=pval_ge)
  outglob<<-data.frame('SNP'=colnames(X_ok),pval=pval_glob)



  sommer_scores=data.frame(t(model$scores))
  py=sommer_scores[,3]
  if(max(py,na.rm=TRUE)>1){
    pvaly<-10^-py
  }else{
    pvaly<-py
  }
  sommer_scores$Y.score=pvaly
  pc=sommer_scores[,4]
  if(max(pc,na.rm=TRUE)>1){
    pvalc<-10^-pc
  }else{
    pvalc<-pc
  }
  sommer_scores$DP_BLUP.score=pvalc
  sommer_scores=cbind(SNP=rownames(sommer_scores),sommer_scores)

  data.out<-merge(MAF_ok,cbind(outfull,outge[,2],outglob[,2]),by='SNP')
  data.outt<-left_join(data.out,sommer_scores,by='SNP')
  colnames(data.outt)[9:11]<-c('pval_full','pval_trait_specific','pval_trait_common')
  data.out_<-data.outt[order(data.outt[,3]),]
  out<-data.out_[order(data.out_[,2]),]
  results<-list(phenotype=Y,pvals=out,statistics=correlation,kinship=K)
  return(results)
}

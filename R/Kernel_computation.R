#Kernel compuation from Montesinos-Lopez
Kernel_computation=function(X,name, degree, nL){
  p=ncol(X)
  x1=X
  x2=X
  d=degree
  ############Polynomial kernel##################
  K.Polynomial=function(x1, x2=x1, gamma=1,
                        b=0, d=3)
  { (gamma*(as.matrix(x1)%*%t(x2))+b)^d}
  ############Sigmoid kernel####################
  K.Sigmoid=function(x1,x2=x1, gamma=1, b=0)
  { tanh(gamma*(as.matrix(x1)%*%t(x2))+b) }
  ############Gaussian kernel##################
  l2norm=function(x){sqrt(sum(x^2))}
  K.Gaussian=function(x1,x2=x1, gamma=1){
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol
                     (x2<- t(x2)),
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j])^2)))}
  ##########Arc-cosine kernel with 1 hidden layer
  K.AK1_Final<-function(x1,x2){
    n1<-nrow(x1)
    n2<-nrow(x2)
    x1tx2<-x1%*%t(x2)
    norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))
    norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))
    costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)
    costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1
    theta<-acos(costheta)
    normx1x2<-norm1%*%t(norm2)
    J = (sin(theta)+(pi-theta)*cos(theta))
    AK1 = 1/pi*normx1x2*J
    AK1<-AK1/median(AK1)
    colnames(AK1)<-rownames(x2)
    rownames(AK1)<-rownames(x1)
    return(AK1)
  }
  ####Kernel Arc-Cosine with deep=L#########
  AK_L_Final<-function(AK1,nL){
    n1<-nrow(AK1)
    n2<-ncol(AK1)
    AKl1 = AK1
    for (l in 1:nL){
      AKAK<-tcrossprod(diag(AKl1),diag(AKl1))
      costheta<-AKl1*(AKAK^(-1/2))
      costheta[which(costheta>1,arr.ind = TRUE)] = 1
      theta<-acos(costheta)
      AKl<-(1/pi)*(AKAK^(1/2))*(sin(theta)+(pi-theta)*cos
                                (theta))
      AKl1 = AKl
    }
    AKl<-AKl/median(AKl)
    rownames(AKl)<-rownames(AK1)
    colnames(AKl)<-colnames(AK1)
    return(AKl)
  }
  ########Exponencial Kernel############
  K.exponential=function(x1,x2=x1, gamma=1){
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol
                     (x2<- t(x2)),
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j]))))}
  if (name=="Linear") {
    K=X%*%t(X)/p
  } else if (name=="Polynomial") {
    K=K.Polynomial(x1=x1, x2=x1, gamma=1/p,
                   b=0, d=d)
  } else if (name=="Sigmoid") {
    K=K.Sigmoid(x1=x1, x2=x1, gamma=1/p, b=0)
  }else if (name=="Gaussian") {
    K=K.Gaussian(x1=x1, x2=x1, gamma=1/p)
  } else if (name=="AK1") {
    K= K.AK1_Final(x1=x1, x2=x1)
  } else if (name=="AKL") {
    AK1=K.AK1_Final(x1=x1, x2=x1)
    K=AK_L_Final(AK1=AK1,nL=nL)
  } else {
    K=K.exponential(x1=x1,x2=x1,gamma=1/p)
  }
}

Sparse_Kernel_Construction=function(m,X,name,degree,nL){
  degree=degree
  nl=nL
  m=m
  XF=X
  p=ncol(XF)
  pos_m=sample(1:nrow(XF),m)
  X_m=XF[pos_m,]
  dim(X_m)
  ########Gaussian Kernel function############
  l2norm=function(x){sqrt(sum(x^2))}
  K.radial=function(x1,x2=x1, gamma=1){
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol
                     (x2<- t(x2)),
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j])^2)))}
  ########Polynomial Kernel############
  K.polynomial=function(x1, x2=x1, gamma=1, b=0,
                        d=degree)
  { (gamma*(as.matrix(x1)%*%t(x2))+b)^d}
  ########Sigmoid Kernel############
  K.sigmoid=function(x1,x2=x1, gamma=1, b=0)
  { tanh(gamma*(as.matrix(x1)%*%t(x2))+b) }
  ########Exponencial Kernel############
  K.exponential=function(x1,x2=x1, gamma=1){
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol
                     (x2<- t(x2)),
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j]))))}
  ############Arcosine kernel with deep=1#######
  K.AK1_Final<-function(x1,x2){
    n1<-nrow(x1)
    n2<-nrow(x2)
    x1tx2<-x1%*%t(x2)
    norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))
    norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))
    costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)
    costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1
    theta<-acos(costheta)
    normx1x2<-norm1%*%t(norm2)
    J = (sin(theta)+(pi-theta)*cos(theta))
    AK1 = 1/pi*normx1x2*J
    AK1<-AK1/median(AK1)
    colnames(AK1)<-rownames(x2)
    rownames(AK1)<-rownames(x1)
    return(AK1)
  }
  ####Kernel Arc-Cosine with deep=L#####
  diagAK_f<-function(dAK1)
  {
    AKAK = dAK1^2
    costheta = dAK1*AKAK^(-1/2)
    costheta[which(costheta>1,arr.ind = TRUE)] = 1
    theta = acos(costheta)
    AKl = (1/pi)*(AKAK^(1/2))*(sin(theta)+(pi-theta)*cos
                               (theta))
    AKl
    AKl<-AKl/median(AKl)
  }
  AK_L_Final<-function(AK1,dAK1,nl){
    n1<-nrow(AK1)
    n2<-ncol(AK1)
    AKl1 = AK1
    for (l in 1:nl){
      AKAK<-tcrossprod(dAK1,diag(AKl1))
      costheta<-AKl1*(AKAK^(-1/2))
      costheta[which(costheta>1,arr.ind = TRUE)] = 1
      theta<-acos(costheta)
      AKl<-(1/pi)*(AKAK^(1/2))*(sin(theta)+(pi-theta)*cos
                                (theta))
      dAKl = diagAK_f(dAK1)
      AKl1 = AKl
      dAK1 = dAKl
    }
    AKl<-AKl/median(AKl)
    rownames(AKl)<-rownames(AK1)
    colnames(AKl)<-colnames(AK1)
    return(AKl)
  }
  AK_ALL=K.AK1_Final(x1=XF,x2=XF)
  AK11=K.AK1_Final(x1=XF,x2=X_m)
  AKL2=AK_L_Final(AK1=AK11,dAK1=diag
                  (AK_ALL),nl=5)
  dim(AKL2)
  ####Step 1 compute K_m###############
  if (name=="Linear") {
    K_m=X_m%*%t(X_m)/p
    ######Step 2 comput K_n_m###########
    K_n_m=XF%*%t(X_m)/p
  } else if (name=="Polynomial") {
    K_m=K.polynomial(x1=X_m,x2=X_m,gamma=1/p)
    ######Step 2 comput K_n_m###########
    K_n_m=K.polynomial(x1=XF,x2=X_m,gamma=1/p)
  } else if (name=="Sigmoid") {
    K_m=K.sigmoid(x1=X_m,x2=X_m,gamma=1/p)
    ######Step 2 comput K_n_m###########
    K_n_m=K.sigmoid(x1=XF,x2=X_m,gamma=1/p)
  }else if (name=="Gaussian") {
    K_m=K.radial(x1=X_m,x2=X_m,gamma=1/p)
    ######Step 2 comput K_n_m###########
    K_n_m=K.radial(x1=XF,x2=X_m,gamma=1/p)
  } else if (name=="AK") {
    K_m=K.AK1_Final(x1=X_m,x2=X_m)
    ######Step 2 comput K_n_m###########
    K_nm1=K.AK1_Final(x1=XF,x2=X_m)
    K_all=K.AK1_Final(x1=XF,x2=XF)
    K_n_m=AK_L_Final(AK1=K_nm1,dAK1=diag
                     (K_all),nl=nl)
  } else {
    K_m=K.exponential(x1=X_m,x2=X_m,gamma=1/p)
    ######Step 2 comput K_n_m###########
    K_n_m=K.exponential(x1=XF,x2=X_m,gamma=1/p)
  }
  ######Step 3 compute Eigen value decomposition of
  K_m######
  EVD_K_m=eigen(K_m)
  ####Eigenvectors
  U=EVD_K_m$vectors
  ###Eigenvalues###
  S=EVD_K_m$values
  S[which(S<0)]=0
  ####Square root of the inverse of eigenvelues #####
  S_0.5_Inv=sqrt(1/S)
  #####Diagonal matrix of square root of inverse of ingenvalues###
  S_mat_Inv=diag(S_0.5_Inv)
  ######Computing matrix P
  P=K_n_m%*%U%*%S_mat_Inv
  return(P)}

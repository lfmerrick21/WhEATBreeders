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
  ########Exponential Kernel############
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


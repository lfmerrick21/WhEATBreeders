#Cullis Heritability (From Shantel Martinez Github)

Cullis_H2=function(model,Line){
  library(arm)
  ses<- se.ranef(model)$Line #where 'm' is your model object from 'lmer' (replace 'genotypes' with whatever you call your individuals in the data)
  v_BLUP<- ses^2
  sigma2_g=VarCorr(model, comp="Variance")$Line[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g)  #where sigma2_g is the genetic variance estimated with the model saved in 'm'
  H2<- round(mean(Reliability),3) #This is equivalent to broad-sense heritability on the line-mean (or family-mean, if your individuals are non-inbred families) basis
  H2
}

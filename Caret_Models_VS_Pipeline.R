if(Package=="caret"){
  Results=sapply(1:Replications, function(i,...){Results=Caret_Models_CV(train_genotypes = get(paste0("GBS_2_",Training,"_",Trait))$geno,
                                                                         train_phenotype = get(paste0("GBS_2_",Training,"_",Trait))$pheno[,c(1,3)],
                                                                         test_genotypes= get(paste0("GBS_2_",Prediction,"_",Trait))$geno,
                                                                         test_phenotype= get(paste0("GBS_2_",Prediction,"_",Trait))$pheno[,c(1,3)],
                                                                         type=type,
                                                                         model=model,
                                                                         Kernel=Kernel,
                                                                         markers=markers,
                                                                         Sparse=Sparse,
                                                                         m=m,
                                                                         degree=degree,
                                                                         nL=nL,
                                                                         transformation=transformation,
                                                                         sampling="up",
                                                                         repeats=5,
                                                                         method="repeatedcv"
  )})
}

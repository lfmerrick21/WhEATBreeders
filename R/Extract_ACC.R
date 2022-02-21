Extract_ACC=function(Results,Replications,Training,Model,Kernel,CV_Message,Trait){
  model_vect <- c("MAE","Pearson","R2","RMSE","Spearman")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(Results[1,])){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(Results[1,i]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  colnames(BGLR_acc_table)<-c("Rep","Metric","Value")
  en=rep(Training,Replications)
  mo=rep(Model,Replications)
  ke=rep(Kernel,Replications)
  if(!is.null(CV_Message)){cvm=rep("CV",Replications)}
  tr=rep(Trait,Replications)
  if(!is.null(CV_Message)){
    ac=data.frame(Environment=en,Model=mo,Kernel=ke,CV=cvm,Trait=tr,BGLR_acc_table)
  }else{
    ac=data.frame(Environment=en,Model=mo,Kernel=ke,Trait=tr,BGLR_acc_table)
  }

  return(ac)
}

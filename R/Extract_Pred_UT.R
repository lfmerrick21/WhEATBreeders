Extract_Pred_UT=function(Results,Replications,Training,Model,Kernel,CV_Message,Trait){
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(Results)){
    results_long <- data.frame(Rep=rep(i, nrow(Results[i][[1]])),Results[i][[1]])
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  en=Training
  mo=Model
  ke=Kernel
  if(!is.null(CV_Message)){cvm="CV"}
  tr=Trait
  if(!is.null(CV_Message)){
    pred=data.frame(Environment=en,Model=mo,Kernel=ke,CV=cvm,Trait=tr,BGLR_acc_table)
  }else{
    pred=data.frame(Environment=en,Model=mo,Kernel=ke,Trait=tr,BGLR_acc_table)
  }

  return(pred)
}

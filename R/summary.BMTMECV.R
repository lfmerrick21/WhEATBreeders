summary.BMTMECV <- function (results, information = 'compact', digits = 4, ...) {
  library(stringr)
  library(rrBLUP)
  library(BGLR)
  library(caret)
  library(Metrics)
  library(mpath)

  library(lsa)
  library(keras)
  library(BMTME)
  #library(plyr)
  #library(tidyr)
  library (dplyr)
  results %>%
    group_by(Environment, Trait, Partition) %>%
    summarise (MSE = mean((Predicted-Observed)^2, na.rm = T),
               MAAPE = mean(atan(abs(Observed-Predicted)/abs(Observed)), na.rm = T),
               ACC=cor(Predicted,Observed, use = "pairwise.complete"),
               SACC=cor(Predicted,Observed, use = "pairwise.complete", method = c("spearman")),
               Rsquare=R2(pred=Predicted,obs=Observed,na.rm=TRUE),
               RMSE=RMSE(pred=Predicted,obs=Observed,na.rm=TRUE),
               MAE=MAE(pred=Predicted,obs=Observed,na.rm=TRUE)) %>%
    select(Environment, Trait, Partition, MSE, MAAPE,ACC,SACC,Rsquare,RMSE,MAE) %>%
    mutate_if(is.numeric, funs (round(., digits))) %>%
    as.data.frame()->presum
  presum %>% group_by(Environment, Trait) %>%
    summarise(SE_MAAPE = sd(MAAPE, na.rm = T)/sqrt(n()),
              MAAPE = mean(MAAPE, na.rm = T),
              SE_MSE = sd(MSE, na.rm = T)/sqrt(n()),
              MSE = mean(MSE, na.rm =T),
              SD=sd(ACC,na.rm=T),
              ACC = mean(ACC, na.rm =T),
              SD_SACC=sd(SACC,na.rm=T),
              SACC = mean(SACC, na.rm =T),
              SD_Rsquare=sd(Rsquare,na.rm=T),
              Rsquare=mean(Rsquare,na.rm=TRUE),
              SD_RMSE=sd(RMSE,na.rm=T),
              RMSE=mean(RMSE,na.rm=TRUE),
              SD_MAE=sd(MAE,na.rm=T),
              MAE=mean(MAE,na.rm=TRUE)) %>%
    dplyr::select(Environment, Trait, MSE, SE_MSE, MAAPE,SE_MAAPE,ACC,SD,Rsquare,SD_Rsquare,RMSE,SD_RMSE,MAE,SD_MAE) %>%
    mutate_if(is.numeric, funs (round(., digits))) %>%
    as.data.frame() -> finalSum
  out <-switch (information,
                compact = finalSum,
                complete = presum,
                extended={
                  finalSum$Partition <- 'All'
                  presum$Partition <- as.character(presum$Partition)
                  presum$SE_MSE < - NA
                  presum$SE_MAAPE <- NA
                  presum$SD <- NA
                  rbind (presum, finalSum)
                }
  )
  return (out)
}

MLP_CV<-function(matrix,trait=c("IT","SEV"),model="ST",digits=4,nCVI=5,K=5,folds=5,Stage=NULL,UN=TRUE){


  if(model=="ST"){
    Phenotype=matrix$MTME$YUN
    #Y <- matrix$BMTME$Y
    #Y[-fold_indices,] <- NA
    #Z_L=matrix$MTME$Z_L
    #Z_E=matrix$MTME$Z_E
    K_expanded=matrix$MTME$K_expanded
    #K_GE=matrix$MTME$K_GE
    #Phenotype$Genotype=as.factor(Phenotype$Genotype)
    #Phenotype$Env=as.factor(Phenotype$Env)
    Y2=Phenotype[complete.cases(Phenotype),]
    X2=K_expanded[complete.cases(Phenotype),]

    results = data.frame()
    for(j in 1:length(trait)){
      ############Selecting the response variable#######################

      #Y=as.matrix(Y2[,-c(1, 2)])
      ####Training testing sets using the BMTMEpackage###############
      pheno=data.frame(Line = Y2[, 1],
                       Env=Y2[, 2],
                       Response=Y2[,trait[j]])
      pheno$Line=as.character(pheno$Line)
      pheno$Env=as.character(pheno$Env)
      CrossV=CV.KFold(pheno, DataSetID = 'Line', K=K)
      #fold_list <- make_CV_sets(length(phenotypes), k = folds)
      #######Final X and y for fitting the model###################
      y=Y2[,trait[j]]
      X=X2
      #########Grid of hyperparameters######## ############# #####################.
      if(is.null(Stage)){
        Stage=expand.grid(units_M=c(16,32,64),epochs_M = 1000, Drop_per = c(0.0,0.05,0.15),
                          alpha=c(0.001,0.01,0.4))
      }
      ############Outer Cross-validation#######################
      digits = digits
      Names_Traits = trait[j]

      t = 1

      if(trait[j]=="IT"){
        y = replace(y, y < 0, 0)
        y = replace(y, y > 9, 9)
      }
      if(trait[j]=="SEV"){
        y = replace(y, y < 0, 0)
        y = replace(y, y > 100, 100)
      }

      for (o in 1:K){
        tst_set = CrossV$CrossValidation_list[[o]]
        #tst_set = which(fold_list[[o]])
        length(tst_set)
        X_trn = (X[-tst_set,])
        X_tst = (X[tst_set,])
        y_trn = sqrt((y[-tst_set]))
        y_tst = sqrt((y [tst_set]))
        nCVI = nCVI ####Number of folds for inner CV
        #i = 1
        #########Matrices for saving the output of inner CV#######################.
        Tab_pred_MSE = matrix (NA,ncol = length(Stage[,1]),
                               nrow = nCVI)
        Tab_pred_Epoch = matrix (NA,ncol = length (Stage[,1]),
                                 nrow = nCVI)
        Tab_pred_Units = matrix (NA,ncol = length (Stage[,1]),
                                 nrow = nCVI)
        Tab_pred_Drop = matrix (NA,ncol = length (Stage[,1]),
                                nrow = nCVI)
        Tab_pred_Alpha = matrix (NA,ncol = length (Stage[,1]),
                                 nrow = nCVI)
        X_trI = X_trn
        y_trI = y_trn
        for(i in 1:nCVI){
          for(stage in seq_len(dim(Stage)[1])){
            X_trII = X_trI
            y_trII = y_trI
            units_M <- Stage[stage, 1]
            epochs_M <- Stage[stage, 2]
            Drop_per = Stage[stage, 3]
            Alpha=Stage[stage, 4]
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per) %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per) %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per) %>%
                layer_dense(units = 1)

              model<-keras_model(input,output) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}

            model<-build_model()
            model %>% summary()
            print_dot_callback <- callback_lambda(
              on_epoch_end = function (epoch, logs) {
                if (epoch %% 20 == 0) cat("\n")
                cat(".")
              })
            ########Fitting the model with Early stopping#######
            early_stop <- callback_early_stopping(monitor = "val_loss",mode = 'min',patience =50)
            ###########Fit of the model for each values of the grid#################
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              X_trII, y_trII,
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
            ###########Saving the output of each hyperparameter###################
            No.Epoch_Min = length(model_fit_Final$metrics$val_mse)
            Min_MSE = model_fit_Final$metrics$val_mse[No.Epoch_Min]
            Tab_pred_MSE[i,stage] = Min_MSE
            Tab_pred_Units[i,stage] = units_M
            Tab_pred_Epoch[i,stage] = No.Epoch_Min[1]
            Tab_pred_Drop[i,stage] = Drop_per
            Tab_pred_Alpha[i,stage] = Alpha
          }
        }
        ##############Selecting the optimal hyperparameters##############
        Median_MSE_Inner = apply(Tab_pred_MSE,2,median)
        Units_Inner = apply(Tab_pred_Units,2,max)
        Drop_Inner = apply(Tab_pred_Drop,2,max)
        Epoch_Inner = apply(Tab_pred_Epoch,2,median)
        Alpha_Inner = apply(Tab_pred_Alpha,2,median)
        Pos_Min_MSE = which(Median_MSE_Inner==min(Median_MSE_Inner))
        Units_O=Units_Inner[Pos_Min_MSE]
        Epoch_O = Epoch_Inner[Pos_Min_MSE]
        Drop_O = Drop_Inner[Pos_Min_MSE]
        Alpha_O = Alpha_Inner[Pos_Min_MSE]
        ###########Refitting the model with the optimal values#################
        input_Sec<- layer_input(shape=dim(X_trn)[2],name="covars")


        output_Sec<-input_Sec %>%
          layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha_O)) %>%
          layer_dropout(rate = Drop_O) %>%
          layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha_O)) %>%
          layer_dropout(rate = Drop_O) %>%
          layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha_O)) %>%
          layer_dropout(rate = Drop_O) %>%
          layer_dense(units =1)

        model_Sec=keras_model(input_Sec,output_Sec) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          X_trn, y_trn,
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
        ####Prediction of testing set##########################.
        Yhat = model_Sec %>% predict(X_tst)
        y_p = Yhat
        y_p_tst = as.numeric(y_p)
        ###########Saving the predicctions of each outer testing set#################
        results=rbind(results,data.frame(Position = tst_set,
                                         Trait = trait[j],
                                         Partition = o,
                                         Genotype = Y2$Genotype[tst_set],
                                         Environment = Y2$ENV[tst_set],
                                         Units = Units_O,
                                         Epochs = Epoch_O,
                                         Drop_Out = Drop_O,
                                         Alpha_Out = Alpha_O,
                                         Observed = round(y_tst, digits), #$response, digits),
                                         Predicted = round(y_p_tst, digits)
        ))
        cat("CV = ",o,"\n")
      }
    }
    Pred_Summary=summary.BMTMECV(results,information = 'compact',digits = digits)
    Pred_Summary1=summary.BMTMECV(results,information = 'complete',digits = digits)
    Results=list(Results=results,Compact=Pred_Summary,Complete=Pred_Summary1)
    #############Average of prediction performance##################################
    return(Results)
  }

  if(model=="MT"){
    Phenotype=matrix$MTME$YUN
    #Y <- matrix$BMTME$Y
    #Y[-fold_indices,] <- NA
    #Z_L=matrix$MTME$Z_L
    #Z_E=matrix$MTME$Z_E
    K_expanded=matrix$MTME$K_expanded
    #K_GE=matrix$MTME$K_GE
    #Phenotype$Genotype=as.factor(Phenotype$Genotype)
    #Phenotype$Env=as.factor(Phenotype$Env)
    Y2=Phenotype[complete.cases(Phenotype),]
    X2=K_expanded[complete.cases(Phenotype),]

    #Y2=Y2[1:50,]
    #X2=X2[1:50,]
    ############Selecting the response variable#######################

    #Y=as.matrix(Y2[,-c(1, 2)])
    ####Training testing sets using the BMTMEpackage###############
    pheno=data.frame(Line = Y2[, 1],
                     Env=Y2[, 2],
                     Response=Y2[,trait[1]])
    pheno$Line=as.character(pheno$Line)
    pheno$Env=as.character(pheno$Env)
    CrossV=CV.KFold(pheno, DataSetID = 'Line', K=K)
    #fold_list <- make_CV_sets(length(phenotypes), k = folds)
    #######Final X and y for fitting the model###################
    y=Y2[,trait]
    X=X2
    dim(X)
    dim(y)
    #########Grid of hyperparameters######## ############# #####################.
    if(is.null(Stage)){
      Stage=expand.grid(units_M=c(16,32,64),epochs_M = 1000, Drop_per = c(0.0,0.05,0.15),
                        alpha=c(0.001,0.01,0.4))
    }
    ############Outer Cross-validation#######################
    digits = digits
    Names_Traits = trait
    results = data.frame()
    t = 1

    y$IT=replace(y$IT, y$IT < 0, 0)
    y$IT=replace(y$IT, y$IT > 9, 9)
    y$SEV=replace(y$SEV, y$SEV < 0, 0)
    y$SEV=replace(y$SEV, y$SEV > 10, 10)


    for (o in 1:K){
      tst_set = CrossV$CrossValidation_list[[o]]
      #tst_set = which(fold_list[[o]])
      length(tst_set)
      X_trn = (X[-tst_set,])
      X_tst = (X[tst_set,])
      y_trn = sqrt((y[-tst_set,]))
      y_tst = sqrt((y[tst_set,]))
      nCVI = nCVI ####Number of folds for inner CV
      #i = 1
      #########Matrices for saving the output of inner CV#######################.
      Tab_pred_MSE = matrix (NA,ncol = length(Stage[,1]),
                             nrow = nCVI)
      Tab_pred_Epoch = matrix (NA,ncol = length (Stage[,1]),
                               nrow = nCVI)
      Tab_pred_Units = matrix (NA,ncol = length (Stage[,1]),
                               nrow = nCVI)
      Tab_pred_Drop = matrix (NA,ncol = length (Stage[,1]),
                              nrow = nCVI)
      Tab_pred_Alpha = matrix (NA,ncol = length (Stage[,1]),
                               nrow = nCVI)
      X_trI = X_trn
      y_trI = y_trn
      for(i in 1:nCVI){
        for(stage in seq_len(dim(Stage)[1])){
          X_trII = X_trI
          y_trII = y_trI
          units_M <- Stage[stage, 1]
          epochs_M <- Stage[stage, 2]
          Drop_per = Stage[stage, 3]
          Alpha = Stage[stage, 4]

          if(length(trait)==2){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")

              model<-keras_model(input,list(yhat1,yhat2)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }

          if(length(trait)==3){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")
              # add output 2
              yhat3 <- output %>%
                layer_dense(units = 1, name="yhat3")

              model<-keras_model(input,list(yhat1,yhat2,yhat3)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }

          if(length(trait)==4){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")
              # add output 2
              yhat3 <- output %>%
                layer_dense(units = 1, name="yhat3")
              # add output 2
              yhat4 <- output %>%
                layer_dense(units = 1, name="yhat4")

              model<-keras_model(input,list(yhat1,yhat2,yhat3,yhat4)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }

          if(length(trait)==5){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")
              # add output 2
              yhat3 <- output %>%
                layer_dense(units = 1, name="yhat3")
              # add output 2
              yhat4 <- output %>%
                layer_dense(units = 1, name="yhat4")
              # add output 2
              yhat5 <- output %>%
                layer_dense(units = 1, name="yhat5")

              model<-keras_model(input,list(yhat1,yhat2,yhat3,yhat4,yhat5)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }


          model<-build_model()
          model %>% summary()
          print_dot_callback <- callback_lambda(
            on_epoch_end = function (epoch, logs) {
              if (epoch %% 20 == 0) cat("\n")
              cat(".")
            })
          ########Fitting the model with Early stopping#######
          early_stop <- callback_early_stopping(monitor = "val_loss",mode = 'min',patience =50)
          ###########Fit of the model for each values of the grid#################
          if(length(trait)==2){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }

          if(length(trait)==3){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2],y_trII[,3]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }

          if(length(trait)==4){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2],y_trII[,3],y_trII[,4]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }

          if(length(trait)==5){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2],y_trII[,3],y_trII[,4],y_trII[,5]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }
          ###########Saving the output of each hyperparameter###################
          No.Epoch_Min = length(model_fit_Final$metrics$val_yhat1_mse)
          Min_MSE = model_fit_Final$metrics$val_yhat1_mse[No.Epoch_Min]
          Tab_pred_MSE[i,stage] = Min_MSE
          Tab_pred_Units[i,stage] = units_M
          Tab_pred_Epoch[i,stage] = No.Epoch_Min[1]
          Tab_pred_Drop[i,stage] = Drop_per
          Tab_pred_Alpha[i,stage] = Alpha
        }
      }
      ##############Selecting the optimal hyperparameters##############
      Median_MSE_Inner = apply(Tab_pred_MSE,2,median)
      Units_Inner = apply(Tab_pred_Units,2,max)
      Drop_Inner = apply(Tab_pred_Drop,2,max)
      Epoch_Inner = apply(Tab_pred_Epoch,2,median)
      Alpha_Inner = apply(Tab_pred_Alpha,2,median)
      Pos_Min_MSE = which(Median_MSE_Inner==min(Median_MSE_Inner))
      Units_O=Units_Inner[Pos_Min_MSE]
      Epoch_O = Epoch_Inner[Pos_Min_MSE]
      Drop_O = Drop_Inner[Pos_Min_MSE]
      Alpha_O = Alpha_Inner[Pos_Min_MSE]
      ###########Refitting the model with the optimal values#################
      input_Sec<- layer_input(shape=dim(X_trn)[2],name="covars")
      #Units_O=33
      #Drop_O=0.05
      #Epoch_O=500
      output_Sec<-input_Sec %>%
        layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l1(l = Alpha_O)) %>%
        layer_dropout(rate = Drop_O)%>%
        layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l1(l = Alpha_O)) %>%
        layer_dropout(rate = Drop_O)%>%
        layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l1(l = Alpha_O)) %>%
        layer_dropout(rate = Drop_O)

      if(length(trait)==2){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }

      if(length(trait)==3){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")
        # add output 3
        yhat3_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat3")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec,yhat3_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2],y_trn[,3]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }

      if(length(trait)==4){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")
        # add output 3
        yhat3_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat3")
        # add output 4
        yhat4_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat4")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec,yhat3_Sec,yhat4_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2],y_trn[,3],y_trn[,4]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }

      if(length(trait)==5){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")
        # add output 3
        yhat3_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat3")
        # add output 4
        yhat4_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat4")
        # add output 4
        yhat5_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat5")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec,yhat3_Sec,yhat4_Sec,yhat5_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2],y_trn[,3],y_trn[,4],y_trn[,5]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }
      ####Prediction of testing set##########################.
      Yhat = model_Sec %>% predict(X_tst)%>%
        data.frame() %>%
        setNames(colnames(y_trn))
      predB=Yhat
      y_p=predB

      #Y_all_tst = data.frame(cbind(y_tst, y_p))

      ########Plots o
      ###########Saving the predicctions of each outer testing set#################
      for(j in 1:length(trait)){
        results=rbind(results,data.frame(Position = tst_set,
                                         Trait = Names_Traits[j],
                                         Partition = o,
                                         Genotype = Y2$Genotype[tst_set],
                                         Environment = Y2$ENV[tst_set],
                                         Units = Units_O,
                                         Epochs = Epoch_O,
                                         Drop_Out = Drop_O,
                                         Alpha_Out = Alpha_O,
                                         Observed = round(y_tst[,j], digits), #$response, digits),
                                         Predicted = round(y_p[,j], digits)
        ))

      }
      cat("CV = ",o,"\n")
    }

    Pred_Summary=summary.BMTMECV(results,information = 'compact',digits = digits)
    Pred_Summary1=summary.BMTMECV(results,information = 'complete',digits = digits)
    Results=list(Results=results,Compact=Pred_Summary,Complete=Pred_Summary1)
    #############Average of prediction performance##################################
    return(Results)
  }

  if(model=="BME"){
    if(UN==TRUE){
      Phenotype=matrix$MTME$YUN
      #Y <- matrix$BMTME$Y
      #Y[-fold_indices,] <- NA
      #Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE

      Genotype=cbind(Z_E,K_expanded,K_GE)
      #Phenotype$Genotype=as.factor(Phenotype$Genotype)
      #Phenotype$Env=as.factor(Phenotype$Env)
      Y2=Phenotype[complete.cases(Phenotype),]
      X2=Genotype[complete.cases(Phenotype),]
    }else{
      Phenotype=complete(matrix$BMTME$YUN, Genotype, ENV)
      X=matrix$BMTME$X
      Z1=matrix$BMTME$Z1
      Z2=matrix$BMTME$Z2
      Genotype=cbind(X,Z1,Z2)
      Y2=Phenotype[complete.cases(Phenotype),]
      X2=Genotype[complete.cases(Phenotype),]
    }
    results = data.frame()
    for(j in 1:length(trait)){
      ############Selecting the response variable#######################

      #Y=as.matrix(Y2[,-c(1, 2)])
      ####Training testing sets using the BMTMEpackage###############
      pheno=data.frame(Line = Y2[, 1],
                       Env=Y2[, 2],
                       Response=Y2[,trait[j]])
      pheno$Line=as.character(pheno$Line)
      pheno$Env=as.character(pheno$Env)
      CrossV=CV.KFold(pheno, DataSetID = 'Line', K=K)
      #fold_list <- make_CV_sets(length(phenotypes), k = folds)
      #######Final X and y for fitting the model###################
      y=Y2[,trait[j]]
      X=X2
      #########Grid of hyperparameters######## ############# #####################.
      if(is.null(Stage)){
        Stage=expand.grid(units_M=c(16,32,64),epochs_M = 1000, Drop_per = c(0.0,0.05,0.15),
                          alpha=c(0.001,0.01,0.4))
      }
      ############Outer Cross-validation#######################
      digits = digits
      Names_Traits = trait[j]

      t = 1

      if(trait[j]=="IT"){
        y = replace(y, y < 0, 0)
        y = replace(y, y > 9, 9)
      }
      if(trait[j]=="SEV"){
        y = replace(y, y < 0, 0)
        y = replace(y, y > 100, 100)
      }

      for (o in 1:K){
        tst_set = CrossV$CrossValidation_list[[o]]
        #tst_set = which(fold_list[[o]])
        length(tst_set)
        X_trn = (X[-tst_set,])
        X_tst = (X[tst_set,])
        y_trn = sqrt((y[-tst_set]))
        y_tst = sqrt((y [tst_set]))
        nCVI = nCVI ####Number of folds for inner CV
        #i = 1
        #########Matrices for saving the output of inner CV#######################.
        Tab_pred_MSE = matrix (NA,ncol = length(Stage[,1]),
                               nrow = nCVI)
        Tab_pred_Epoch = matrix (NA,ncol = length (Stage[,1]),
                                 nrow = nCVI)
        Tab_pred_Units = matrix (NA,ncol = length (Stage[,1]),
                                 nrow = nCVI)
        Tab_pred_Drop = matrix (NA,ncol = length (Stage[,1]),
                                nrow = nCVI)
        Tab_pred_Alpha = matrix (NA,ncol = length (Stage[,1]),
                                 nrow = nCVI)
        X_trI = X_trn
        y_trI = y_trn
        for(i in 1:nCVI){
          for(stage in seq_len(dim(Stage)[1])){
            X_trII = X_trI
            y_trII = y_trI
            units_M <- Stage[stage, 1]
            epochs_M <- Stage[stage, 2]
            Drop_per = Stage[stage, 3]
            Alpha=Stage[stage, 4]
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per) %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per) %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per) %>%
                layer_dense(units = 1)

              model<-keras_model(input,output) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}

            model<-build_model()
            model %>% summary()
            print_dot_callback <- callback_lambda(
              on_epoch_end = function (epoch, logs) {
                if (epoch %% 20 == 0) cat("\n")
                cat(".")
              })
            ########Fitting the model with Early stopping#######
            early_stop <- callback_early_stopping(monitor = "val_loss",mode = 'min',patience =50)
            ###########Fit of the model for each values of the grid#################
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              X_trII, y_trII,
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
            ###########Saving the output of each hyperparameter###################
            No.Epoch_Min = length(model_fit_Final$metrics$val_mse)
            Min_MSE = model_fit_Final$metrics$val_mse[No.Epoch_Min]
            Tab_pred_MSE[i,stage] = Min_MSE
            Tab_pred_Units[i,stage] = units_M
            Tab_pred_Epoch[i,stage] = No.Epoch_Min[1]
            Tab_pred_Drop[i,stage] = Drop_per
            Tab_pred_Alpha[i,stage] = Alpha
          }
        }
        ##############Selecting the optimal hyperparameters##############
        Median_MSE_Inner = apply(Tab_pred_MSE,2,median)
        Units_Inner = apply(Tab_pred_Units,2,max)
        Drop_Inner = apply(Tab_pred_Drop,2,max)
        Epoch_Inner = apply(Tab_pred_Epoch,2,median)
        Alpha_Inner = apply(Tab_pred_Alpha,2,median)
        Pos_Min_MSE = which(Median_MSE_Inner==min(Median_MSE_Inner))
        Units_O=Units_Inner[Pos_Min_MSE]
        Epoch_O = Epoch_Inner[Pos_Min_MSE]
        Drop_O = Drop_Inner[Pos_Min_MSE]
        Alpha_O = Alpha_Inner[Pos_Min_MSE]
        ###########Refitting the model with the optimal values#################
        input_Sec<- layer_input(shape=dim(X_trn)[2],name="covars")


        output_Sec<-input_Sec %>%
          layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha_O)) %>%
          layer_dropout(rate = Drop_O) %>%
          layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha_O)) %>%
          layer_dropout(rate = Drop_O) %>%
          layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha_O)) %>%
          layer_dropout(rate = Drop_O) %>%
          layer_dense(units =1)

        model_Sec=keras_model(input_Sec,output_Sec) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          X_trn, y_trn,
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
        ####Prediction of testing set##########################.
        Yhat = model_Sec %>% predict(X_tst)
        y_p = Yhat
        y_p_tst = as.numeric(y_p)
        ###########Saving the predicctions of each outer testing set#################
        results=rbind(results,data.frame(Position = tst_set,
                                         Trait = trait[j],
                                         Partition = o,
                                         Genotype = Y2$Genotype[tst_set],
                                         Environment = Y2$ENV[tst_set],
                                         Units = Units_O,
                                         Epochs = Epoch_O,
                                         Drop_Out = Drop_O,
                                         Alpha_Out = Alpha_O,
                                         Observed = round(y_tst, digits), #$response, digits),
                                         Predicted = round(y_p_tst, digits)
        ))
        cat("CV = ",o,"\n")
      }
    }
    Pred_Summary=summary.BMTMECV(results,information = 'compact',digits = digits)
    Pred_Summary1=summary.BMTMECV(results,information = 'complete',digits = digits)
    Results=list(Results=results,Compact=Pred_Summary,Complete=Pred_Summary1)
    #############Average of prediction performance##################################
    return(Results)
  }

  if(model=="BMTME"){
    if(UN==TRUE){
      Phenotype=matrix$MTME$YUN
      #Y <- matrix$BMTME$Y
      #Y[-fold_indices,] <- NA
      #Z_L=matrix$MTME$Z_L
      Z_E=matrix$MTME$Z_E
      K_expanded=matrix$MTME$K_expanded
      K_GE=matrix$MTME$K_GE

      Genotype=cbind(Z_E,K_expanded,K_GE)
      #Phenotype$Genotype=as.factor(Phenotype$Genotype)
      #Phenotype$Env=as.factor(Phenotype$Env)
      Y2=Phenotype[complete.cases(Phenotype),]
      X2=Genotype[complete.cases(Phenotype),]
    }else{
      Phenotype=complete(matrix$BMTME$YUN, Genotype, ENV)
      X=matrix$BMTME$X
      Z1=matrix$BMTME$Z1
      Z2=matrix$BMTME$Z2
      Genotype=cbind(X,Z1,Z2)
      Y2=Phenotype[complete.cases(Phenotype),]
      X2=Genotype[complete.cases(Phenotype),]
    }
    #Y2=Y2[1:50,]
    #X2=X2[1:50,]
    ############Selecting the response variable#######################

    #Y=as.matrix(Y2[,-c(1, 2)])
    ####Training testing sets using the BMTMEpackage###############
    pheno=data.frame(Line = Y2[, 1],
                     Env=Y2[, 2],
                     Response=Y2[,trait[1]])
    pheno$Line=as.character(pheno$Line)
    pheno$Env=as.character(pheno$Env)
    CrossV=CV.KFold(pheno, DataSetID = 'Line', K=K)
    #fold_list <- make_CV_sets(length(phenotypes), k = folds)
    #######Final X and y for fitting the model###################
    y=Y2[,trait]
    X=X2
    dim(X)
    dim(y)
    #########Grid of hyperparameters######## ############# #####################.
    if(is.null(Stage)){
      Stage=expand.grid(units_M=c(16,32,64),epochs_M = 1000, Drop_per = c(0.0,0.05,0.15),
                        alpha=c(0.001,0.01,0.4))
    }
    ############Outer Cross-validation#######################
    digits = digits
    Names_Traits = trait
    results = data.frame()
    t = 1

    y$IT=replace(y$IT, y$IT < 0, 0)
    y$IT=replace(y$IT, y$IT > 9, 9)
    y$SEV=replace(y$SEV, y$SEV < 0, 0)
    y$SEV=replace(y$SEV, y$SEV > 10, 10)


    for (o in 1:K){
      tst_set = CrossV$CrossValidation_list[[o]]
      #tst_set = which(fold_list[[o]])
      length(tst_set)
      X_trn = (X[-tst_set,])
      X_tst = (X[tst_set,])
      y_trn = sqrt((y[-tst_set,]))
      y_tst = sqrt((y[tst_set,]))
      nCVI = nCVI ####Number of folds for inner CV
      #i = 1
      #########Matrices for saving the output of inner CV#######################.
      Tab_pred_MSE = matrix (NA,ncol = length(Stage[,1]),
                             nrow = nCVI)
      Tab_pred_Epoch = matrix (NA,ncol = length (Stage[,1]),
                               nrow = nCVI)
      Tab_pred_Units = matrix (NA,ncol = length (Stage[,1]),
                               nrow = nCVI)
      Tab_pred_Drop = matrix (NA,ncol = length (Stage[,1]),
                              nrow = nCVI)
      Tab_pred_Alpha = matrix (NA,ncol = length (Stage[,1]),
                               nrow = nCVI)
      X_trI = X_trn
      y_trI = y_trn
      for(i in 1:nCVI){
        for(stage in seq_len(dim(Stage)[1])){
          X_trII = X_trI
          y_trII = y_trI
          units_M <- Stage[stage, 1]
          epochs_M <- Stage[stage, 2]
          Drop_per = Stage[stage, 3]
          Alpha = Stage[stage, 4]

          if(length(trait)==2){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")

              model<-keras_model(input,list(yhat1,yhat2)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }

          if(length(trait)==3){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")
              # add output 2
              yhat3 <- output %>%
                layer_dense(units = 1, name="yhat3")

              model<-keras_model(input,list(yhat1,yhat2,yhat3)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }

          if(length(trait)==4){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")
              # add output 2
              yhat3 <- output %>%
                layer_dense(units = 1, name="yhat3")
              # add output 2
              yhat4 <- output %>%
                layer_dense(units = 1, name="yhat4")

              model<-keras_model(input,list(yhat1,yhat2,yhat3,yhat4)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }

          if(length(trait)==5){
            build_model<-function() {
              input<-layer_input(shape=dim(X_trII)[2],name="covars")

              output <- input %>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)%>%
                layer_dense(units = units_M, activation = "relu",kernel_regularizer =regularizer_l2(l = Alpha)) %>%
                layer_dropout(rate = Drop_per)

              # add output 1
              yhat1 <- output %>%
                layer_dense(units = 1, name="yhat1")
              # add output 2
              yhat2 <- output %>%
                layer_dense(units = 1, name="yhat2")
              # add output 2
              yhat3 <- output %>%
                layer_dense(units = 1, name="yhat3")
              # add output 2
              yhat4 <- output %>%
                layer_dense(units = 1, name="yhat4")
              # add output 2
              yhat5 <- output %>%
                layer_dense(units = 1, name="yhat5")

              model<-keras_model(input,list(yhat1,yhat2,yhat3,yhat4,yhat5)) %>% compile(
                loss = "mse",
                #optimizer = optimizer_adam(),
                optimizer = "rmsprop",
                metrics = c("mse","accuracy"))
              model}
          }


          model<-build_model()
          model %>% summary()
          print_dot_callback <- callback_lambda(
            on_epoch_end = function (epoch, logs) {
              if (epoch %% 20 == 0) cat("\n")
              cat(".")
            })
          ########Fitting the model with Early stopping#######
          early_stop <- callback_early_stopping(monitor = "val_loss",mode = 'min',patience =50)
          ###########Fit of the model for each values of the grid#################
          if(length(trait)==2){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }

          if(length(trait)==3){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2],y_trII[,3]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }

          if(length(trait)==4){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2],y_trII[,3],y_trII[,4]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }

          if(length(trait)==5){
            model_Final<-build_model()
            model_fit_Final<-model_Final %>% fit(
              x=X_trII, y=list(y_trII[,1],y_trII[,2],y_trII[,3],y_trII[,4],y_trII[,5]),
              epochs = epochs_M,
              batch_size =72,
              ###shuffled = F,
              validation_split = 0.2,
              verbose = 0,
              #verbose = 1,
              #view_metrics = TRUE,
              callbacks = list(early_stop, print_dot_callback))
          }
          ###########Saving the output of each hyperparameter###################
          No.Epoch_Min = length(model_fit_Final$metrics$val_yhat1_mse)
          Min_MSE = model_fit_Final$metrics$val_yhat1_mse[No.Epoch_Min]
          Tab_pred_MSE[i,stage] = Min_MSE
          Tab_pred_Units[i,stage] = units_M
          Tab_pred_Epoch[i,stage] = No.Epoch_Min[1]
          Tab_pred_Drop[i,stage] = Drop_per
          Tab_pred_Alpha[i,stage] = Alpha
        }
      }
      ##############Selecting the optimal hyperparameters##############
      Median_MSE_Inner = apply(Tab_pred_MSE,2,median)
      Units_Inner = apply(Tab_pred_Units,2,max)
      Drop_Inner = apply(Tab_pred_Drop,2,max)
      Epoch_Inner = apply(Tab_pred_Epoch,2,median)
      Alpha_Inner = apply(Tab_pred_Alpha,2,median)
      Pos_Min_MSE = which(Median_MSE_Inner==min(Median_MSE_Inner))
      Units_O=Units_Inner[Pos_Min_MSE]
      Epoch_O = Epoch_Inner[Pos_Min_MSE]
      Drop_O = Drop_Inner[Pos_Min_MSE]
      Alpha_O = Alpha_Inner[Pos_Min_MSE]
      ###########Refitting the model with the optimal values#################
      input_Sec<- layer_input(shape=dim(X_trn)[2],name="covars")
      #Units_O=33
      #Drop_O=0.05
      #Epoch_O=500
      output_Sec<-input_Sec %>%
        layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l1(l = Alpha_O)) %>%
        layer_dropout(rate = Drop_O)%>%
        layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l1(l = Alpha_O)) %>%
        layer_dropout(rate = Drop_O)%>%
        layer_dense(units = Units_O, activation = "relu",kernel_regularizer =regularizer_l1(l = Alpha_O)) %>%
        layer_dropout(rate = Drop_O)

      if(length(trait)==2){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec)) %>% compile(
          loss = "mean_squared_error",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mean_squared_error","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }

      if(length(trait)==3){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")
        # add output 3
        yhat3_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat3")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec,yhat3_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2],y_trn[,3]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }

      if(length(trait)==4){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")
        # add output 3
        yhat3_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat3")
        # add output 4
        yhat4_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat4")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec,yhat3_Sec,yhat4_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2],y_trn[,3],y_trn[,4]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }

      if(length(trait)==5){
        # add output 1
        yhat1_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat1")
        # add output 2
        yhat2_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat2")
        # add output 3
        yhat3_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat3")
        # add output 4
        yhat4_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat4")
        # add output 4
        yhat5_Sec <- output_Sec %>%
          layer_dense(units = 1, name="yhat5")

        model_Sec=keras_model(input_Sec,list(yhat1_Sec,yhat2_Sec,yhat3_Sec,yhat4_Sec,yhat5_Sec)) %>% compile(
          loss = "mse",
          #optimizer = optimizer_adam(),
          optimizer = "rmsprop",
          metrics = c("mse","accuracy"))


        ModelFited <-model_Sec %>% fit(
          x=X_trn, y=list(y_trn[,1],y_trn[,2],y_trn[,3],y_trn[,4],y_trn[,5]),
          epochs = Epoch_O, batch_size =30,
          ####### validation_split = 0.2,early_stop,
          verbose = 0,
          #verbose = 1,
          #view_metrics = TRUE,
          callbacks=list(print_dot_callback))
      }
      ####Prediction of testing set##########################.
      Yhat = model_Sec %>% predict(X_tst)%>%
        data.frame() %>%
        setNames(colnames(y_trn))
      predB=Yhat
      y_p=predB

      #Y_all_tst = data.frame(cbind(y_tst, y_p))

      ########Plots o
      ###########Saving the predicctions of each outer testing set#################
      for(j in 1:length(trait)){
        results=rbind(results,data.frame(Position = tst_set,
                                         Trait = Names_Traits[j],
                                         Partition = o,
                                         Genotype = Y2$Genotype[tst_set],
                                         Environment = Y2$ENV[tst_set],
                                         Units = Units_O,
                                         Epochs = Epoch_O,
                                         Drop_Out = Drop_O,
                                         Alpha_Out = Alpha_O,
                                         Observed = round(y_tst[,j], digits), #$response, digits),
                                         Predicted = round(y_p[,j], digits)
        ))

      }
      cat("CV = ",o,"\n")
    }

    Pred_Summary=summary.BMTMECV(results,information = 'compact',digits = digits)
    Pred_Summary1=summary.BMTMECV(results,information = 'complete',digits = digits)
    Results=list(Results=results,Compact=Pred_Summary,Complete=Pred_Summary1)
    #############Average of prediction performance##################################
    return(Results)
  }
}

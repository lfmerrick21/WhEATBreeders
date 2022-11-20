WHEAT<-function(Phenotype,
                Genotype,
                QC=TRUE,
                GS=TRUE,
                #QC Info
                Geno_Type="Hapmap",
                Imputation="Beagle",
                Filter=TRUE,
                Missing_Rate=0.20,
                MAF=0.05,
                Filter_Ind=TRUE,
                Missing_Rate_Ind=0.80,
                #If QC Info is FALSE
                GDre=NULL,
                GT=NULL,
                GIre=NULL,
                GBS_Train=NULL,
                GBS_Predict=NULL,
                Matrix=NULL,
                #GS Info
                Type="Regression",
                Replications=1,
                Training="F5_2015",
                Prediction="DH_2020",
                CV=NULL,
                PC=NULL,
                Trait="EM",
                Study="Tutorial",
                Outcome="Tested", #Tested or Untested
                Trial=c("F5_2015","DH_2020","BL_2015_2020"),
                Scheme="K-Fold",
                Method="Two-Step",
                Package="rrBLUP",
                model="rrBLUP",
                Kernel="Markers",
                markers=NULL,
                folds = 5,
                nIter = 1500,
                burnIn = 500,
                Sparse=FALSE,
                m=NULL,
                degree=NULL,
                nL=NULL,
                transformation="none",
                #GAGS Info
                GAGS=FALSE,
                PCA.total=3,
                QTN=10,
                GWAS=c("BLINK"),
                alpha=0.05,
                threshold=NULL,
                GE=TRUE,
                UN=FALSE,
                GE_model="MTME",
                sampling="up",
                repeats=5,
                method="repeatedcv",
                digits=4,
                nCVI=5,
                Messages=TRUE){
                colnames(Phenotype)[1:2]<-c("Genotype","Env")
                Phenotype=Phenotype %>% filter(Env %in% c(Trial))
                  if(QC==TRUE){
                    #############PandG#########################
                    ###############################################################################
                    #Get Taxa
                    #Phenotype=Phenotype1
                    Pheno_Names=Phenotype
                    names(Pheno_Names)[1]<-"Taxa"
                    Pheno_Names$Taxa=as.character(Pheno_Names$Taxa)
                    gname=as.vector(Pheno_Names$Taxa)
                    gname<-clean_names(gname)
                    #Filter Hapmap for Taxa
                    if(Geno_Type=="VCF"){
                      #Genotype1=Genotype
                      Genotype1=data.frame(Genotype@fix,Genotype@gt)
                      genotypes_num<- VCF_2_numeric(Genotype1)
                      GD=genotypes_num$genotypes
                      rownames(GD)<-clean_names_op(rownames(GD))
                      names.use <- rownames(GD)[(rownames(GD) %in% gname)]
                      GD <-GD[row.names(GD) %in% c(names.use),]
                      #GD <- GD[, names.use, with = FALSE]
                      GI=genotypes_num$marker_map
                      #setnames(GI, old = c('rs','chr','pos'),
                               #new = c('SNP','Chromosome','Position'))
                      GI <- GI %>%
                        rename(SNP=rs,
                               Chromosome=chr,
                               Position=pos)
                      GT=data.frame(V1=rownames(GD))
                      #Save Raw Files
                      fwrite(GT,paste0(Study,"_GT_taxa.csv"))
                      fwrite(GD,paste0(Study,"_GD_numeric_Raw.csv"))
                      fwrite(GI,paste0(Study,"_GI_map_Raw.csv"))
                      save(GD,GI,GT,file=paste0(Study,"_Raw.RData"))
                    }

                    if(Geno_Type=="Hapmap"){
                      names.use <- names(Genotype)[(names(Genotype) %in% gname)]
                      hapmap <- Genotype[, names.use, with = FALSE]
                      hapmap=cbind(Genotype[,1:11],hapmap)
                      #str(output)
                      hapmap[hapmap=="NA"]<-"NA"
                      hapmap$`assembly#`=as.character(hapmap$`assembly#`)
                      hapmap$center=as.character(hapmap$center)
                      hapmap$protLSID=as.character(hapmap$protLSID)
                      hapmap$assayLSID=as.character(hapmap$assayLSID)
                      hapmap$panelLSID=as.character(hapmap$panelLSID)
                      hapmap$QCcode=as.character(hapmap$QCcode)
                      dim(hapmap)
                      #Save filtered hapmap
                      save(hapmap,file=paste0(Study,"_Filtered_Hapmap.RData"))
                      fwrite(hapmap,paste0(Study,"_Filtered_Hapmap.hmp.txt"),sep="\t",row.names = FALSE,col.names = TRUE)
                      #Convert to Numeric
                      if(Imputation=="Middle"){
                        outG=GAPIT.HapMap(hapmap,SNP.impute="Middle")
                      }else{
                        outG=GAPIT.HapMap(hapmap,SNP.impute="None")
                      }
                      GT=outG$GT
                      GT=data.frame(V1=rownames(GT),GT)
                      GD=outG$GD
                      GI=outG$GI
                      rownames(GD)=rownames(GT)
                      colnames(GD)=GI$SNP
                      #Save Raw Files
                      fwrite(GT,paste0(Study,"_GT_taxa.csv"))
                      fwrite(GD,paste0(Study,"_GD_numeric_Raw.csv"))
                      fwrite(GI,paste0(Study,"_GI_map_Raw.csv"))
                      save(GD,GI,GT,file=paste0(Study,"_Raw.RData"))
                    }

                    #### Filter
                    ##### Remove makers based on missing data
                    if(Filter==FALSE){
                      GDmf = GD
                      GImf = GI
                    }else{
                      mr=calc_missrate(GD)
                      mr_indices <- which(mr > Missing_Rate)
                      if(length(mr_indices)!=0){
                        GDmr = GD[,-mr_indices]
                        GImr = GI[-mr_indices,]
                      }else{
                        GDmr = GD
                        GImr = GI
                      }
                      dim(GDmr)
                      dim(GImr)

                      ##### Remove Monomorphic Markers
                      maf <- calc_maf_apply(GDmr, encoding = c(0, 1, 2))
                      mono_indices <- which(maf ==0)
                      if(length(mono_indices)!=0){
                        GDmo = GDmr[,-mono_indices]
                        GImo = GImr[-mono_indices,]
                      }else{
                        GDmo = GDmr
                        GImo = GImr
                      }
                      dim(GDmo)
                      dim(GImo)

                      ##### Remove MAF
                      maf <- calc_maf_apply(GDmo, encoding = c(0, 1, 2))
                      mono_indices <- which(maf < MAF)
                      if(length(mono_indices)!=0){
                        GDmf = GDmo[,-mono_indices]
                        GImf = GImo[-mono_indices,]
                      }else{
                        GDmf = GDmo
                        GImf = GImo
                      }
                      dim(GDmf)
                      dim(GImf)
                    }

                    if(Filter_Ind==FALSE){
                      GDmf = GDmf
                    }else{
                      mr=calc_missrate(t(GDmf))
                      mr_indices <- which(mr > Missing_Rate_Ind)
                      if(length(mr_indices)!=0){
                        GDmf = GDmf[-mr_indices,]
                        GT=GT[-mr_indices,]
                      }else{
                        GDmf = GDmf
                        GT=GT
                      }
                      dim(GDmf)
                      dim(GImf)
                    }
                    #### Imputation
                    ##### Imputation Using BEAGLE
                    if(Imputation=="None"){
                      GDEre = GDmf
                      GIEre = GImf
                      dim(GDEre)
                      dim(GIEre)

                      fwrite(GDre,paste0(Study,"_GD_Filt_Imputed.csv"))
                      fwrite(GIre[,1:3],paste0(Study,"_GI_Filt_Imputed.csv"))
                      save(GDre,GIre,GT,file=paste0(Study,"_Filt_Imputed.RData"))
                    }else{

                      if(Imputation=="Beagle"){
                        GImf$chr <- GImf$Chromosome
                        GImf$rs <- GImf$SNP
                        GImf$pos <- GImf$Position

                        LD_file <- paste0(Study)
                        vcf_LD <- numeric_2_VCF(GDmf, GImf)
                        write_vcf(vcf_LD, outfile = LD_file)#Exports vcf
                        # Assign parameters
                        genotype_file = paste0(Study,".vcf")
                        outfile = paste0(Study,"_imp")

                        # Define a system command
                        command1_prefix <- "java -jar beagle.25Nov19.28d.jar"
                        command_args <- paste(" gt=", genotype_file, " out=", outfile, sep = "")
                        command1 <- paste(command1_prefix, command_args)
                        test <- system(command1)
                        # Run BEAGLE using the system function, this will produce a gzip .vcf file
                        Xvcf=read.vcfR(paste0(Study,"_imp",".vcf.gz"))
                        Xbgl=data.frame(Xvcf@fix,Xvcf@gt)
                        genotypes_imp <- VCF_2_numeric(Xbgl)[[1]]
                        #Pre-imputed genotype matrix
                        dim(GDmf)
                        sum(is.na(GDmf))
                        #Imputed genotype matrix
                        dim(genotypes_imp)
                        sum(is.na(genotypes_imp))

                        #Remove MAF <5%
                        maf <- calc_maf_apply(genotypes_imp, encoding = c(0, 1, 2))
                        mono_indices <- which(maf < MAF)

                        if(length(mono_indices)!=0){
                          GDre = genotypes_imp[,-mono_indices]
                          GIre = GImf[-mono_indices,]
                        }else{
                          GDre = genotypes_imp
                          GIre = GImf
                        }

                        GIre=GIre[,1:3]
                        dim(GDre)
                        dim(GIre)

                        fwrite(GDre,paste0(Study,"_GD_Filt_Imputed.csv"))
                        fwrite(GIre[,1:3],paste0(Study,"_GI_Filt_Imputed.csv"))
                        save(GDre,GIre,GT,file=paste0(Study,"_Filt_Imputed.RData"))
                        dim(GDre)
                        dim(GIre)
                        dim(GT)

                      }

                      if(Imputation=="KNN"){
                        #Section for downloading package
                        #if (!requireNamespace("BiocManager", quietly = TRUE))
                          #install.packages("BiocManager")

                        #BiocManager::install("impute")

                        x=impute::impute.knn(as.matrix(t(GDmf)))
                        myGD_imp=t(x$data)
                        myGD_imp<-round(myGD_imp,0)
                        sum(is.na(myGD_imp))
                        #Filter again
                        #Remove MAF <5%
                        maf <- calc_maf_apply(myGD_imp, encoding = c(0, 1, 2))
                        mono_indices <- which(maf < MAF)
                        if(length(mono_indices)!=0){
                          GDre = myGD_imp[,-mono_indices]
                          GIre = GImf[-mono_indices,]
                        }else{
                          GDre = myGD_imp
                          GIre = GImf
                        }
                        dim(GDre)
                        dim(GIre)
                        GIre=GIre[,1:3]

                        fwrite(GDre,paste0(Study,"_GD_Filt_Imputed.csv"))
                        fwrite(GIre[,1:3],paste0(Study,"_GI_Filt_Imputed.csv"))
                        save(GDre,GIre,GT,file=paste0(Study,"_Filt_Imputed.RData"))
                      }
                    }
                    ##########PG################################
                    if(Method=="Two-Step"){
                      if(Outcome=="Tested"){
                      if(!is.null(CV)){
                        GBS.list=c()
                        for(i in 1:length(Trial)){
                          Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                          Pheno=data.frame(Pheno)
                          Namescol=c("Genotype","Env",Trait)
                          Pheno=Pheno[,Namescol]
                          Namescole=c("Genotype",Trait)
                          #for(j in 1:length(Trait)){
                            #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),]
                            #GBS<<-PGandCV(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre,CV)
                            #mv(from = "GBS", to = paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),envir = globalenv())
                            #GBS<<-PGandCV(Pheno[,c("Genotype",Trait)],GDre,GT,GIre,CV)
                            GBS<<-PGandCV(Pheno[,Namescole],GDre,GT,GIre,CV)
                            mv(from = "GBS", to = paste0("GBS_2_CV_",Trial[i]),envir = globalenv())

                            #GBS.list=c(GBS.list,paste0("GBS_2_CV_",Trial[i],"_",Trait[j]))
                            GBS.list=c(GBS.list,paste0("GBS_2_CV_",Trial[i]))
                            #save(list=paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_CV_",Trial[i],"_",Trait[j],".RData"))
                          #}

                        }
                        save(list = GBS.list, file=paste0("GBS_2_CV_",Study,".RData"))
                      }else{
                        GBS.list=c()
                        for(i in 1:length(Trial)){
                          Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                          Pheno=data.frame(Pheno)
                          Namescol=c("Genotype","Env",Trait)
                          Pheno=Pheno[,Namescol]
                          Namescole=c("Genotype",Trait)
                          #for(j in 1:length(Trait)){
                            #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),]
                            #GBS<<-PandG(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre)
                            #mv(from = "GBS", to = paste0("GBS_2_",Trial[i],"_",Trait[j]),envir = globalenv())
                            #GBS<<-PandG(Pheno[,c("Genotype",Trait)],GDre,GT,GIre)
                            GBS<<-PandG(Pheno[,Namescole],GDre,GT,GIre)
                            mv(from = "GBS", to = paste0("GBS_2_",Trial[i]),envir = globalenv())

                            #GBS.list=c(GBS.list,paste0("GBS_2_",Trial[i],"_",Trait[j]))
                            GBS.list=c(GBS.list,paste0("GBS_2_",Trial[i]))
                            #save(list=paste0("GBS_2_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_",Trial[i],"_",Trait[j],".RData"))
                          #}

                        }
                        save(list = GBS.list, file=paste0("GBS_2_",Study,".RData"))
                      }
                      }

                      if(Outcome=="Untested"){
                        if(!is.null(CV)){
                          GBS.list=c()
                          for(i in 1:length(Trial)){
                            Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                            Pheno=data.frame(Pheno)
                            Namescol=c("Genotype","Env",Trait)
                            Pheno=Pheno[,Namescol]
                            Namescole=c("Genotype",Trait)
                            #for(j in 1:length(Trait)){
                              #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
                              #GBS<<-PGandCV(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre,CV)
                              #mv(from = "GBS", to = paste0("GBS_2_CV_Untested_",Trial[i],"_",Trait[j]),envir = globalenv())
                              #GBS<<-PGandCV(Pheno[,c("Genotype",Trait)],GDre,GT,GIre,CV)
                              GBS<<-PGandCV(Pheno[,Namescole],GDre,GT,GIre,CV)
                              mv(from = "GBS", to = paste0("GBS_2_CV_Untested_",Trial[i]),envir = globalenv())

                              #GBS.list=c(GBS.list,paste0("GBS_2_CV_Untested_",Trial[i],"_",Trait[j]))
                              GBS.list=c(GBS.list,paste0("GBS_2_CV_Untested_",Trial[i]))
                              #save(list=paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_CV_",Trial[i],"_",Trait[j],".RData"))
                            #}

                          }
                          save(list = GBS.list, file=paste0("GBS_2_CV_Untested_",Study,".RData"))
                        }else{
                          GBS.list=c()
                          for(i in 1:length(Trial)){
                            Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                            Pheno=data.frame(Pheno)
                            Namescol=c("Genotype","Env",Trait)
                            Pheno=Pheno[,Namescol]
                            Namescole=c("Genotype",Trait)
                            #for(j in 1:length(Trait)){
                              #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
                              #GBS<<-PandG(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre)
                              #mv(from = "GBS", to = paste0("GBS_2_Untested_",Trial[i],"_",Trait[j]),envir = globalenv())
                              #GBS<<-PandG(Pheno[,c("Genotype",Trait)],GDre,GT,GIre)
                              GBS<<-PandG(Pheno[,Namescole],GDre,GT,GIre)
                              mv(from = "GBS", to = paste0("GBS_2_Untested_",Trial[i]),envir = globalenv())

                              #GBS.list=c(GBS.list,paste0("GBS_2_Untested_",Trial[i],"_",Trait[j]))
                              GBS.list=c(GBS.list,paste0("GBS_2_Untested_",Trial[i]))
                              #save(list=paste0("GBS_2_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_",Trial[i],"_",Trait[j],".RData"))
                            #}

                          }
                          save(list = GBS.list, file=paste0("GBS_2_Untested_",Study,".RData"))
                        }
                      }
                    }

                    if(Method=="One-Step"){
                      if(!is.null(CV)){
                        Pheno=Phenotype %>% filter(Env %in% c(Trial))
                        Pheno=Phenotype[,c("Genotype","Env",Trait)]
                        GBS<<-PGEandCV(Pheno,GDre,GT,GIre,CV)
                        mv(from = "GBS", to = paste0("GBS_1_CV_",Study),envir = globalenv())
                        save(list=paste0("GBS_1_CV_",Study),file=paste0("GBS_1_CV_",Study,".RData"))
                      }else{
                        Pheno=Phenotype %>% filter(Env %in% c(Trial))
                        Pheno=Phenotype[,c("Genotype","Env",Trait)]
                        GBS<<-PGandE(Pheno,GDre,GT,GIre)
                        mv(from = "GBS", to = paste0("GBS_1_",Study),envir = globalenv())
                        save(list=paste0("GBS_1_",Study),file=paste0("GBS_1_",Study,".RData"))
                      }
                    }
                    if(Method=="One-Step"){
                      if(Scheme=="K-Fold"){
                        Matrix<<-GE_Matrix_IND(genotypes=get(paste0("GBS_1_",Study))$geno, phenotype=get(paste0("GBS_1_",Study))$pheno,trait=Trait,Kernel=Kernel,GE=GE,model=GE_model,Sparse=Sparse,m=m,degree=degree, nL=nL)
                        mv(from = "Matrix", to = paste0("Matrix_",Study,GE_model),envir = globalenv())
                        save(list=paste0("Matrix_",Study),file=paste0("Matrix_",Study,".RData"))
                      }else{
                        Matrix<<-GE_Matrix_IND(train_genotypes=get(paste0("GBS_1_",Study))$geno, train_phenotype=get(paste0("GBS_1_",Study))$pheno,test_genotypes=get(paste0("GBS_1_",Study))$geno, test_phenotype=get(paste0("GBS_1_",Study))$pheno,trait=Trait,Kernel=Kernel,GE=GE,model=GE_model,Sparse=Sparse,m=m,degree=degree, nL=nL)
                        mv(from = "Matrix", to = paste0("Matrix_",Study,GE_model),envir = globalenv())
                        save(list=paste0("Matrix_",Study),file=paste0("Matrix_",Study,".RData"))
                      }

                    }


                    #############GS#########################
                    #Load files
                    #if(Method=="Two-Step"){
                      #if(!is.null(CV)){
                        #load(file=paste0("GBS_2_CV_",Study,".RData"))
                      #}else{
                        #load(file=paste0("GBS_2_",Study,".RData"))
                      #}
                    #}

                    #if(Method=="One-Step"){
                      #if(!is.null(CV)){
                        #load(file=paste0("GBS_1_CV_",Study,".RData"))
                      #}else{
                        #load(file=paste0("GBS_1_",Study,".RData"))
                      #}
                    #}


                  }
                  if(GS==TRUE){
                    #####################GS######################################################
                    if(QC==FALSE){
                      if(is.null(GBS_Train)){

                      if(Method=="Two-Step"){
                        if(Outcome=="Tested"){
                        if(!is.null(CV)){
                          GBS.list=c()
                          for(i in 1:length(Trial)){
                            Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                            Pheno=data.frame(Pheno)
                            Namescol=c("Genotype","Env",Trait)
                            Pheno=Pheno[,Namescol]
                            Namescole=c("Genotype",Trait)
                            #for(j in 1:length(Trait)){
                              #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),]
                              #GBS<<-PGandCV(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre,CV)
                              #mv(from = "GBS", to = paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),envir = globalenv())
                              #GBS<<-PGandCV(Pheno[,c("Genotype",Trait)],GDre,GT,GIre,CV)
                              GBS<<-PGandCV(Pheno[,Namescole],GDre,GT,GIre,CV)
                              mv(from = "GBS", to = paste0("GBS_2_CV_",Trial[i]),envir = globalenv())

                              #GBS.list=c(GBS.list,paste0("GBS_2_CV_",Trial[i],"_",Trait[j]))
                              GBS.list=c(GBS.list,paste0("GBS_2_CV_",Trial[i]))
                              #save(list=paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_CV_",Trial[i],"_",Trait[j],".RData"))
                            #}

                          }
                          save(list = GBS.list, file=paste0("GBS_2_CV_",Study,".RData"))
                        }else{
                          GBS.list=c()
                          for(i in 1:length(Trial)){
                            Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                            Pheno=data.frame(Pheno)
                            Namescol=c("Genotype","Env",Trait)
                            Pheno=Pheno[,Namescol]
                            Namescole=c("Genotype",Trait)
                            #for(j in 1:length(Trait)){
                              #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),]
                              #GBS<<-PandG(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre)
                              #mv(from = "GBS", to = paste0("GBS_2_",Trial[i],"_",Trait[j]),envir = globalenv())
                              #GBS<<-PandG(Pheno[,c("Genotype",Trait)],GDre,GT,GIre)
                              GBS<<-PandG(Pheno[,Namescole],GDre,GT,GIre)
                              mv(from = "GBS", to = paste0("GBS_2_",Trial[i]),envir = globalenv())

                              #GBS.list=c(GBS.list,paste0("GBS_2_",Trial[i],"_",Trait[j]))
                              GBS.list=c(GBS.list,paste0("GBS_2_",Trial[i]))
                              #save(list=paste0("GBS_2_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_",Trial[i],"_",Trait[j],".RData"))
                            #}

                          }
                          save(list = GBS.list, file=paste0("GBS_2_",Study,".RData"))
                        }}
                        if(Outcome=="Untested"){
                          if(!is.null(CV)){
                            GBS.list=c()
                            for(i in 1:length(Trial)){
                              Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                              Pheno=data.frame(Pheno)
                              Namescol=c("Genotype","Env",Trait)
                              Pheno=Pheno[,Namescol]
                              Namescole=c("Genotype",Trait)
                              #for(j in 1:length(Trait)){
                                #Pheno<-Pheno[complete.cases(Pheno[,Trait[j]]),c(1,Trait[j])]
                                #GBS<<-PGandCV(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre,CV)
                                #mv(from = "GBS", to = paste0("GBS_2_CV_Untested_",Trial[i],"_",Trait[j]),envir = globalenv())
                                #GBS<<-PGandCV(Pheno[,c("Genotype",Trait)],GDre,GT,GIre,CV)
                                GBS<<-PGandCV(Pheno[,Namescole],GDre,GT,GIre,CV)
                                mv(from = "GBS", to = paste0("GBS_2_CV_Untested_",Trial[i]),envir = globalenv())

                                #GBS.list=c(GBS.list,paste0("GBS_2_CV_Untested_",Trial[i],"_",Trait[j]))
                                GBS.list=c(GBS.list,paste0("GBS_2_CV_Untested_",Trial[i]))
                                #save(list=paste0("GBS_2_CV_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_CV_",Trial[i],"_",Trait[j],".RData"))
                              #}

                            }
                            save(list = GBS.list, file=paste0("GBS_2_CV_Untested_",Study,".RData"))
                          }else{
                            GBS.list=c()
                            for(i in 1:length(Trial)){
                              Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                              Pheno=data.frame(Pheno)
                              Namescol=c("Genotype","Env",Trait)
                              Pheno=Pheno[,Namescol]
                              Namescole=c("Genotype",Trait)
                              #for(j in 1:length(Trait)){
                                #Pheno=Phenotype %>% filter(Env %in% c(Trial[i]))
                                #GBS<<-PandG(Pheno[,c("Genotype",Trait[j])],GDre,GT,GIre)
                                #mv(from = "GBS", to = paste0("GBS_2_Untested_",Trial[i],"_",Trait[j]),envir = globalenv())
                                #GBS<<-PandG(Pheno[,c("Genotype",Trait)],GDre,GT,GIre)
                                GBS<<-PandG(Pheno[,Namescole],GDre,GT,GIre)
                                mv(from = "GBS", to = paste0("GBS_2_Untested_",Trial[i]),envir = globalenv())

                                #GBS.list=c(GBS.list,paste0("GBS_2_Untested_",Trial[i],"_",Trait[j]))
                                GBS.list=c(GBS.list,paste0("GBS_2_Untested_",Trial[i]))
                                #save(list=paste0("GBS_2_",Trial[i],"_",Trait[j]),file=paste0("GBS_2_",Trial[i],"_",Trait[j],".RData"))
                              #}

                            }
                            save(list = GBS.list, file=paste0("GBS_2_Untested_",Study,".RData"))
                          }
                        }
                      }

                      if(Method=="One-Step"){
                        if(!is.null(CV)){
                          Pheno=Phenotype %>% filter(Env %in% c(Trial))
                          Pheno=Phenotype[,c("Genotype","Env",Trait)]
                          GBS<<-PGEandCV(Pheno,GDre,GT,GIre,CV)
                          mv(from = "GBS", to = paste0("GBS_1_CV_",Study),envir = globalenv())
                          save(list=paste0("GBS_1_CV_",Study),file=paste0("GBS_1_CV_",Study,".RData"))
                        }else{
                          Pheno=Phenotype %>% filter(Env %in% c(Trial))
                          Pheno=Phenotype[,c("Genotype","Env",Trait)]
                          GBS<<-PGandE(Pheno,GDre,GT,GIre)
                          mv(from = "GBS", to = paste0("GBS_1_",Study),envir = globalenv())
                          save(list=paste0("GBS_1_",Study),file=paste0("GBS_1_",Study,".RData"))
                        }
                      }
                        if(Method=="One-Step"){
                          if(Scheme=="K-Fold"){
                            Matrix<<-GE_Matrix_IND(genotypes=get(paste0("GBS_1_",Study))$geno, phenotype=get(paste0("GBS_1_",Study))$pheno,trait=Trait,Kernel=Kernel,GE=GE,model=GE_model,Sparse=Sparse,m=m,degree=degree, nL=nL)
                            mv(from = "Matrix", to = paste0("Matrix_",Study,GE_model),envir = globalenv())
                            save(list=paste0("Matrix_",Study),file=paste0("Matrix_",Study,".RData"))
                          }else{
                            Matrix<<-GE_Matrix_IND(train_genotypes=get(paste0("GBS_1_",Study))$geno, train_phenotype=get(paste0("GBS_1_",Study))$pheno,test_genotypes=get(paste0("GBS_1_",Study))$geno, test_phenotype=get(paste0("GBS_1_",Study))$pheno,trait=Trait,Kernel=Kernel,GE=GE,model=GE_model,Sparse=Sparse,m=m,degree=degree, nL=nL)
                            mv(from = "Matrix", to = paste0("Matrix_",Study,GE_model),envir = globalenv())
                            save(list=paste0("Matrix_",Study),file=paste0("Matrix_",Study,".RData"))
                          }

                        }

                    }
                    }


                    if(Method=="Two-Step"){
                    Results_All=list()
                    for(j in 1:length(Trait)){
                      if(Outcome=="Tested"){
                        if(Scheme=="K-Fold"){
                          if(is.null(GBS_Train)){
                            #GBS_Train=get(paste0("GBS_2_",Training,"_",Trait[i]))
                            GBS_Train=get(paste0("GBS_2_",Training))
                          }
                          if(Package=="rrBLUP"){
                            if(GAGS==TRUE){
                              if(!is.null(PC)){
                                #No CV

                                Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                      phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      Kernel=Kernel,
                                                                                                      PCA=GBS_Train$PC[,1:PC],
                                                                                                      CV=CV,
                                                                                                      markers=markers,
                                                                                                      folds = folds,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation,
                                                                                                      Y=GBS_Train$pheno,
                                                                                                      GM=GBS_Train$map,
                                                                                                      GD=GBS_Train$numeric,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total
                                )})
                              }else{
                                #No CV, No PC
                                Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                      phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      Kernel=Kernel,
                                                                                                      PCA=PC,
                                                                                                      CV=CV,
                                                                                                      markers=markers,
                                                                                                      folds = folds,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation,
                                                                                                      Y=GBS_Train$pheno,
                                                                                                      GM=GBS_Train$map,
                                                                                                      GD=GBS_Train$numeric,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total
                                )})
                              }
                            }else{
                              if(!is.null(CV)){
                                if(!is.null(PC)){
                                  #No GAGS
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=GBS_Train$PC[,1:PC],
                                                                                                   CV=GBS_Train$CV[,-1],
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=PC,
                                                                                                   CV=GBS_Train$CV[,-1],
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }
                              }else{
                                if(!is.null(PC)){
                                  #No GAGS, No CV
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=GBS_Train$PC[,1:PC],
                                                                                                   CV=CV,
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No CV, No PC
                                  if(Messages==TRUE){
                                    print(paste0("Peforming ",Method," ",Scheme," ",Type," using the package ",Package," with the model ",model," and ",Kernel," on ",Trait[j]," using ",Training,"."))
                                  }
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=PC,
                                                                                                   CV=CV,
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }
                              }

                            }


                          }

                          if(Package=="MAS"){
                            if(GAGS==TRUE){
                              if(!is.null(PC)){
                                #No CV
                                Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   PCA=GBS_Train$PC[,1:PC],
                                                                                                   CV=CV,
                                                                                                   folds = folds,
                                                                                                   transformation=transformation,
                                                                                                   Y=GBS_Train$pheno,
                                                                                                   GM=GBS_Train$map,
                                                                                                   GD=GBS_Train$numeric,
                                                                                                   GWAS=GWAS,
                                                                                                   alpha=alpha,
                                                                                                   threshold=threshold,
                                                                                                   QTN=QTN,
                                                                                                   PCA.total=PCA.total
                                )})
                              }else{
                                #No CV, No PC
                                Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   PCA=PC,
                                                                                                   CV=CV,
                                                                                                   folds = folds,
                                                                                                   transformation=transformation,
                                                                                                   Y=GBS_Train$pheno,
                                                                                                   GM=GBS_Train$map,
                                                                                                   GD=GBS_Train$numeric,
                                                                                                   GWAS=GWAS,
                                                                                                   alpha=alpha,
                                                                                                   threshold=threshold,
                                                                                                   QTN=QTN,
                                                                                                   PCA.total=PCA.total
                                )})
                              }
                            }else{
                              if(!is.null(CV)){
                                if(!is.null(PC)){
                                  #No GAGS
                                  Results=sapply(1:Replications, function(i,...){Results=MAS_CV(phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                PCA=GBS_Train$PC[,1:PC],
                                                                                                CV=GBS_Train$CV[,-1],
                                                                                                folds = folds,
                                                                                                transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=MAS_CV(phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                PCA=PC,
                                                                                                CV=GBS_Train$CV[,-1],
                                                                                                folds = folds,
                                                                                                transformation=transformation
                                  )})
                                }
                              }

                            }


                          }

                          if(Package=="BGLR"){
                            if(Type=="Ordinal"){
                              if(GAGS==TRUE){
                                if(!is.null(PC)){
                                  #No CV
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                              phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              PCA=GBS_Train$PC[,1:PC],
                                                                                                              CV=CV,
                                                                                                              markers=markers,
                                                                                                              folds = folds,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL,
                                                                                                              Y=GBS_Train$pheno,
                                                                                                              GM=GBS_Train$map,
                                                                                                              GD=GBS_Train$numeric,
                                                                                                              GWAS=GWAS,
                                                                                                              alpha=alpha,
                                                                                                              threshold=threshold,
                                                                                                              QTN=QTN,
                                                                                                              PCA.total=PCA.total
                                  )})
                                }else{
                                  #No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                              phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              PCA=PC,
                                                                                                              CV=CV,
                                                                                                              markers=markers,
                                                                                                              folds = folds,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL,
                                                                                                              Y=GBS_Train$pheno,
                                                                                                              GM=GBS_Train$map,
                                                                                                              GD=GBS_Train$numeric,
                                                                                                              GWAS=GWAS,
                                                                                                              alpha=alpha,
                                                                                                              threshold=threshold,
                                                                                                              QTN=QTN,
                                                                                                              PCA.total=PCA.total
                                  )})
                                }
                              }else{
                                if(!is.null(CV)){
                                  if(!is.null(PC)){
                                    #No GAGS
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = GBS_Train$geno,
                                                                                                           phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           PCA=GBS_Train$PC[,1:PC],
                                                                                                           CV=GBS_Train$CV[,-1],
                                                                                                           markers=markers,
                                                                                                           folds = folds,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }else{
                                    #No GAGS, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = GBS_Train$geno,
                                                                                                           phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           PCA=PC,
                                                                                                           CV=GBS_Train$CV[,-1],
                                                                                                           markers=markers,
                                                                                                           folds = folds,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }
                                }else{
                                  if(!is.null(PC)){
                                    #No GAGS, No CV
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = GBS_Train$geno,
                                                                                                           phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           PCA=GBS_Train$PC[,1:PC],
                                                                                                           CV=CV,
                                                                                                           markers=markers,
                                                                                                           folds = folds,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }else{
                                    #No GAGS, No CV, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_CV(genotypes = GBS_Train$geno,
                                                                                                           phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           PCA=PC,
                                                                                                           CV=CV,
                                                                                                           markers=markers,
                                                                                                           folds = folds,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }
                                }

                              }
                            }else{

                              if(GAGS==TRUE){
                                if(!is.null(PC)){
                                  #No CV
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                      phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      PCA=GBS_Train$PC[,1:PC],
                                                                                                      CV=CV,
                                                                                                      markers=markers,
                                                                                                      folds = folds,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation,
                                                                                                      Y=GBS_Train$pheno,
                                                                                                      GM=GBS_Train$map,
                                                                                                      GD=GBS_Train$numeric,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total
                                  )})
                                }else{
                                  #No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_CV(genotypes = GBS_Train$geno,
                                                                                                      phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      PCA=PC,
                                                                                                      CV=CV,
                                                                                                      markers=markers,
                                                                                                      folds = folds,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation,
                                                                                                      Y=GBS_Train$pheno,
                                                                                                      GM=GBS_Train$map,
                                                                                                      GD=GBS_Train$numeric,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total
                                  )})
                                }
                              }else{
                                if(!is.null(CV)){
                                  if(!is.null(PC)){
                                    #No GAGS
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=GBS_Train$PC[,1:PC],
                                                                                                   CV=GBS_Train$CV[,-1],
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }else{
                                    #No GAGS, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=PC,
                                                                                                   CV=GBS_Train$CV[,-1],
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }
                                }else{
                                  if(!is.null(PC)){
                                    #No GAGS, No CV
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=GBS_Train$PC[,1:PC],
                                                                                                   CV=CV,
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }else{
                                    #No GAGS, No CV, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   PCA=PC,
                                                                                                   CV=CV,
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }
                                }

                              }

                            }
                          }

                          if(Package=="GAPIT"){
                            if(!is.null(CV)){
                              #No GAGS
                              Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_CV(genotypes = GBS_Train$numeric,
                                                                                                 phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                 myGD=GBS_Train$map,
                                                                                                 model=model,
                                                                                                 PCA.total=PC,
                                                                                                 CV=GBS_Train$CV[,-1],
                                                                                                 kinship=kinship,
                                                                                                 markers=markers,
                                                                                                 folds = folds,
                                                                                                 transformation=transformation
                              )})

                            }else{
                              #No GAGS, No CV
                              Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_CV(genotypes = GBS_Train$numeric,
                                                                                                 phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                 myGD=GBS_Train$map,
                                                                                                 model=model,
                                                                                                 PCA.total=PC,
                                                                                                 CV=CV,
                                                                                                 kinship=kinship,
                                                                                                 markers=markers,
                                                                                                 folds = folds,
                                                                                                 transformation=transformation
                              )})

                            }
                          }

                          if(Package=="caret"){
                            Results=sapply(1:Replications, function(i,...){Results=Caret_Models_CV(genotypes = GBS_Train$geno,
                                                                                                   phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   type=type,
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   folds = folds,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation,
                                                                                                   sampling=sampling,
                                                                                                   repeats=repeats,
                                                                                                   method=method
                            )})
                          }

                          if(Package=="GLM"){
                            Results=sapply(1:Replications, function(i,...){Results=GLM_CV(genotypes = GBS_Train$geno,
                                                                                          phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                          fam=fam,
                                                                                          Kernel=Kernel,
                                                                                          markers=markers,
                                                                                          folds = folds,
                                                                                          Sparse=Sparse,
                                                                                          m=m,
                                                                                          degree=degree,
                                                                                          nL=nL
                            )})
                          }

                        }
                        if(Scheme=="VS"){
                          if(is.null(GBS_Train)){
                            #GBS_Train=get(paste0("GBS_2_",Training,"_",Trait[j]))
                            #GBS_Predict=get(paste0("GBS_2_",Prediction,"_",Trait[j]))
                            GBS_Train=get(paste0("GBS_2_",Training))
                            GBS_Predict=get(paste0("GBS_2_",Prediction))
                          }
                          if(Package=="rrBLUP"){
                            if(GAGS==TRUE){
                              if(!is.null(PC)){
                                #No CV
                                Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_GM=GBS_Train$map,
                                                                                                      train_GD=GBS_Train$numeric,
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total,
                                                                                                      transformation=transformation
                                )})

                              }else{
                                #No CV, No PC
                                Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_GM=GBS_Train$map,
                                                                                                      train_GD=GBS_Train$numeric,
                                                                                                      train_PCA=PC,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total,
                                                                                                      transformation=transformation
                                )})
                              }
                            }else{
                              if(!is.null(CV)){
                                if(!is.null(PC)){
                                  #No GAGS
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=GBS_Train$PC[,1:PC],
                                                                                                   train_CV=GBS_Train$CV[,-1],
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                   test_CV=GBS_Predict$CV[,-1],
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=PC,
                                                                                                   train_CV=GBS_Train$CV[,-1],
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=PC,
                                                                                                   test_CV=GBS_Predict$CV[,-1],
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }
                              }else{
                                if(!is.null(PC)){
                                  #No GAGS, No CV
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=GBS_Train$PC[,1:PC],
                                                                                                   train_CV=CV,
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                   test_CV=CV,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=PC,
                                                                                                   train_CV=CV,
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=PC,
                                                                                                   test_CV=CV,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                  )})
                                }
                              }

                            }


                          }

                          if(Package=="MAS"){
                            if(GAGS==TRUE){
                              if(!is.null(PC)){
                                #No CV
                                Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_GM=GBS_Train$map,
                                                                                                   train_GD=GBS_Train$numeric,
                                                                                                   train_PCA=GBS_Train$PC[,1:PC],
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                   markers=markers,
                                                                                                   GWAS=GWAS,
                                                                                                   alpha=alpha,
                                                                                                   threshold=threshold,
                                                                                                   QTN=QTN,
                                                                                                   PCA.total=PCA.total,
                                                                                                   transformation=transformation
                                )})
                              }else{
                                #No CV, No PC
                                Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_GM=GBS_Train$map,
                                                                                                   train_GD=GBS_Train$numeric,
                                                                                                   train_PCA=PC,
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=PC,
                                                                                                   markers=markers,
                                                                                                   GWAS=GWAS,
                                                                                                   alpha=alpha,
                                                                                                   threshold=threshold,
                                                                                                   QTN=QTN,
                                                                                                   PCA.total=PCA.total,
                                                                                                   transformation=transformation
                                )})
                              }
                            }else{
                              if(!is.null(CV)){
                                if(!is.null(PC)){
                                  #No GAGS
                                  Results=sapply(1:Replications, function(i,...){Results=MAS_VS(train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                train_PCA=GBS_Train$PC[,1:PC],
                                                                                                train_CV=GBS_Train$CV[,-1],
                                                                                                test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                test_CV=GBS_Predict$CV[,-1],
                                                                                                transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=MAS_VS(train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                train_PCA=PC,
                                                                                                train_CV=GBS_Train$CV[,-1],
                                                                                                test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                test_PCA=PC,
                                                                                                test_CV=GBS_Predict$CV[,-1],
                                                                                                transformation=transformation
                                  )})
                                }
                              }

                            }


                          }
                          if(Package=="BGLR"){
                            if(Type=="Ordinal"){
                              if(GAGS==TRUE){
                                if(!is.null(PC)){
                                  #No CV
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                              train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              train_GM=GBS_Train$map,
                                                                                                              train_GD=GBS_Train$numeric,
                                                                                                              train_PCA=GBS_Train$PC[,1:PC],
                                                                                                              test_genotypes= GBS_Predict$geno,
                                                                                                              test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                              test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              markers=markers,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL,
                                                                                                              GWAS=GWAS,
                                                                                                              alpha=alpha,
                                                                                                              threshold=threshold,
                                                                                                              QTN=QTN,
                                                                                                              PCA.total=PCA.total
                                  )})
                                }else{
                                  #No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                              train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              train_GM=GBS_Train$map,
                                                                                                              train_GD=GBS_Train$numeric,
                                                                                                              train_PCA=PC,
                                                                                                              test_genotypes= GBS_Predict$geno,
                                                                                                              test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                              test_PCA=PC,
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              markers=markers,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL,
                                                                                                              GWAS=GWAS,
                                                                                                              alpha=alpha,
                                                                                                              threshold=threshold,
                                                                                                              QTN=QTN,
                                                                                                              PCA.total=PCA.total
                                  )})
                                }
                              }else{
                                if(!is.null(CV)){
                                  if(!is.null(PC)){
                                    #No GAGS
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = GBS_Train$geno,
                                                                                                           train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           train_PCA=GBS_Train$PC[,1:PC],
                                                                                                           train_CV=GBS_Train$CV[,-1],
                                                                                                           test_genotypes= GBS_Predict$geno,
                                                                                                           test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                           test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                           test_CV=GBS_Predict$CV[,-1],
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           markers=markers,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }else{
                                    #No GAGS, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = GBS_Train$geno,
                                                                                                           train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           train_PCA=PC,
                                                                                                           train_CV=GBS_Train$CV[,-1],
                                                                                                           test_genotypes= GBS_Predict$geno,
                                                                                                           test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                           test_PCA=PC,
                                                                                                           test_CV=GBS_Predict$CV[,-1],
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           markers=markers,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }
                                }else{
                                  if(!is.null(PC)){
                                    #No GAGS, No CV
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = GBS_Train$geno,
                                                                                                           train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           train_PCA=GBS_Train$PC[,1:PC],
                                                                                                           train_CV=CV,
                                                                                                           test_genotypes= GBS_Predict$geno,
                                                                                                           test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                           test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                           test_CV=CV,
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           markers=markers,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }else{
                                    #No GAGS, No CV, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS(train_genotypes = GBS_Train$geno,
                                                                                                           train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                           train_PCA=PC,
                                                                                                           train_CV=CV,
                                                                                                           test_genotypes= GBS_Predict$geno,
                                                                                                           test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                           test_PCA=PC,
                                                                                                           test_CV=CV,
                                                                                                           model=model,
                                                                                                           Kernel=Kernel,
                                                                                                           markers=markers,
                                                                                                           nIter = nIter,
                                                                                                           burnIn = burnIn,
                                                                                                           Sparse=Sparse,
                                                                                                           m=m,
                                                                                                           degree=degree,
                                                                                                           nL=nL
                                    )})
                                  }
                                }

                              }
                            }else{

                              if(GAGS==TRUE){
                                if(!is.null(PC)){
                                  #No CV
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_GM=GBS_Train$map,
                                                                                                      train_GD=GBS_Train$numeric,
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total
                                  )})
                                }else{
                                  #No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_VS(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_GM=GBS_Train$map,
                                                                                                      train_GD=GBS_Train$numeric,
                                                                                                      train_PCA=PC,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total
                                  )})
                                }
                              }else{
                                if(!is.null(CV)){
                                  if(!is.null(PC)){
                                    #No GAGS
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=GBS_Train$PC[,1:PC],
                                                                                                   train_CV=GBS_Train$CV[,-1],
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                   test_CV=GBS_Predict$CV[,-1],
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }else{
                                    #No GAGS, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=PC,
                                                                                                   train_CV=GBS_Train$CV[,-1],
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=PC,
                                                                                                   test_CV=GBS_Predict$CV[,-1],
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }
                                }else{
                                  if(!is.null(PC)){
                                    #No GAGS, No CV
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=GBS_Train$PC[,1:PC],
                                                                                                   train_CV=CV,
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                   test_CV=CV,
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }else{
                                    #No GAGS, No CV, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=PC,
                                                                                                   train_CV=CV,
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=PC,
                                                                                                   test_CV=CV,
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   nIter = nIter,
                                                                                                   burnIn = burnIn,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation
                                    )})
                                  }
                                }

                              }

                            }
                          }
                          if(Package=="caret"){
                            Results=sapply(1:Replications, function(i,...){Results=Caret_Models_VS(train_genotypes = GBS_Train$geno,
                                                                                                   train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_genotypes= GBS_Predict$geno,
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   type=type,
                                                                                                   model=model,
                                                                                                   Kernel=Kernel,
                                                                                                   markers=markers,
                                                                                                   Sparse=Sparse,
                                                                                                   m=m,
                                                                                                   degree=degree,
                                                                                                   nL=nL,
                                                                                                   transformation=transformation,
                                                                                                   sampling=sampling,
                                                                                                   repeats=repeats,
                                                                                                   method=method
                            )})
                          }


                          if(Package=="GAPIT"){
                            if(!is.null(CV)){
                              #No GAGS
                              Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_VS(genotypes = GBS_Train$numeric,
                                                                                                 phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                 train_GM=GBS_Train$map,
                                                                                                 train_CV=GBS_Train$CV[,-1],
                                                                                                 test_genotypes= GBS_Predict$numeric,
                                                                                                 test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                 test_GM=GBS_Predict$map,
                                                                                                 test_CV=GBS_Predict$CV[,-1],
                                                                                                 model=model,
                                                                                                 PCA.total=PC,
                                                                                                 kinship=kinship,
                                                                                                 markers=markers,
                                                                                                 transformation=transformation
                              )})

                            }else{
                              #No GAGS, No CV
                              Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_VS(genotypes = GBS_Train$numeric,
                                                                                                 phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                 train_GM=GBS_Train$map,
                                                                                                 train_CV=CV,
                                                                                                 test_genotypes= GBS_Predict$numeric,
                                                                                                 test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                 test_GM=GBS_Predict$map,
                                                                                                 test_CV=CV,
                                                                                                 model=model,
                                                                                                 PCA.total=PC,
                                                                                                 kinship=kinship,
                                                                                                 markers=markers,
                                                                                                 transformation=transformation
                              )})

                            }
                          }
                          if(Package=="GLM"){
                            Results=sapply(1:Replications, function(i,...){Results=GLM_VS(train_genotypes = GBS_Train$geno,
                                                                                          train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                          test_genotypes= GBS_Predict$geno,
                                                                                          test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                          fam=fam,
                                                                                          Kernel=Kernel,
                                                                                          markers=markers,
                                                                                          Sparse=Sparse,
                                                                                          m=m,
                                                                                          degree=degree,
                                                                                          nL=nL
                            )})
                          }

                        }
                      }
                      if(Outcome=="Untested"){
                        if(Scheme=="K-Fold"){
                          if(is.null(GBS_Train)){
                            GBS_Train=get(paste0("GBS_2_",Training,"_",Trait[j]))
                          }
                          #Is there even a point for K-fold for untested?
                          #We can do this if we include blank values, but not in rrBLUP or DL
                          if(package=="rrBLUP"){
                            if(Kernel=="Markers"){
                              NULL
                            }else{}
                          }
                        }
                        if(Scheme=="VS"){
                          if(is.null(GBS_Train)){
                            #GBS_Train=get(paste0("GBS_2_Untested_",Training,"_",Trait[j]))
                            #GBS_Predict=get(paste0("GBS_2_Untested_",Prediction,"_",Trait[j]))
                            GBS_Train=get(paste0("GBS_2_Untested_",Training))
                            GBS_Predict=get(paste0("GBS_2_Untested_",Prediction))
                          }
                          if(Package=="rrBLUP"){
                            if(GAGS==TRUE){
                              if(!is.null(PC)){
                                #No CV
                                Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                         train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                         train_GM=GBS_Train$map,
                                                                                                         train_GD=GBS_Train$numeric,
                                                                                                         train_PCA=GBS_Train$PC[,1:PC],
                                                                                                         test_genotypes= GBS_Predict$geno,
                                                                                                         test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                         test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                         Kernel=Kernel,
                                                                                                         markers=markers,
                                                                                                         Sparse=Sparse,
                                                                                                         m=m,
                                                                                                         degree=degree,
                                                                                                         nL=nL,
                                                                                                         GWAS=GWAS,
                                                                                                         alpha=alpha,
                                                                                                         threshold=threshold,
                                                                                                         QTN=QTN,
                                                                                                         PCA.total=PCA.total,
                                                                                                         transformation=transformation
                                )})

                              }else{
                                #No CV, No PC
                                Results=sapply(1:Replications, function(i,...){Results=rrBLUP_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                         train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                         train_GM=GBS_Train$map,
                                                                                                         train_GD=GBS_Train$numeric,
                                                                                                         train_PCA=PC,
                                                                                                         test_genotypes= GBS_Predict$geno,
                                                                                                         test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                         test_PCA=PC,
                                                                                                         Kernel=Kernel,
                                                                                                         markers=markers,
                                                                                                         Sparse=Sparse,
                                                                                                         m=m,
                                                                                                         degree=degree,
                                                                                                         nL=nL,
                                                                                                         GWAS=GWAS,
                                                                                                         alpha=alpha,
                                                                                                         threshold=threshold,
                                                                                                         QTN=QTN,
                                                                                                         PCA.total=PCA.total,
                                                                                                         transformation=transformation
                                )})
                              }
                            }else{
                              if(!is.null(CV)){
                                if(!is.null(PC)){
                                  #No GAGS
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      train_CV=GBS_Train$CV[,-1],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      test_CV=GBS_Predict$CV[,-1],
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=PC,
                                                                                                      train_CV=GBS_Train$CV[,-1],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      test_CV=GBS_Predict$CV[,-1],
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                  )})
                                }
                              }else{
                                if(!is.null(PC)){
                                  #No GAGS, No CV
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      train_CV=CV,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      test_CV=CV,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                  )})
                                }else{
                                  if(Messages==TRUE){
                                    print(paste0("Peforming ",Method," ",Scheme," ",Type," using the package ",Package," with the model ",model," and ",Kernel," on ",Trait[j]," using ",Training,"."))
                                  }
                                  #No GAGS, No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=rrBLUP_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=PC,
                                                                                                      train_CV=CV,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      test_CV=CV,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                  )})
                                }
                              }

                            }


                          }

                          if(Package=="MAS"){
                            if(GAGS==TRUE){
                              if(!is.null(PC)){
                                #No CV
                                Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_GM=GBS_Train$map,
                                                                                                      train_GD=GBS_Train$numeric,
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      markers=markers,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total,
                                                                                                      transformation=transformation
                                )})
                              }else{
                                #No CV, No PC
                                Results=sapply(1:Replications, function(i,...){Results=MAS_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_GM=GBS_Train$map,
                                                                                                      train_GD=GBS_Train$numeric,
                                                                                                      train_PCA=PC,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      markers=markers,
                                                                                                      GWAS=GWAS,
                                                                                                      alpha=alpha,
                                                                                                      threshold=threshold,
                                                                                                      QTN=QTN,
                                                                                                      PCA.total=PCA.total,
                                                                                                      transformation=transformation
                                )})
                              }
                            }else{
                              if(!is.null(CV)){
                                if(!is.null(PC)){
                                  #No GAGS
                                  Results=sapply(1:Replications, function(i,...){Results=MAS_VS_UT(train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=GBS_Train$PC[,1:PC],
                                                                                                   train_CV=GBS_Train$CV[,-1],
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                   test_CV=GBS_Predict$CV[,-1],
                                                                                                   transformation=transformation
                                  )})
                                }else{
                                  #No GAGS, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=MAS_VS_UT(train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                   train_PCA=PC,
                                                                                                   train_CV=GBS_Train$CV[,-1],
                                                                                                   test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                   test_PCA=PC,
                                                                                                   test_CV=GBS_Predict$CV[,-1],
                                                                                                   transformation=transformation
                                  )})
                                }
                              }

                            }


                          }
                          if(Package=="BGLR"){
                            if(Type=="Ordinal"){
                              if(GAGS==TRUE){
                                if(!is.null(PC)){
                                  #No CV
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                                 train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                                 train_GM=GBS_Train$map,
                                                                                                                 train_GD=GBS_Train$numeric,
                                                                                                                 train_PCA=GBS_Train$PC[,1:PC],
                                                                                                                 test_genotypes= GBS_Predict$geno,
                                                                                                                 test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                                 test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                                 model=model,
                                                                                                                 Kernel=Kernel,
                                                                                                                 markers=markers,
                                                                                                                 nIter = nIter,
                                                                                                                 burnIn = burnIn,
                                                                                                                 Sparse=Sparse,
                                                                                                                 m=m,
                                                                                                                 degree=degree,
                                                                                                                 nL=nL,
                                                                                                                 GWAS=GWAS,
                                                                                                                 alpha=alpha,
                                                                                                                 threshold=threshold,
                                                                                                                 QTN=QTN,
                                                                                                                 PCA.total=PCA.total
                                  )})
                                }else{
                                  #No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                                 train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                                 train_GM=GBS_Train$map,
                                                                                                                 train_GD=GBS_Train$numeric,
                                                                                                                 train_PCA=PC,
                                                                                                                 test_genotypes= GBS_Predict$geno,
                                                                                                                 test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                                 test_PCA=PC,
                                                                                                                 model=model,
                                                                                                                 Kernel=Kernel,
                                                                                                                 markers=markers,
                                                                                                                 nIter = nIter,
                                                                                                                 burnIn = burnIn,
                                                                                                                 Sparse=Sparse,
                                                                                                                 m=m,
                                                                                                                 degree=degree,
                                                                                                                 nL=nL,
                                                                                                                 GWAS=GWAS,
                                                                                                                 alpha=alpha,
                                                                                                                 threshold=threshold,
                                                                                                                 QTN=QTN,
                                                                                                                 PCA.total=PCA.total
                                  )})
                                }
                              }else{
                                if(!is.null(CV)){
                                  if(!is.null(PC)){
                                    #No GAGS
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                              train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              train_PCA=GBS_Train$PC[,1:PC],
                                                                                                              train_CV=GBS_Train$CV[,-1],
                                                                                                              test_genotypes= GBS_Predict$geno,
                                                                                                              test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                              test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                              test_CV=GBS_Predict$CV[,-1],
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              markers=markers,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL
                                    )})
                                  }else{
                                    #No GAGS, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                              train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              train_PCA=PC,
                                                                                                              train_CV=GBS_Train$CV[,-1],
                                                                                                              test_genotypes= GBS_Predict$geno,
                                                                                                              test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                              test_PCA=PC,
                                                                                                              test_CV=GBS_Predict$CV[,-1],
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              markers=markers,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL
                                    )})
                                  }
                                }else{
                                  if(!is.null(PC)){
                                    #No GAGS, No CV
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                              train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              train_PCA=GBS_Train$PC[,1:PC],
                                                                                                              train_CV=CV,
                                                                                                              test_genotypes= GBS_Predict$geno,
                                                                                                              test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                              test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                              test_CV=CV,
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              markers=markers,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL
                                    )})
                                  }else{
                                    #No GAGS, No CV, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_Ordinal_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                              train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                              train_PCA=PC,
                                                                                                              train_CV=CV,
                                                                                                              test_genotypes= GBS_Predict$geno,
                                                                                                              test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                              test_PCA=PC,
                                                                                                              test_CV=CV,
                                                                                                              model=model,
                                                                                                              Kernel=Kernel,
                                                                                                              markers=markers,
                                                                                                              nIter = nIter,
                                                                                                              burnIn = burnIn,
                                                                                                              Sparse=Sparse,
                                                                                                              m=m,
                                                                                                              degree=degree,
                                                                                                              nL=nL
                                    )})
                                  }
                                }

                              }
                            }else{

                              if(GAGS==TRUE){
                                if(!is.null(PC)){
                                  #No CV
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                         train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                         train_GM=GBS_Train$map,
                                                                                                         train_GD=GBS_Train$numeric,
                                                                                                         train_PCA=GBS_Train$PC[,1:PC],
                                                                                                         test_genotypes= GBS_Predict$geno,
                                                                                                         test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                         test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                         model=model,
                                                                                                         Kernel=Kernel,
                                                                                                         markers=markers,
                                                                                                         nIter = nIter,
                                                                                                         burnIn = burnIn,
                                                                                                         Sparse=Sparse,
                                                                                                         m=m,
                                                                                                         degree=degree,
                                                                                                         nL=nL,
                                                                                                         transformation=transformation,
                                                                                                         GWAS=GWAS,
                                                                                                         alpha=alpha,
                                                                                                         threshold=threshold,
                                                                                                         QTN=QTN,
                                                                                                         PCA.total=PCA.total
                                  )})
                                }else{
                                  #No CV, No PC
                                  Results=sapply(1:Replications, function(i,...){Results=BGLR_GAGS_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                         train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                         train_GM=GBS_Train$map,
                                                                                                         train_GD=GBS_Train$numeric,
                                                                                                         train_PCA=PC,
                                                                                                         test_genotypes= GBS_Predict$geno,
                                                                                                         test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                         test_PCA=PC,
                                                                                                         model=model,
                                                                                                         Kernel=Kernel,
                                                                                                         markers=markers,
                                                                                                         nIter = nIter,
                                                                                                         burnIn = burnIn,
                                                                                                         Sparse=Sparse,
                                                                                                         m=m,
                                                                                                         degree=degree,
                                                                                                         nL=nL,
                                                                                                         transformation=transformation,
                                                                                                         GWAS=GWAS,
                                                                                                         alpha=alpha,
                                                                                                         threshold=threshold,
                                                                                                         QTN=QTN,
                                                                                                         PCA.total=PCA.total
                                  )})
                                }
                              }else{
                                if(!is.null(CV)){
                                  if(!is.null(PC)){
                                    #No GAGS
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      train_CV=GBS_Train$CV[,-1],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      test_CV=GBS_Predict$CV[,-1],
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                    )})
                                  }else{
                                    #No GAGS, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=PC,
                                                                                                      train_CV=GBS_Train$CV[,-1],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      test_CV=GBS_Predict$CV[,-1],
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                    )})
                                  }
                                }else{
                                  if(!is.null(PC)){
                                    #No GAGS, No CV
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=GBS_Train$PC[,1:PC],
                                                                                                      train_CV=CV,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=GBS_Predict$PC[,1:PC],
                                                                                                      test_CV=CV,
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                    )})
                                  }else{
                                    #No GAGS, No CV, No PC
                                    Results=sapply(1:Replications, function(i,...){Results=BGLR_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      train_PCA=PC,
                                                                                                      train_CV=CV,
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_PCA=PC,
                                                                                                      test_CV=CV,
                                                                                                      model=model,
                                                                                                      Kernel=Kernel,
                                                                                                      markers=markers,
                                                                                                      nIter = nIter,
                                                                                                      burnIn = burnIn,
                                                                                                      Sparse=Sparse,
                                                                                                      m=m,
                                                                                                      degree=degree,
                                                                                                      nL=nL,
                                                                                                      transformation=transformation
                                    )})
                                  }
                                }

                              }

                            }
                          }
                          if(Package=="caret"){
                            Results=sapply(1:Replications, function(i,...){Results=Caret_Models_CV_UT(train_genotypes = GBS_Train$geno,
                                                                                                      train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                      test_genotypes= GBS_Predict$geno,
                                                                                                      test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
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


                          if(Package=="GAPIT"){
                            if(!is.null(CV)){
                              #No GAGS
                              Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_VS_UT(genotypes = GBS_Train$numeric,
                                                                                                    phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                    train_GM=GBS_Train$map,
                                                                                                    train_CV=GBS_Train$CV[,-1],
                                                                                                    test_genotypes= GBS_Predict$numeric,
                                                                                                    test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                    test_GM=GBS_Predict$map,
                                                                                                    test_CV=GBS_Predict$CV[,-1],
                                                                                                    model=model,
                                                                                                    PCA.total=PC,
                                                                                                    kinship=kinship,
                                                                                                    markers=markers,
                                                                                                    transformation=transformation
                              )})

                            }else{
                              #No GAGS, No CV
                              Results=sapply(1:Replications, function(i,...){Results=GAPIT_GS_VS_UT(genotypes = GBS_Train$numeric,
                                                                                                    phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                                    train_GM=GBS_Train$map,
                                                                                                    train_CV=CV,
                                                                                                    test_genotypes= GBS_Predict$numeric,
                                                                                                    test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                                    test_GM=GBS_Predict$map,
                                                                                                    test_CV=CV,
                                                                                                    model=model,
                                                                                                    PCA.total=PC,
                                                                                                    kinship=kinship,
                                                                                                    markers=markers,
                                                                                                    transformation=transformation
                              )})

                            }
                          }
                          if(Package=="GLM"){
                            Results=sapply(1:Replications, function(i,...){Results=GLM_VS_UT(train_genotypes = GBS_Train$geno,
                                                                                             train_phenotype = GBS_Train$pheno[,c("Genotype",Trait[j])],
                                                                                             test_genotypes= GBS_Predict$geno,
                                                                                             test_phenotype= GBS_Predict$pheno[,c("Genotype",Trait[j])],
                                                                                             fam=fam,
                                                                                             Kernel=Kernel,
                                                                                             markers=markers,
                                                                                             Sparse=Sparse,
                                                                                             m=m,
                                                                                             degree=degree,
                                                                                             nL=nL
                            )})
                          }

                        }
                      }

                      if(!is.null(CV)){
                        CV_Message=TRUE
                      }else{
                        CV_Message=NULL
                      }

                      if(!is.null(PC)){
                        PCA_Message=TRUE
                      }else{
                        PCA_Message=NULL
                      }

                      if(GAGS==TRUE){
                        GAGS_Message=TRUE
                      }else{
                        GAGS_Message=NULL
                      }

                      if(Outcome=="Tested"){
                        Results_Accuracy=Extract_ACC(Results,Replications,Training,model,Kernel,CV_Message,Trait[j])

                        Results_Predictions=Extract_Pred(Results,Replications,Training,model,Kernel,CV_Message,Trait[j])

                        Results_Both=list(Accuracy=Results_Accuracy,Predictions=Results_Predictions)
                        #Results_All[[j]]=Results_Both
                        Results_All[[Trait[j]]]=Results_Both
                        #names(Results_All[[j]])<-Trait[j]
                      }

                      if(Outcome=="Untested"){
                        Results_Predictions=Extract_Pred_UT(Results,Replications,Prediction,model,Kernel,CV_Message,Trait[j])

                        Results_Both=list(Predictions=Results_Predictions)
                        #Results_All[[j]]=Results_Both
                        Results_All[[Trait[j]]]=Results_Both
                        #names(Results_All[[j]])<-Trait[j]
                      }

                    }
                    save(Results_All, file=paste0(Study,"_Results",".RData"))
                    return(Results_All)
                    }
                    if(Method=="One-Step"){
                      if(is.null(Matrix)){
                        Matrix=get(paste0("Matrix_",Study))
                      }

                      #We can do this with rrBLUP if we use Kernel/kin.blup
                      if(Outcome=="Tested"){
                        if(Scheme=="K-Fold"){
                          if(Package=="rrBLUP"){


                          }
                          if(Package=="BGLR"){
                            if(UN==TRUE){
                              Results=sapply(1:Replications, function(i,...){MTME(Matrix=Matrix,
                                                                                  trait=trait,
                                                                                  model=GE_model,
                                                                                  nIter = nIter,
                                                                                  burnIn = burnIn,
                                                                                  folds = folds,
                                                                                  UN=UN
                              )})
                            }else{
                              Results=sapply(1:Replications, function(i,...){MTME_UN(Matrix=Matrix,
                                                                                     trait=trait,
                                                                                     model=GE_model,
                                                                                     nIter = nIter,
                                                                                     burnIn = burnIn,
                                                                                     folds = folds
                              )})

                            }


                          }
                          if(Package=="tensorflow"){
                            if(UN==TRUE){
                              Results=sapply(1:Replications, function(i,...){MLP_CV(Matrix=Matrix,
                                                                                    trait=trait,
                                                                                    model=GE_model,
                                                                                    digits=digits,
                                                                                    nCVI=nCVI,
                                                                                    K=folds,
                                                                                    Sparse=Sparse,
                                                                                    folds = folds,
                                                                                    UN=UN
                              )})


                            }else{
                              Results=sapply(1:Replications, function(i,...){MLP_CV_UN(Matrix=Matrix,
                                                                                       trait=trait,
                                                                                       model=GE_model,
                                                                                       digits=digits,
                                                                                       nCVI=nCVI,
                                                                                       K=folds,
                                                                                       Sparse=Sparse,
                                                                                       folds = folds
                              )})

                            }


                          }


                        }


                        if(Scheme=="VS"){




                        }
                      }
                      if(Outcome=="Untested"){
                        if(Scheme=="K-Fold"){
                          #We can do this if we include blank values, but not in rrBLUP or DL
                          if(package=="rrBLUP"){
                            if(Kernel=="Markers"){


                            }else{


                            }
                          }
                        }
                        if(Scheme=="VS"){




                        }
                      }
                    }
                  }
}

Location_Processing <- function(Location,JD,Year,Indices=c("NDVI","SR","NDRE_1","NDRE_2","TCARI","SAVI","NWI1","NWI2","GNDVI","MTVI","ChlRE"),Band=c("900","975","BLU","GRN","RED","RE1","RE2","NIR"),Trials="All",Modifed=FALSE,Full=TRUE,Messages=TRUE)
{
  if(Messages==TRUE){
    print("Loading Packages")
  }
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(raster,imager,rgdal,foreign,utils,dplyr,FIELDimageR,stringr,doParallel,foreach,ggplot2,tictoc,tabularaster,gdata,svDialogs)
  #######################Commented out section is for Troubleshootin###########################
  #Location <- "Davenport"
  #JD <- 154
  #Date<-061721
  #Year <- 2021
  #Indices=c("NDVI","SR","NDRE_1","NDRE_2","TCARI","SAVI","NWI1","NWI2","GNDVI","MTVI","ChlRE")
  #Band=c("900","975","BLU","GRN","RED","RE1","RE2","NIR")
  #Trials="All"
  #i=1
  #Modifed=FALSE
  #Messages=TRUE
  #Full=TRUE

  #Let's get started
  wd=getwd()
  if(Messages==TRUE){
    print(paste0("You're working directory is set to: ",wd))
  }
  
  #Load in Target Shape file
  if(Messages==TRUE){
    print(paste0("Loading Target Shape file: ","Targets-",JD,".shp"))
  }
  Target <- readOGR(paste0("Targets-",JD,".shp"))
  
  #load in Band files
  if(Messages==TRUE){
    print(paste0("Bands Selected:",Band))
  }
  if(Messages==TRUE){
    if(Modifed==FALSE){
      print("Loading tif files")
    }else{
      print("Loading modified tif files")
    }
  }
  if(Modifed==FALSE){
    for(i in 1:length(Band)){
      if(Band[i]=="900"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","band900",".tif"), band = 2)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="975"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","band975",".tif"), band = 2)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="BLU"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group2",".tif"), band = 3)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="GRN"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group2",".tif"), band = 2)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="RED"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group2",".tif"), band = 1)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="RE1"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group1",".tif"), band = 1)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="RE2"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group1",".tif"), band = 2)
        mv(from = "band", to = Band[i])
      }
      if(Band[i]=="NIR"){
        band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group1",".tif"), band = 3)
        mv(from = "band", to = Band[i])
      }
    }
  }else{
    if(Band[i]=="900"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","band900","_modified.tif"), band = 2)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="975"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","band975","_modified.tif"), band = 2)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="BLU"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group2","_modified.tif"), band = 3)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="GRN"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group2","_modified.tif"), band = 2)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="RED"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group2","_modified.tif"), band = 1)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="RE1"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group1","_modified.tif"), band = 1)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="RE2"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group1","_modified.tif"), band = 2)
      mv(from = "band", to = Band[i])
    }
    if(Band[i]=="NIR"){
      band <- raster(paste0(Location,"-",JD,"_transparent_reflectance_","group1","_modified.tif"), band = 3)
      mv(from = "band", to = Band[i])
    }
  }

  #Set Raster Processing Options
  rasterOptions(maxmemory = 1e+11)
  rasterOptions(chunksize = 1e+11)  
  
  if(Messages==TRUE){
    print("Extract Raster Means Data from Overlying Targets")
  }
  # Extract Raster Means Data from Overlying Targets
  target.means.list=list()
  for(i in 1:length(Band)){
    target.mean <- extract(get(Band[i]), Target, fun = mean, na.rm=F, df=TRUE)/255
    target.mean$ID <- Target@data$id
    colnames(target.mean) <- c("ID",Band[i])
    target.mean <- rbind(target.mean, c(0,0))
    mv(from = "target.mean", to = paste0("target.mean.",Band[i]))
    target.means.list[[i]]=get(paste0("target.mean.",Band[i]))
  }
  #Merge Extracted Means
  target.means <- Reduce(function(...) merge(..., all=TRUE), target.means.list)  
  
  if(Messages==TRUE){
    print("Calibrating Values")
  }
  # Calibration Values
  ID <- c(0,1,2,3,4,5)
  Ref.Table <- data.frame(ID)
  for(i in 1:length(Band)){
    ID <- c(0,1,2,3,4,5)
    if(Band[i]=="BLU"){
      ref.BLU <- data.frame(ID,c(0,0.017,0.072,0.217,0.355,0.672))
      colnames(ref.BLU) <- c("ID","ref.BLU")
      Ref.Table=merge(Ref.Table,ref.BLU,by.x="ID",by.y="ID")
    }
    if(Band[i]=="GRN"){
      ref.GRN <- data.frame(ID,c(0,0.018,0.068,0.217,0.451,0.851))
      colnames(ref.GRN) <- c("ID","ref.GRN")
      Ref.Table=merge(Ref.Table,ref.GRN,by.x="ID",by.y="ID")
    }
    if(Band[i]=="RED"){
      ref.RED <- data.frame(ID,c(0,0.019,0.078,0.223,0.45,0.87))
      colnames(ref.RED) <- c("ID","ref.RED")
      Ref.Table=merge(Ref.Table,ref.RED,by.x="ID",by.y="ID")
    }
    if(Band[i]=="RE1"){
      ref.RE1 <- data.frame(ID,c(0,0.019,0.084,0.234,0.467,0.874))
      colnames(ref.RE1) <- c("ID","ref.RE1")
      Ref.Table=merge(Ref.Table,ref.RE1,by.x="ID",by.y="ID")
    }
    if(Band[i]=="RE2"){
      ref.RE2 <- data.frame(ID,c(0,0.02,0.095,0.244,0.482,0.875))
      colnames(ref.RE2) <- c("ID","ref.RE2")
      Ref.Table=merge(Ref.Table,ref.RE2,by.x="ID",by.y="ID")
    }
    if(Band[i]=="NIR"){
      ref.NIR <- data.frame(ID,c(0,0.02,0.099,0.25,0.495,0.872))
      colnames(ref.NIR) <- c("ID","ref.NIR")
      Ref.Table=merge(Ref.Table,ref.NIR,by.x="ID",by.y="ID")
    }
    if(Band[i]=="900"){
      ref.900 <- data.frame(ID,c(0,0.02,0.102,0.254,0.502,0.865))
      colnames(ref.900) <- c("ID","ref.900")
      Ref.Table=merge(Ref.Table,ref.900,by.x="ID",by.y="ID")
    }
    if(Band[i]=="975"){
      ref.975 <- data.frame(ID,c(0,0.021,0.104,0.258,0.503,0.845))
      colnames(ref.975) <- c("ID","ref.975")
      Ref.Table=merge(Ref.Table,ref.975,by.x="ID",by.y="ID")
    }
  }

  if(Messages==TRUE){
    print("Removing Reference Column 6")
  }
  #Fix Values
  Ref.Table <- Ref.Table[-c(6),]
  target.means <- target.means[-c(6),]

  if(Messages==TRUE){
    print("Calculate Regression")
  }
  #Calculate Regression
  for(i in 1:length(Band)){
    reg <- lm(Ref.Table[,paste0("ref.",Band[i])]~target.means[,Band[i]])
    int <- summary(reg)$coefficients[1]
    slp <- summary(reg)$coefficients[2]
    ad <- (get(Band[i])*slp-int)
    
    mv(from = "reg", to = paste0("reg.",Band[i]))
    mv(from = "int", to = paste0("int.",Band[i]))
    mv(from = "slp", to = paste0("slp.",Band[i]))
    mv(from = "ad", to = paste0("ad.",Band[i]))
  }
  if(Messages==TRUE){
    print("Identify Soil")
  }
  #Identify Soil
  NDVI.mask <- (ad.NIR-ad.RED)/(ad.NIR+ad.RED)
  m <- hist(NDVI.mask, breaks = 100)
  hista <- data.frame(counts= m$counts,breaks = m$mids)
  p=ggplot(hista, aes(x = breaks, y = counts, fill = counts)) + 
    geom_bar(stat = "identity",fill='black',alpha = 0.8)+
    scale_x_continuous(breaks = seq(-1,1,0.1),  ## without this you will get the same scale
                       labels = seq(-1,1,0.1))    ## as hist (question picture)
  #print(plot(p))
  
  png(paste0(Location,"_",JD,"_histogram",".png"))
  plot(p)
  dev.off()
  
  system(paste0('open ',Location,"_",JD,"_histogram",".png"))
  #Prompt
  user.input <- dlgInput("Enter a number to mask.", default = 0.1)$res
  mask <- NDVI.mask>as.numeric(user.input)
  #print(plot(mask))
  
  png(paste0(Location,"_",JD,"_mask",".png"))
  plot(mask)
  dev.off()
  
  system(paste0('open ',Location,"_",JD,"_mask",".png"))
  
  
  tic("Total")
  if(Messages==TRUE){
    print("Adjust Raster Layers")
  }
  # Adjust Raster Layers
  for(i in 1:length(Band)){
    #if(Band[i]=="900"){
      #adj.900 <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="975"){
      #adj.975 <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="BLU"){
      #adj.BLU <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="GRN"){
      #adj.GRN <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="RED"){
      #adj.RED <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="RE1"){
      #adj.RE1 <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="RE2"){
      #adj.RE2 <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    #if(Band[i]=="NIR"){
      #adj.NIR <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    #}
    
    
    adj <<- (get(Band[i])*get(paste0("slp.",Band[i]))-get(paste0("int.",Band[i])))*mask
    mv(from = "adj", to = paste0("adj.",Band[i]),envir = globalenv())
  }
  
  if(Messages==TRUE){
    print(paste0("Indices Selected:",Indices))
  }
  if(Messages==TRUE){
    print("Calculate Indices")
  }
  #Calculate Indices
  for(i in 1:length(Indices)){
    if(Indices[i]=="NDVI"){
      NDVI <<- ((adj.NIR-adj.RED)/(adj.NIR+adj.RED))*mask
    }
    if(Indices[i]=="SR"){
      SR <<- (adj.NIR/adj.RED)*mask
    }
    if(Indices[i]=="NDRE_1"){
      NDRE_1 <<- ((adj.NIR-adj.RE1)/(adj.NIR+adj.RE1))*mask
    }
    if(Indices[i]=="NDRE_2"){
      NDRE_2 <<- ((adj.NIR-adj.RE2)/(adj.NIR+adj.RE2))*mask
    }
    if(Indices[i]=="TCARI"){
      TCARI <<- (3*((adj.RE1-adj.RED)-0.2*(adj.RE1-adj.GRN)*(adj.RE1/adj.RED)))*mask
    }
    if(Indices[i]=="SAVI"){
      SAVI <<- (((adj.NIR-adj.RED)/(adj.NIR+adj.RED+0.5))*1.5)*mask
    }
    if(Indices[i]=="NWI1"){
      NWI1 <<- ((adj.975-adj.900)/(adj.975+adj.900))*mask
    }
    if(Indices[i]=="NWI2"){
      NWI2 <<- ((adj.975-adj.NIR)/(adj.975+adj.NIR))*mask
    }
    if(Indices[i]=="GNDVI"){
      GNDVI <<- ((adj.NIR-adj.GRN)/(adj.NIR+adj.GRN))*mask
    }
    if(Indices[i]=="MTVI"){
      MTVI <<- ((adj.RE1-adj.GRN)/(sqrt(((2*adj.NIR+1)^2)-(6*adj.NIR-5*sqrt(adj.RED))-0.5)))*mask
    }
    if(Indices[i]=="ChlRE"){
      ChlRE <<- ((adj.NIR/adj.RE1)-1)*mask
      }
  }
  base::save.image(file="Location_Processing.RData")
  if(Full==FALSE){stop("Part 1 has been completed and the function has stopped.")}
  
  if(Messages==TRUE){
    print(paste0("Trials selected: ", Trials))
  }
  #Identify Trials
  if(Trials=="All"){
    trial.names <- list.files(wd,pattern="*.shp", full.names=T, recursive=FALSE)
    trial.names <- trial.names[!str_detect(trial.names,pattern = "Target")]
  }else{
    trial.names <- list.files(wd,pattern="*.shp", full.names=T, recursive=FALSE)
    trial.names <- trial.names[str_detect(trial.names,pattern = Trials)]
  }
  base::save.image(file="Location_Processing.RData")
  #Define Cores
  UseCores <- detectCores() -1
  if(Messages==TRUE){
    print(paste0("Cores being used: ", UseCores))
  }
  
  
  #Register Core Cluster
  cl <- makeCluster(UseCores)

  clusterEvalQ(cl,{
    load(file = "Location_Processing.RData")
  })
  
  registerDoParallel(cl)
  
  #Set Capacity
  rasterOptions(maxmemory = 1e+11)
  rasterOptions(chunksize = 1e+11)
  
  if(Messages==TRUE){
    print("Extract Trial Data")
  }
  # Extract Trial Data
  tic("Data Extrtaction")
  foreach(i = 1:length(trial.names)) %dopar% {
    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(raster,imager,rgdal,foreign,utils,dplyr,FIELDimageR,stringr,doParallel,foreach,ggplot2,tictoc,tabularaster,gdata,svDialogs)
    
    mean.plot.list=list()
    plot.means <- NULL
    trials <- readOGR(trial.names[i])
    
    #Fix file Naming
    name <- basename(trial.names[i])
    name <- substr(name,1,nchar(name)-4)
    
    #Calculate Canopy Cover 
    plot.canopy <- fieldArea(mosaic = mask, fieldShape = trials, areaValue = 1)
    cc <- plot.canopy$areaPorcent
    can.cover <- select(cc, objArea) 
    id <- rownames(can.cover)
    can.cover <- cbind(id=id,can.cover)
    can.cover$objArea <- can.cover$objArea*100
    colnames(can.cover) <- c("Plot_ID",paste0("%Can.Cover_",JD))
    mean.plot.list[[i]]=can.cover
    
    for(i in 1:length(Indices)){
    mean.plot <- extract(get(Indices[i]), trials, fun = mean, na.rm=T, df=TRUE)
    mean.plot$ID <- trials@data$id
    colnames(mean.plot) <- c("Plot_ID",paste0(Indices[i],"_",JD))
    mv(from = "mean.plot", to = paste0("mean.plot.",Indices[i]))
    mean.plot.list[[i]]=get(paste0("mean.plot.",Indices[i]))
    }
    
    for(i in 1:length(Band)){
      mean.plot <- extract(get(paste0("adj.",Band[i])), trials, fun = mean, na.rm=T, df=TRUE)
      mean.plot$ID <- trials@data$id
      colnames(mean.plot) <- c("Plot_ID",paste0(Band[i],"_",JD))
      mv(from = "mean.plot", to = paste0("mean.plot.",Band[i]))
      mean.plot.list[[i]]=get(paste0("mean.plot.",Band[i]))
    }
    
    plot.means <- Reduce(function(...) merge(..., all=TRUE), mean.plot.list)
    #Add Columns for Location and Trial
    plot.means['Trial'] = name
    plot.means['Loc'] = Location
    plot.means['Year'] = Year
    plot.means <- plot.means %>% relocate(Trial, Loc, Year, .after = Plot_ID )
    
    #Create CSV
    write.csv(plot.means,paste0(Location,"_",name,"_",JD,".csv"), row.names = FALSE)
  }
  if(Messages==TRUE){
    print("Stopping Cluster")
  }
  stopCluster(cl)
  
  if(Messages==TRUE){
    print("Merging and Saving Files")
  }
  #Merge CSV Files
  all_trials <- list.files(pattern = paste0("*_",JD,".csv")) %>%
    lapply(read.csv) %>%
    bind_rows()
  
  #Save File for Location
  write.csv(all_trials,paste0(Location,"_",Year,"_",JD,".csv"),row.names = FALSE)
  tic("Save Files")
  
  # Save Adjusted Rasters
  for(i in 1:length(Band)){
    writeRaster(get(paste0("adj.",Band[i])), filename= file.path(wd,paste0(Band[i],"_",JD,".tif")), format="GTiff", overwrite=TRUE )
  }
  # Save Index Rasters
  for(i in 1:length(Indices)){
    writeRaster(get(Indices[i]), filename= file.path(wd,paste0(Indices[i],"_",JD,".tif")), format="GTiff", overwrite=TRUE )
  }
  toc()
  toc()
}


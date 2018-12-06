library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)



#### A function that defines CV blocks and folds (using BlockCV package) and then fits BRT model (from dismo package) on pre-defined folds, and saves outputs



blockCV_brt <- function(data = datcombo, pred.stack = pred_abs_2011, seed = 1222, blocking.iterations = 250, pred.variables, autocorr.variables = NULL,output.folder="D://CHID regional Alberta BRT/BRT_outputs/CV_output/",keep.out = TRUE, tc=3,lr=0.001,bf=0.5){ 
  # Arguments for this function
  ## data: data.frame object containing data for model fitting
  ## pred.stack: the raster stack/brick used as prediction dataset
  ## seed: random seed for blocking
  ## blocking.iterations: iterations for spatiaclBlock call
  ## pred.variables: a character vector giving the names of predictor variables that will be included in BRT models
  ## autocorr.variables: a character vector giving the names of predictor variables that will be used for autocorrelation assessment that will inform block size. defaults to 'pred.variables'. user can specify a subset by excluding variables, as in pred.variables[-c(57:63)]
  ## output.folder: path to output folder (if folder does not exist it is created)
  ## keep.BRT: logical, whether to keep the outputs from 'spatialBlocks' and 'gbm.step' functions in the workspace. both are automatically saved in the output folder regardless
  ## tc: BRT tree complexity
  ## lr: BRT learning rate
  ## bf: BRT bag fraction
  
  library(gstat)
  library(sp)
  library(blockCV)
  library(dismo)
  
  # resolve autocorrelation variables
  if (is.null(autocorr.variables) == TRUE) {
    autocorr.variables <- pred.variables
  }
  
  # convert data into SpatialPointDataFrame class
  datcombo_sp <-
    SpatialPointsDataFrame(coords = datcombo[, 33:34],
                           data = data,
                           proj4string = LCC) # where columns 33:34 are the X and Y coordinates
  
  # create a column converting abundance to occupancy
  datcombo_sp$OCCU <- 0 
  datcombo_sp$OCCU[which(datcombo_sp$ABUND > 0)] <- 1
  
  # create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
  ind <- which(names(pred.stack) %in% autocorr.variables)
  
  # Calculate autocorrelation range
  start_time <- Sys.time()
  sp.auto.arr1 <- spatialAutoRange(pred.stack[[ind]])
  end_time <- Sys.time()
  end_time - start_time
  
  # Use spatial blocks to separate train and test folds
  set.seed(seed)
  sp_blocks <-
    spatialBlock(
      speciesData = datcombo_sp,
      species = "OCCU",
      rasterLayer = pred.stack[[1]],
      iteration = blocking.iterations,
      theRange = sp.auto.arr1$range,
      selection = "random"
    )
  
  # fit BRT models using pre-determined folds for CV
  x1 <-
    try(brt1 <-
          gbm.step(
            datcombo,
            gbm.y = "ABUND",
            gbm.x = pred.variables,
            fold.vector = sp_blocks$foldID,
            n.folds = sp_blocks$k,
            family = "poisson",
            tree.complexity = tc,
            learning.rate = lr,
            bag.fraction = bf,
            offset = datcombo$logoffset,
            site.weights = datcombo$wt,
            keep.fold.models = T,
            keep.fold.fit = T
          ))
  
  # Define/creat folders for storing outputs
  if (class(x1) != "try-error") {
    z <- output.folder
    
    if (file.exists(z) == FALSE) {
      dir.create(z)
    }
    
    save(sp_blocks, file = paste(z, speclist[j], "blocks", sep = ""))
    save(brt1, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    
    ## Model evaluation
    varimp <- as.data.frame(brt1$contributions)
    write.csv(varimp, file = paste(z, speclist[j], "varimp.csv", sep = ""))
    cvstats <- t(as.data.frame(brt1$cv.statistics))
    write.csv(cvstats, file = paste(z, speclist[j], "cvstats.csv", sep =
                                      ""))
    pdf(paste(z, speclist[j], "_plot.pdf", sep = ""))
    gbm.plot(
      brt1,
      n.plots = 32,
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = T
    )
    dev.off()
    pdf(paste(z, speclist[j], "_plot.var-scale.pdf", sep = ""))
    gbm.plot(
      brt1,
      n.plots = 32,
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = F,
      write.title = F
    )
    dev.off()
    
    
    ## Model prediction
    
    rast <-
      predict(pred.stack,
              brt1,
              type = "response",
              n.trees = brt1$n.trees)
    writeRaster(
      rast,
      filename = paste(z, speclist[j], "_pred1km", sep = ""),
      format = "GTiff",
      overwrite = TRUE
    )
    png(paste(z, speclist[j], "_pred1km.png", sep = ""))
    plot(rast, zlim = c(0, 1))
    points(datcombo_sp$X, datcombo_sp$Y, cex = 0.05)
    dev.off()
    
    writeOGR(
      datcombo_sp,
      dsn = paste(z, "surveypoints.shp", sep = ""),
      layer = "datcombo_sp",
      driver = "ESRI Shapefile"
    )
    
    out<-list(brt1, sp_blocks)
    if(keep.out==T) 
      out
  }
}

#### load data and prepare objects
load("D:/CHID regional Alberta BRT/AB_BRT_Rproject/data_pack.RData")

j<-which(speclist=="CAWA") 
specoff <- offl[offl$SPECIES==as.character(speclist[j]),]
specdat2001 <- ABPC2001[ABPC2001$SPECIES == as.character(speclist[j]),] #n=423 for CAWA
dat1 <- right_join(specdat2001[,c(1:4,10)],survey2001[,1:3],by=c("SS","PCODE","PKEY")) #n=31864
dat1$SPECIES <- as.character(speclist[j])
dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
d2001 <- left_join(s2001, dat2001, by=c("SS","PCODE")) 
specdat2011 <- ABPC2011[ABPC2011$SPECIES == as.character(speclist[j]),] #n=989 for CAWA
dat1 <- right_join(specdat2011[,c(1:4,10)],survey2011[,1:3],by=c("SS","PCODE","PKEY")) #n=26314
dat1$SPECIES <- as.character(speclist[j])
dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
s2011 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
d2011 <- left_join(s2011, dat2011, by=c("SS","PCODE")) 

datcombo <- rbind(d2001,d2011)
datcombo$eco <- as.factor(datcombo$eco)

rm(list=setdiff(ls(),c("datcombo","pred_abs_2011","w","LCC","speclist","randomCV_brt","j","blockCV_brt")))
gc()


# Running blockCV_brt function for full model

pred.variables<-c( 
  "Species_Abie_Bal_v1",
  "Species_Abie_Bal_v1_Gauss250",
  "Species_Abie_Bal_v1_Gauss750",
  "Species_Betu_Pap_v1",
  "Species_Betu_Pap_v1_Gauss250",
  "Species_Betu_Pap_v1_Gauss750",
  "Species_Lari_Lar_v1",
  "Species_Lari_Lar_v1_Gauss250",
  "Species_Lari_Lar_v1_Gauss750",
  "Species_Pice_Gla_v1",
  "Species_Pice_Gla_v1_Gauss250",
  "Species_Pice_Gla_v1_Gauss750",
  "Species_Pice_Mar_v1",
  "Species_Pice_Mar_v1_Gauss250",
  "Species_Pice_Mar_v1_Gauss750",
  "Species_Pinu_Ban_v1",
  "Species_Pinu_Ban_v1_Gauss250",
  "Species_Pinu_Ban_v1_Gauss750",
  "Species_Pinu_Con_v1",
  "Species_Pinu_Con_v1_Gauss250",
  "Species_Pinu_Con_v1_Gauss750",
  "Species_Popu_Bal_v1",
  "Species_Popu_Bal_v1_Gauss250",
  "Species_Popu_Bal_v1_Gauss750",
  "Species_Popu_Tre_v1",
  "Species_Popu_Tre_v1_Gauss250",
  "Species_Popu_Tre_v1_Gauss750",
  "Species_Pseu_Men_v1",
  "Species_Pseu_Men_v1_Gauss250",
  "Species_Pseu_Men_v1_Gauss750",
  "Species_Thuj_Pli_v1",
  "Species_Thuj_Pli_v1_Gauss250",
  "Species_Thuj_Pli_v1_Gauss750",
  "Species_Tsug_Het_v1",
  "Species_Tsug_Het_v1_Gauss250",
  "Species_Tsug_Het_v1_Gauss750",
  "Structure_Biomass_TotalLiveAboveGround_v1",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  "Structure_Stand_Age_v1",
  "Structure_Stand_Age_v1_Gauss250",
  "Structure_Stand_Age_v1_Gauss750",
  "cti250",
  "cti250_Gauss250",
  "cti250_Gauss750",
  "CultivationCrop",
  "IndustrialSiteRural",
  "MineSite",
  "Pipeline",
  "RoadHardSurface",
  "SeismicLineNarrow",
  "SeismicLineWide",
  "TransmissionLine",
  "AHM",
  "DD18",
  "MAT",
  "MAP",
  "FFP",
  "MWMT",
  "MCMT"
)

start_time<-Sys.time()
brt1<- blockCV_brt(datcombo,pred.variables = pred.variables,autocorr.variables=pred.variables[-c(54:60)])  
end_time<-Sys.time()
end_time-start_time




# reordering varimp according to variable name (without Gauss filer) and relative influence
varimp$ordervar<-gsub("_Gauss500","",varimp$var)
varimp$ordervar<-gsub("_Gauss250","",varimp$ordervar)
varimp[order(varimp$ordervar,-varimp$rel.inf),]





# Assessing sampling representativeness ### OUTDATED!!!
load("D:/CHID regional Alberta BRT/BRT_outputs/model1/CAWAbrtAB.R") # load example brt
brt1$var.names
mess1<-mess(abs2011_1km[[brt1$var.names]], extract(abs2011_1km[[brt1$var.names]],cbind(datcombo$X,datcombo$Y)))

writeRaster(mess1, filename="mess_brt_preds.tif", format="GTiff",overwrite=TRUE)



# mean and coefficient of variation among the 20 brt models ### OUTDATED!!!
load("D:/CHID regional Alberta BRT/BRT_outputs/model1/CAWAbrtAB.R")
brt_preds<-stack(raster("D:/CHID regional Alberta BRT/BRT_outputs/model1/CAWA_pred1km.tif"))
names(brt_preds)<-"model1"
for(i in 2:20){
  brt_preds <- addLayer(brt_preds, raster(paste("D:/CHID regional Alberta BRT/BRT_outputs/model",i,"/CAWA_pred1km.tif",sep="")))
  names(brt_preds)[[i]] <- paste("model",i,sep="")
}

brt_preds<-addLayer(brt_preds,overlay(brt_preds,fun=mean, na.rm=T))
names(brt_preds[[21]])<-"mean"

brt_preds<-addLayer(brt_preds,overlay(brt_preds,fun=sd, na.rm=T))
names(brt_preds[[22]])<-"SD"

brt_preds<-addLayer(brt_preds,overlay(brt_preds,fun=function(x) sd(x, na.rm=T)/mean(x,na.rm=T)))
names(brt_preds[[23]])<-"CV"


setwd("D:/CHID regional Alberta BRT/BRT_outputs")
writeRaster(brt_preds$mean, filename="mean_brt_preds.tif", format="GTiff",overwrite=TRUE)
writeRaster(brt_preds$SD, filename="SD_brt_preds.tif", format="GTiff",overwrite=TRUE)
writeRaster(brt_preds$CV, filename="CV_brt_preds.tif", format="GTiff",overwrite=TRUE)

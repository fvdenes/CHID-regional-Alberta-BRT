library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(blockCV)
library(gbm)

# A function that fits the BRT model ('gbm.step' from dismo package) on pre-defined folds, and saves outputs ####
brt_blocks <- function(data = datcombo, pred.stack = pred_abs_2011, seed = 1222, pred.variables ,output.folder, blocks=NULL, keep.out = TRUE, tc=3,lr=0.001,bf=0.5, save.points.shp=FALSE){ 
  # Arguments for this function
  ## data: data.frame object containing data for model fitting
  ## pred.stack: the raster stack/brick used as prediction dataset
  ## pred.variables: a character vector giving the names of predictor variables that will be included in BRT models
  ## blocks: object resulting from 'spatialBlocks' function that contains classification of sample points into folds
  ## output.folder: path to output folder (if folder does not exist it is created)
  ## keep.out: logical, whether to keep the output in the workspace. both blocks object and the brt output are automatically saved in the output folder regardless
  ## tc: BRT tree complexity
  ## lr: BRT learning rate
  ## bf: BRT bag fraction
  ## save.points.shp: logical, whether to save survey points as a shapefile in output folder
  
  # fit BRT models using pre-determined folds for CV
  if (is.null(blocks)){
    folds<-NULL
    n.folds<-10
  }
  
  else {
    folds<-blocks$foldID
    n.folds<-blocks$k
  }
  
  x1 <-
    try(brt <-
          gbm.step(
            datcombo,
            gbm.y = "ABUND",
            gbm.x = pred.variables,
            fold.vector = folds,
            n.folds = n.folds,
            family = "poisson",
            tree.complexity = tc,
            learning.rate = lr,
            bag.fraction = bf,
            offset = datcombo$logoffset,
            site.weights = datcombo$wt,
            keep.fold.models = T,
            keep.fold.fit = T
          ))
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/10,
              bag.fraction = bf,
              offset = datcombo$logoffset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/100,
              bag.fraction = bf,
              offset = datcombo$logoffset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){
    stop("Restart model with even smaller learning rate, or other predictors!")
  }
  
  # Define/create folders for storing outputs
  if (class(x1) != "try-error") {
    z <- output.folder
    
    if (file.exists(z) == FALSE) {
      dir.create(z)
    }
    
    if (is.null(blocks)){
      save(brt, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    }
    
    else {
      save(blocks, file = paste(z, speclist[j], "blocks.R", sep = ""))
      save(brt, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    }  
    
    
    ## Model evaluation
    varimp <- as.data.frame(brt$contributions)
    write.csv(varimp, file = paste(z, speclist[j], "varimp.csv", sep = ""))
    cvstats <- t(as.data.frame(brt$cv.statistics))
    write.csv(cvstats, file = paste(z, speclist[j], "cvstats.csv", sep =
                                      ""))
    pdf(paste(z, speclist[j], "_plot.pdf", sep = ""))
    gbm.plot(
      brt,
      n.plots = length(pred.variables),
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = T
    )
    dev.off()
    pdf(paste(z, speclist[j], "_plot.var-scale.pdf", sep = ""))
    gbm.plot(
      brt,
      n.plots = length(pred.variables),
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = F,
      write.title = F
    )
    dev.off()
    
    
    ## Model prediction
    
    rast <-
      predict(pred.stack,
              brt,
              type = "response",
              n.trees = brt$n.trees)
    writeRaster(
      rast,
      filename = paste(z, speclist[j], "_pred1km", sep = ""),
      format = "GTiff",
      overwrite = TRUE
    )
    
    data_sp <-SpatialPointsDataFrame(coords = data[, 34:35], data = data, proj4string = LCC)
    png(paste(z, speclist[j], "_pred1km.png", sep = ""))
    plot(rast, zlim = c(0, 1))
    points(data_sp$X, data_sp$Y, cex = 0.05)
    dev.off()
    
    if(save.points.shp==T){
      writeOGR(
        data_sp,
        dsn = paste(z, "surveypoints.shp", sep = ""),
        layer = "data_sp",
        driver = "ESRI Shapefile"
      )
    }
    
    if(keep.out==T) {return(brt)}
  }
}


#load data and prepare objects ####
load("D:/CHID regional Alberta BRT/AB_BRT_Rproject/data_pack.RData")

#pred_abs_2011<-brick("D:/CHID regional Alberta BRT/prediction dataset/abs2011_250m.grd")


j<-which(speclist=="CAWA") 

specoff <- filter(offl, SPECIES==as.character(speclist[j]))
specoff <- distinct(specoff)

specdat2001 <- filter(ABPC2001, SPECIES == as.character(speclist[j]))
specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
names(specdat2001x)[3] <- "ABUND"
dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) #n=31864
dat1$SPECIES <- as.character(speclist[j])
dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
d2001 <- left_join(s2001, dat2001, by=c("SS")) 

specdat2011 <- filter(ABPC2011, SPECIES == as.character(speclist[j]))
specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
names(specdat2011x)[3] <- "ABUND"  
dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
dat2$SPECIES <- as.character(speclist[j])
dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
s2011 <- left_join(dat2,specoff, by=c("SPECIES","PKEY"))
d2011 <- left_join(s2011, dat2011, by=c("SS")) 


datcombo <- rbind(d2001,d2011)
datcombo$eco <- as.factor(datcombo$eco)



rm(list=setdiff(ls(),c("datcombo","datcombo_sp","pred_abs_2011","w","LCC","speclist","randomCV_brt","j","blockCV_brt","brt_blocks")))
gc()

# convert data into SpatialPointDataFrame class
datcombo_sp <-SpatialPointsDataFrame(coords = datcombo[, 34:35], data = datcombo, proj4string = LCC)
# create a column converting abundance to occupancy
datcombo_sp$OCCU <- 0 
datcombo_sp$OCCU[which(datcombo_sp$ABUND > 0)] <- 1

# Full model ####
# list variables for full model
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
  "watAB",
  "watAB_Gauss250",
  "watAB_Gauss750"
)





# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, etc).
ind <- which(names(pred_abs_2011) %in% pred.variables[-55])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr1 <- spatialAutoRange(pred_abs_2011[[ind]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_full <-  spatialBlock(
                      speciesData = datcombo_sp,
                      species = "OCCU",
                      rasterLayer = pred_abs_2011[[1]],
                      iteration = 250,
                      theRange = sp.auto.arr1$range,
                      selection = "random",
                      maskBySpecies = FALSE
  )
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt1<- brt_blocks(data=datcombo,pred.variables = pred.variables, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/full_model/", blocks=sp_block_full, save.points.shp = TRUE, lr=0.01)  
end_time<-Sys.time()
end_time-start_time

# Simplified model (selecting most influential scales) ####
pred.variables2<-c( 
  #"Species_Abie_Bal_v1",
  "Species_Abie_Bal_v1_Gauss250",
  #"Species_Abie_Bal_v1_Gauss750",
  #"Species_Betu_Pap_v1",
  #"Species_Betu_Pap_v1_Gauss250",
  "Species_Betu_Pap_v1_Gauss750",
  #"Species_Lari_Lar_v1",
  "Species_Lari_Lar_v1_Gauss250",
  #"Species_Lari_Lar_v1_Gauss750",
  #"Species_Pice_Gla_v1",
  "Species_Pice_Gla_v1_Gauss250",
  #"Species_Pice_Gla_v1_Gauss750",
  #"Species_Pice_Mar_v1",
  #"Species_Pice_Mar_v1_Gauss250",
  "Species_Pice_Mar_v1_Gauss750",
  #"Species_Pinu_Ban_v1",
  #"Species_Pinu_Ban_v1_Gauss250",
  "Species_Pinu_Ban_v1_Gauss750",
  #"Species_Pinu_Con_v1",
  #"Species_Pinu_Con_v1_Gauss250",
  "Species_Pinu_Con_v1_Gauss750",
  #"Species_Popu_Bal_v1",
  #"Species_Popu_Bal_v1_Gauss250",
  "Species_Popu_Bal_v1_Gauss750",
  #"Species_Popu_Tre_v1",
  "Species_Popu_Tre_v1_Gauss250",
  #"Species_Popu_Tre_v1_Gauss750",
  "Species_Pseu_Men_v1",
  #"Species_Pseu_Men_v1_Gauss250",
  #"Species_Pseu_Men_v1_Gauss750",
  "Species_Thuj_Pli_v1",
  #"Species_Thuj_Pli_v1_Gauss250",
  #"Species_Thuj_Pli_v1_Gauss750",
  "Species_Tsug_Het_v1",
  #"Species_Tsug_Het_v1_Gauss250",
  #"Species_Tsug_Het_v1_Gauss750",
  #"Structure_Biomass_TotalLiveAboveGround_v1",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  #"Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
  "Structure_Stand_Age_v1_Gauss750",
  #"cti250",
  "cti250_Gauss250",
  #"cti250_Gauss750",
  "CultivationCrop",
  "IndustrialSiteRural",
  "MineSite",
  "Pipeline",
  "RoadHardSurface",
  "SeismicLineNarrow",
  "SeismicLineWide",
  "TransmissionLine",
  #"watAB",
  #"watAB_Gauss250",
  "watAB_Gauss750"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind2 <- which(names(pred_abs_2011) %in% pred.variables2)

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr2 <- spatialAutoRange(pred_abs_2011[[ind2]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(142)
sp_block_select <- spatialBlock(
                      speciesData = datcombo_sp,
                      species = "OCCU",
                      rasterLayer = pred_abs_2011[[1]],
                      iteration = 250,
                      k=5,
                      theRange = sp.auto.arr2$range,
                      selection = "random",
                      maskBySpecies = FALSE
  )
end_time <- Sys.time()
end_time - start_time


start_time<-Sys.time()
brt2<- brt_blocks(data=datcombo,pred.variables = pred.variables2, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/selected_scales/", blocks=sp_block_select, save.points.shp = F, lr=0.01)  
end_time<-Sys.time()
end_time-start_time

# Cell level (250m, no Gaussian filter) ####

pred.variables3<-c( 
  "Species_Abie_Bal_v1",
  #"Species_Abie_Bal_v1_Gauss250",
  #"Species_Abie_Bal_v1_Gauss750",
  "Species_Betu_Pap_v1",
  #"Species_Betu_Pap_v1_Gauss250",
  #"Species_Betu_Pap_v1_Gauss750",
  "Species_Lari_Lar_v1",
  #"Species_Lari_Lar_v1_Gauss250",
  #"Species_Lari_Lar_v1_Gauss750",
  "Species_Pice_Gla_v1",
  #"Species_Pice_Gla_v1_Gauss250",
  #"Species_Pice_Gla_v1_Gauss750",
  "Species_Pice_Mar_v1",
  #"Species_Pice_Mar_v1_Gauss250",
  #"Species_Pice_Mar_v1_Gauss750",
  "Species_Pinu_Ban_v1",
  #"Species_Pinu_Ban_v1_Gauss250",
  #"Species_Pinu_Ban_v1_Gauss750",
  "Species_Pinu_Con_v1",
  #"Species_Pinu_Con_v1_Gauss250",
  #"Species_Pinu_Con_v1_Gauss750",
  "Species_Popu_Bal_v1",
  #"Species_Popu_Bal_v1_Gauss250",
  #"Species_Popu_Bal_v1_Gauss750",
  "Species_Popu_Tre_v1",
  #"Species_Popu_Tre_v1_Gauss250",
  #"Species_Popu_Tre_v1_Gauss750",
  "Species_Pseu_Men_v1",
  #"Species_Pseu_Men_v1_Gauss250",
  #"Species_Pseu_Men_v1_Gauss750",
  "Species_Thuj_Pli_v1",
  #"Species_Thuj_Pli_v1_Gauss250",
  #"Species_Thuj_Pli_v1_Gauss750",
  "Species_Tsug_Het_v1",
  #"Species_Tsug_Het_v1_Gauss250",
  #"Species_Tsug_Het_v1_Gauss750",
  "Structure_Biomass_TotalLiveAboveGround_v1",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  "Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
  #"Structure_Stand_Age_v1_Gauss750",
  "cti250",
  #"cti250_Gauss250",
  #"cti250_Gauss750",
  "CultivationCrop",
  "IndustrialSiteRural",
  "MineSite",
  "Pipeline",
  "RoadHardSurface",
  "SeismicLineNarrow",
  "SeismicLineWide",
  "TransmissionLine",
  "watAB"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind3 <- which(names(pred_abs_2011) %in% pred.variables3[-24])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr3 <- spatialAutoRange(pred_abs_2011[[ind3]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_cell <- spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr3$range,
  selection = "random",
  maskBySpecies = FALSE
)
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt3<- brt_blocks(data=datcombo,pred.variables = pred.variables3, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/cell_level/", blocks=sp_block_cell, save.points.shp = F, lr=0.01)  
end_time<-Sys.time()
end_time-start_time
 
save.image("D:/CHID regional Alberta BRT/BRT_outputs/outputs.RData")

# Only Gaussian filter with sigma = 250m ####
pred.variables4<-c( 
  #"Species_Abie_Bal_v1",
  "Species_Abie_Bal_v1_Gauss250",
  #"Species_Abie_Bal_v1_Gauss750",
  #"Species_Betu_Pap_v1",
  "Species_Betu_Pap_v1_Gauss250",
  #"Species_Betu_Pap_v1_Gauss750",
  #"Species_Lari_Lar_v1",
  "Species_Lari_Lar_v1_Gauss250",
  #"Species_Lari_Lar_v1_Gauss750",
  #"Species_Pice_Gla_v1",
  "Species_Pice_Gla_v1_Gauss250",
  #"Species_Pice_Gla_v1_Gauss750",
  #"Species_Pice_Mar_v1",
  "Species_Pice_Mar_v1_Gauss250",
  #"Species_Pice_Mar_v1_Gauss750",
  #"Species_Pinu_Ban_v1",
  "Species_Pinu_Ban_v1_Gauss250",
  #"Species_Pinu_Ban_v1_Gauss750",
  #"Species_Pinu_Con_v1",
  "Species_Pinu_Con_v1_Gauss250",
  #"Species_Pinu_Con_v1_Gauss750",
  #"Species_Popu_Bal_v1",
  "Species_Popu_Bal_v1_Gauss250",
  #"Species_Popu_Bal_v1_Gauss750",
  #"Species_Popu_Tre_v1",
  "Species_Popu_Tre_v1_Gauss250",
  #"Species_Popu_Tre_v1_Gauss750",
  #"Species_Pseu_Men_v1",
  "Species_Pseu_Men_v1_Gauss250",
  #"Species_Pseu_Men_v1_Gauss750",
  #"Species_Thuj_Pli_v1",
  "Species_Thuj_Pli_v1_Gauss250",
  #"Species_Thuj_Pli_v1_Gauss750",
  #"Species_Tsug_Het_v1",
  "Species_Tsug_Het_v1_Gauss250",
  #"Species_Tsug_Het_v1_Gauss750",
  #"Structure_Biomass_TotalLiveAboveGround_v1",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  #"Structure_Stand_Age_v1",
  "Structure_Stand_Age_v1_Gauss250",
  #"Structure_Stand_Age_v1_Gauss750",
  #"cti250",
  "cti250_Gauss250",
  #"cti250_Gauss750",
  "CultivationCrop",
  "IndustrialSiteRural",
  "MineSite",
  "Pipeline",
  "RoadHardSurface",
  "SeismicLineNarrow",
  "SeismicLineWide",
  "TransmissionLine",
  #"watAB",
  "watAB_Gauss250"
  #,
  #"watAB_Gauss750",
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind4 <- which(names(pred_abs_2011) %in% pred.variables4)

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr4 <- spatialAutoRange(pred_abs_2011[[ind4]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_GF250 <- spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr4$range,
  selection = "random",
  maskBySpecies = FALSE
)
end_time <- Sys.time()
end_time - start_time


start_time<-Sys.time()
brt4<- brt_blocks(data=datcombo,pred.variables = pred.variables4, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/GFsigma250m/", blocks=sp_block_GF250, save.points.shp = F, lr=0.01)  
end_time<-Sys.time()
end_time-start_time

save.image("D:/CHID regional Alberta BRT/BRT_outputs/outputs.RData")

# Only Gaussian filter with sigma = 750m ####
pred.variables5<-c( 
  #"Species_Abie_Bal_v1",
  #"Species_Abie_Bal_v1_Gauss250",
  "Species_Abie_Bal_v1_Gauss750",
  #"Species_Betu_Pap_v1",
  #"Species_Betu_Pap_v1_Gauss250",
  "Species_Betu_Pap_v1_Gauss750",
  #"Species_Lari_Lar_v1",
  #"Species_Lari_Lar_v1_Gauss250",
  "Species_Lari_Lar_v1_Gauss750",
  #"Species_Pice_Gla_v1",
  #"Species_Pice_Gla_v1_Gauss250",
  "Species_Pice_Gla_v1_Gauss750",
  #"Species_Pice_Mar_v1",
  #"Species_Pice_Mar_v1_Gauss250",
  "Species_Pice_Mar_v1_Gauss750",
  #"Species_Pinu_Ban_v1",
  #"Species_Pinu_Ban_v1_Gauss250",
  "Species_Pinu_Ban_v1_Gauss750",
  #"Species_Pinu_Con_v1",
  #"Species_Pinu_Con_v1_Gauss250",
  "Species_Pinu_Con_v1_Gauss750",
  #"Species_Popu_Bal_v1",
  #"Species_Popu_Bal_v1_Gauss250",
  "Species_Popu_Bal_v1_Gauss750",
  #"Species_Popu_Tre_v1",
  #"Species_Popu_Tre_v1_Gauss250",
  "Species_Popu_Tre_v1_Gauss750",
  #"Species_Pseu_Men_v1",
  #"Species_Pseu_Men_v1_Gauss250",
  "Species_Pseu_Men_v1_Gauss750",
  #"Species_Thuj_Pli_v1",
  #"Species_Thuj_Pli_v1_Gauss250",
  "Species_Thuj_Pli_v1_Gauss750",
  #"Species_Tsug_Het_v1",
  #"Species_Tsug_Het_v1_Gauss250",
  "Species_Tsug_Het_v1_Gauss750",
  #"Structure_Biomass_TotalLiveAboveGround_v1",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  #"Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
  "Structure_Stand_Age_v1_Gauss750",
  #"cti250",
  #"cti250_Gauss250",
  "cti250_Gauss750",
  "CultivationCrop",
  "IndustrialSiteRural",
  "MineSite",
  "Pipeline",
  "RoadHardSurface",
  "SeismicLineNarrow",
  "SeismicLineWide",
  "TransmissionLine",
  #"watAB",
  #"watAB_Gauss250",
  "watAB_Gauss750"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind5 <- which(names(pred_abs_2011) %in% pred.variables5)

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr5 <- spatialAutoRange(pred_abs_2011[[ind5]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_GF750 <- spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr5$range,
  selection = "random",
  maskBySpecies = FALSE
)
end_time <- Sys.time()
end_time - start_time


start_time<-Sys.time()
brt5<- brt_blocks(data=datcombo,pred.variables = pred.variables5, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/GFsigma750m/", blocks=sp_block_GF750, save.points.shp = F, lr=0.01)  
end_time<-Sys.time()
end_time-start_time




save.image("D:/CHID regional Alberta BRT/BRT_outputs/outputs.RData")


load("D:/CHID regional Alberta BRT/BRT_outputs/outputs.RData")

rast <-
  predict(pred_abs_2011,
          brt,
          type = "response",
          n.trees = brt$n.trees)



# Comparison of the different models ####
library(ggplot2)


# Model deviance
brt1$cv.statistics$deviance.mean
df<- data.frame(grp=c("Full model", "Selected scales", "cell level", "250m Gaussian smooth", "750m Gaussian smooth"),
                deviance=c(brt1$cv.statistics$deviance.mean,
                           brt2$cv.statistics$deviance.mean,
                           brt3$cv.statistics$deviance.mean,
                           brt4$cv.statistics$deviance.mean,
                           brt5$cv.statistics$deviance.mean),
                se=c(brt1$cv.statistics$deviance.se,
                     brt2$cv.statistics$deviance.se,
                     brt3$cv.statistics$deviance.se,
                     brt4$cv.statistics$deviance.se,
                     brt5$cv.statistics$deviance.se)
                )
df
k<-ggplot(df, aes(grp,deviance,ymin=deviance-se,ymax=deviance+se))
k+geom_pointrange()

# correlation
df2<- data.frame(grp=c("Full model", "Selected scales", "cell level", "250m Gaussian smooth", "750m Gaussian smooth"),
                 correlation=c(brt1$cv.statistics$correlation.mean,
                               brt2$cv.statistics$correlation.mean,
                               brt3$cv.statistics$correlation.mean,
                               brt4$cv.statistics$correlation.mean,
                               brt5$cv.statistics$correlation.mean),
                 se=c(brt1$cv.statistics$correlation.se,
                      brt2$cv.statistics$correlation.se,
                      brt3$cv.statistics$correlation.se,
                      brt4$cv.statistics$correlation.se,
                      brt5$cv.statistics$correlation.se)
                 )
df2
k2<-ggplot(df2, aes(grp,correlation,ymin=correlation-se,ymax=correlation+se))
k2+geom_pointrange()

# Calibration:  The first two statistics were the estimated intercepts and slopes of linear regression models of predictions against observations. The intercept measures the magnitude and direction of bias, with values close to 0 indicating low or no bias. The slope yields information about the consistency in the bias as a function of the mean, with a value of 1 indicating a consistent bias if the intercept is a nonzero value.
## intercept
df3<- data.frame(grp=c("Full model", "Selected scales", "cell level", "250m Gaussian smooth", "750m Gaussian smooth"),
                 calibration.intercept=c(brt1$cv.statistics$calibration.mean[1],
                                         brt2$cv.statistics$calibration.mean[1],
                                         brt3$cv.statistics$calibration.mean[1],
                                         brt4$cv.statistics$calibration.mean[1],
                                         brt5$cv.statistics$calibration.mean[1]),
                 se=c(brt1$cv.statistics$calibration.se[1],
                      brt2$cv.statistics$calibration.se[1],
                      brt3$cv.statistics$calibration.se[1],
                      brt4$cv.statistics$calibration.se[1],
                      brt5$cv.statistics$calibration.se[1])
                 )
df3
k3<-ggplot(df3, aes(grp,calibration.intercept,ymin=calibration.intercept-se,ymax=calibration.intercept+se))
k3+geom_pointrange()

## slope
df4<- data.frame(grp=c("Full model", "Selected scales", "cell level", "250m Gaussian smooth", "750m Gaussian smooth"),
                 calibration.slope=c(brt1$cv.statistics$calibration.mean[2],
                                         brt2$cv.statistics$calibration.mean[2],
                                         brt3$cv.statistics$calibration.mean[2],
                                         brt4$cv.statistics$calibration.mean[2],
                                         brt5$cv.statistics$calibration.mean[2]),
                 se=c(brt1$cv.statistics$calibration.se[2],
                      brt2$cv.statistics$calibration.se[2],
                      brt3$cv.statistics$calibration.se[2],
                      brt4$cv.statistics$calibration.se[2],
                      brt5$cv.statistics$calibration.se[2])
)
df4
k4<-ggplot(df4, aes(grp,calibration.slope,ymin=calibration.slope-se,ymax=calibration.slope+se))
k4+geom_pointrange()


# Plotting deviance for each brt model ####
dev_plot<-function(brt){
  m<-brt[[1]]
  y.bar <- min(m$cv.values) 
  y.min <- min(m$cv.values - m$cv.loss.ses)
  y.max <- max(m$cv.values + m$cv.loss.ses)
  
  plot(m$trees.fitted, m$cv.values, type = 'l', ylab = "Holdout deviance", xlab = "no. of trees", ylim = c(y.min,y.max))
  abline(h = y.bar, col = 3)
  
  lines(m$trees.fitted, m$cv.values + m$cv.loss.ses, lty=2)  
  lines(m$trees.fitted, m$cv.values - m$cv.loss.ses, lty=2)  
  
  target.trees <- m$trees.fitted[match(TRUE,m$cv.values == y.bar)]
  abline(v = target.trees, col=4)
}

dev_plot(brt2)



# Population size from density predictions ####
pred_abs_2011<-brick("D:/CHID regional Alberta BRT/prediction dataset/abs2011_250m.grd")

library(gbm)

rast <-  predict(pred_abs_2011,
                  brt2,
                  type = "response",
                  n.trees = brt2$n.trees)
plot(rast)

preds_brt2<-getValues(rast)[!is.na(getValues(rast))]

density_brt2_250m<-preds_brt2*6.25

abundance_brt2<-matrix(NA,length(density_brt2_250m),100)

for(i in 1:length(density_brt2_250m)){
  abundance_brt2[i,]<- rpois(100,density_brt2_250m[i])
}


mean(colSums(abundance_brt2))
sd(colSums(abundance_brt2))

sum(density_brt2_250m) 




####
# A function to obtain confidence intervals of model predictions ####
rm(list=setdiff(ls(),c("datcombo","sp_block_select","pred_abs_2011","blockCV_brt","brt2")))
gc()

boot_brt<-function(data,brtmodel,blocks,pred.data,nsamples=100,output.folder){
  #rast <-  predict(pred.data,
  #                 brtmodel,
  #                 type = "response",
  #                 n.trees = brtmodel$n.trees)
  #stack<-stack(rast)
  #names(stack)[1]<-"model estimate"
  
  z <- output.folder
  
  if (file.exists(z) == FALSE) {
    dir.create(z)
  }
  
  for(i in 1:nsamples){
    
    
    cat("loop",i,"\n") # this prints the loop number on console to track function progress
    
    
    sample<-sample(1:nrow(data),size=nrow(data),replace=T)
    data2<-data[sample,]
    x1<- try(brt <-
               gbm.step(
                 data = data2,
                 gbm.y = brtmodel$gbm.call$gbm.y,
                 gbm.x = brtmodel$gbm.call$gbm.x,
                 fold.vector = blocks$foldID[sample],
                 n.folds = brtmodel$gbm.call$cv.folds,
                 family = brtmodel$gbm.call$family,
                 tree.complexity = brtmodel$gbm.call$tree.complexity,
                 learning.rate = brtmodel$gbm.call$learning.rate,
                 bag.fraction = brtmodel$gbm.call$bag.fraction,
                 offset = brtmodel$gbm.call$offset[sample],
                 site.weights = data2$wt
               ))
    
    if(class(x1)=="NULL"){
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      x1<- try(brt <-
                 gbm.step(
                   data = data2,
                   gbm.y = brtmodel$gbm.call$gbm.y,
                   gbm.x = brtmodel$gbm.call$gbm.x,
                   fold.vector = blocks$foldID[sample],
                   n.folds = brtmodel$gbm.call$cv.folds,
                   family = brtmodel$gbm.call$family,
                   tree.complexity = brtmodel$gbm.call$tree.complexity,
                   learning.rate = brtmodel$gbm.call$learning.rate,
                   bag.fraction = brtmodel$gbm.call$bag.fraction,
                   offset = brtmodel$gbm.call$offset[sample],
                   site.weights = data2$wt
                 ))
    }
    
    if(class(x1)=="NULL"){
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      x1<- try(brt <-
                 gbm.step(
                   data = data2,
                   gbm.y = brtmodel$gbm.call$gbm.y,
                   gbm.x = brtmodel$gbm.call$gbm.x,
                   fold.vector = blocks$foldID[sample],
                   n.folds = brtmodel$gbm.call$cv.folds,
                   family = brtmodel$gbm.call$family,
                   tree.complexity = brtmodel$gbm.call$tree.complexity,
                   learning.rate = brtmodel$gbm.call$learning.rate,
                   bag.fraction = brtmodel$gbm.call$bag.fraction,
                   offset = brtmodel$gbm.call$offset[sample],
                   site.weights = data2$wt
                 ))
    }
    
    if(brt$n.trees>9999){
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      x1<- try(brt <-
                 gbm.step(
                   data = data2,
                   gbm.y = brtmodel$gbm.call$gbm.y,
                   gbm.x = brtmodel$gbm.call$gbm.x,
                   fold.vector = blocks$foldID[sample],
                   n.folds = brtmodel$gbm.call$cv.folds,
                   family = brtmodel$gbm.call$family,
                   tree.complexity = brtmodel$gbm.call$tree.complexity,
                   learning.rate = brtmodel$gbm.call$learning.rate,
                   bag.fraction = brtmodel$gbm.call$bag.fraction,
                   offset = brtmodel$gbm.call$offset[sample],
                   site.weights = data2$wt
                 ))
    }
    
    
    rast <- predict(pred.data,
                    brt,
                    type = "response",
                    n.trees = brt$n.trees)
    
    
    
    #stack <- addLayer(stack, rast)
    #names(stack)[i+1]<-paste0("sample",i) 
    
    writeRaster(rast, filename=paste0(z, "sample",i,".grd"), format="raster",overwrite=TRUE)
    rm(rast)
    gc()
    
  }
  
  #fun0.05 <- function(x) {quantile(x, probs = 0.05, na.rm = TRUE)}
  #lower<- calc(stack[[-1]],fun0.05)
  #fun0.95 <- function(x) {quantile(x, probs = 0.95, na.rm = TRUE)}
  #upper<- calc(stack[[-1]],fun0.95)
  
  #writeRaster(lower, filename=paste0(z, " confint_lower.tif"), format="GTiff",overwrite=TRUE)
  #writeRaster(upper, filename=paste0(z, " confint_upper.tif"), format="GTiff",overwrite=TRUE)
  
  #return(stack)
}

start_time <- Sys.time()
confintBRT2<-boot_brt(datcombo,brt2,sp_block_select,pred_abs_2011,nsamples=50,output.folder = "D://CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/run2/")
end_time <- Sys.time()
end_time - start_time

save.image("D:/CHID regional Alberta BRT/BRT_outputs/outputs2.RData")

# first run of the function made it to 66 samples, then crashed. these are in the object below:


start_time <- Sys.time()
confintBRT3<-boot_brt(datcombo,brt2,sp_block_select,pred_abs_2011,nsamples=134,output.folder = "D://CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/run3/")
end_time <- Sys.time()
end_time - start_time


# put together all bootstrap samples in one stack:

confintBRT<-brick("D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/samplesA.grd")

run2 <- list.files("D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/run2",pattern="gri$")
setwd("D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/run2/")
for (i in 1:length(run2)) { confintBRT <- addLayer(confintBRT, raster(run2[i]))}

run3 <- list.files("D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/run3",pattern="gri$")
setwd("D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/run3/")
for (i in 1:length(run3)) { confintBRT <- addLayer(confintBRT, raster(run3[i]))}

names(confintBRT)[68:251] <- paste0("sample",67:250)

writeRaster(confintBRT, filename="D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/all_samples.grd", format="raster",overwrite=TRUE)

# Calculate confidence intervals from bootstrap samples
fun0.05 <- function(x) {quantile(x, probs = 0.05, na.rm = TRUE)}
lower<- calc(confintBRT[[-1]],fun0.05)

fun0.95 <- function(x) {quantile(x, probs = 0.95, na.rm = TRUE)}
start_time <- Sys.time()
upper<- calc(confintBRT[[-1]],fun0.95)
end_time <- Sys.time()
end_time - start_time

writeRaster(lower, filename="D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/confint_lower.tif", format="GTiff",overwrite=TRUE)
writeRaster(upper, filename="D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/confint/confint_upper.tif", format="GTiff",overwrite=TRUE)


# Population size from density predictions ####
#pred_abs_2011<- brick("D:/CHID regional Alberta BRT/prediction dataset/abs2011_250m.grd")

library(gbm)

#rast <-  predict(pred_abs_2011,
#                 brt2,
#                 type = "response",
#                 n.trees = brt2$n.trees)
#plot(rast)


preds_brt2<-getValues(confintBRT[[1]])[!is.na(getValues(confintBRT[[1]]))]

density_brt2_250m<-preds_brt2*6.25

abundance_brt2<-matrix(NA,length(density_brt2_250m),100)

for(i in 1:length(density_brt2_250m)){
  abundance_brt2[i,]<- rpois(100,density_brt2_250m[i])
}


mean(colSums(abundance_brt2))
mean(colSums(abundance_brt2))*2
sd(colSums(abundance_brt2))

sum(density_brt2_250m) 

# 95% CI pop size
upper_brt2<-getValues(upper)[!is.na(getValues(upper))]
upper_brt2_250m<-upper_brt2*6.25
abundance_upper_brt2<-matrix(NA,length(upper_brt2_250m),100)
for(i in 1:length(upper_brt2_250m)){
  abundance_upper_brt2[i,]<- rpois(100,upper_brt2_250m[i])
}
mean(colSums(abundance_upper_brt2))
mean(colSums(abundance_upper_brt2))*2

# 5% CI pop size
lower_brt2<-getValues(lower)[!is.na(getValues(lower))]
lower_brt2_250m<-lower_brt2*6.25
abundance_lower_brt2<-matrix(NA,length(lower_brt2_250m),100)
for(i in 1:length(lower_brt2_250m)){
  abundance_lower_brt2[i,]<- rpois(100,lower_brt2_250m[i])
}
mean(colSums(abundance_lower_brt2))
mean(colSums(abundance_lower_brt2))*2


# Assessing sampling representativeness ### OUTDATED!!! ####

brt2$var.names
mess1<-mess(pred_abs_2011[[brt2$var.names]], extract(pred_abs_2011[[brt2$var.names]],cbind(datcombo$X,datcombo$Y)))
plot(mess1)
writeRaster(mess1, filename="D:/CHID regional Alberta BRT/BRT_outputs/selected_scales/mess.tif", format="GTiff",overwrite=TRUE)

gc()

# interrogate and find interactions for 'selected scales' model ####
find.int <- gbm.interactions(brt)
find.int$interactions
find.int$rank.list

gbm.perspec(brt, 2, 9, theta=310)

gbm.perspec(brt, 13, 9, theta=310)

gbm.perspec(brt, 13, 2, theta=310)

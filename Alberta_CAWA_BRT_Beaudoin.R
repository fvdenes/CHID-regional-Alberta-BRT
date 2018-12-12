library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(blockCV)

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
    try(brt1 <-
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
  
  # Define/create folders for storing outputs
  if (class(x1) != "try-error") {
    z <- output.folder
    
    if (file.exists(z) == FALSE) {
      dir.create(z)
    }
    
    if (is.null(blocks)){
      save(brt1, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    }
    
    else {
      save(blocks, file = paste(z, speclist[j], "blocks", sep = ""))
      save(brt1, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    }  
    
    
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
    
    data_sp <-SpatialPointsDataFrame(coords = data[, 33:34], data = data, proj4string = LCC)
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
    
    if(keep.out==T) {return(brt1)}
  }
}

#load data and prepare objects ####
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



rm(list=setdiff(ls(),c("datcombo","datcombo_sp","pred_abs_2011","w","LCC","speclist","randomCV_brt","j","blockCV_brt")))
gc()


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
  "AHM",
  "DD18",
  "MAT",
  "MAP",
  "FFP",
  "MWMT",
  "MCMT"
)

# convert data into SpatialPointDataFrame class
datcombo_sp <-SpatialPointsDataFrame(coords = datcombo[, 33:34], data = datcombo, proj4string = LCC)

# create a column converting abundance to occupancy
datcombo_sp$OCCU <- 0 
datcombo_sp$OCCU[which(datcombo_sp$ABUND > 0)] <- 1

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind <- which(names(pred_abs_2011) %in% pred.variables[-c(54:60)])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr1 <- spatialAutoRange(pred_abs_2011[[ind]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(seed)
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
brt1<- brt_blocks(data=datcombo,pred.variables = pred.variables, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/full_model/", blocks=sp_block_full, save.points.shp = TRUE)  
end_time<-Sys.time()
end_time-start_time

# Simplified model (selecting most influential scales) ####
pred.variables2<-c( 
  #"Species_Abie_Bal_v1",
  #"Species_Abie_Bal_v1_Gauss250",
  "Species_Abie_Bal_v1_Gauss750",
  #"Species_Betu_Pap_v1",
  #"Species_Betu_Pap_v1_Gauss250",
  "Species_Betu_Pap_v1_Gauss750",
  #"Species_Lari_Lar_v1",
  "Species_Lari_Lar_v1_Gauss250",
  #"Species_Lari_Lar_v1_Gauss750",
  #"Species_Pice_Gla_v1",
  #"Species_Pice_Gla_v1_Gauss250",
  "Species_Pice_Gla_v1_Gauss750",
  "Species_Pice_Mar_v1",
  #"Species_Pice_Mar_v1_Gauss250",
  #"Species_Pice_Mar_v1_Gauss750",
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
  "Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
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
  "AHM",
  "DD18",
  "MAT",
  "MAP",
  "FFP",
  "MWMT",
  "MCMT"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind2 <- which(names(pred_abs_2011) %in% pred.variables2[-c(24:30)])

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
brt2<- brt_blocks(data=datcombo,pred.variables = pred.variables2, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/selected_scales/", blocks=sp_block_select, save.points.shp = F)  
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
  "AHM",
  "DD18",
  "MAT",
  "MAP",
  "FFP",
  "MWMT",
  "MCMT"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind3 <- which(names(pred_abs_2011) %in% pred.variables3[-c(24:30)])

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
brt3<- brt_blocks(data=datcombo,pred.variables = pred.variables3, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/cell_level/", blocks=sp_block_cell, save.points.shp = F)  
end_time<-Sys.time()
end_time-start_time
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
  "AHM",
  "DD18",
  "MAT",
  "MAP",
  "FFP",
  "MWMT",
  "MCMT"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind4 <- which(names(pred_abs_2011) %in% pred.variables4[-c(24:30)])

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
brt4<- brt_blocks(data=datcombo,pred.variables = pred.variables4, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/GFsigma250m/", blocks=sp_block_GF250, save.points.shp = F)  
end_time<-Sys.time()
end_time-start_time

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
  "AHM",
  "DD18",
  "MAT",
  "MAP",
  "FFP",
  "MWMT",
  "MCMT"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind5 <- which(names(pred_abs_2011) %in% pred.variables5[-c(XX:YY)])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr5 <- spatialAutoRange(pred_abs_2011[[ind5]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(seed)
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
brt5<- brt_blocks(data=datcombo,pred.variables = pred.variables5, output.folder = "D://CHID regional Alberta BRT/BRT_outputs/GFsigma750m/", blocks=sp_block_GF750, save.points.shp = F)  
end_time<-Sys.time()
end_time-start_time

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


# Assessing sampling representativeness ### OUTDATED!!! ####
load("D:/CHID regional Alberta BRT/BRT_outputs/model1/CAWAbrtAB.R") # load example brt
brt1$var.names
mess1<-mess(abs2011_1km[[brt1$var.names]], extract(abs2011_1km[[brt1$var.names]],cbind(datcombo$X,datcombo$Y)))

writeRaster(mess1, filename="mess_brt_preds.tif", format="GTiff",overwrite=TRUE)



# mean and coefficient of variation among the 20 brt models ### OUTDATED!!! ####
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

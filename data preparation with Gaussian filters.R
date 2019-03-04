library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

# Preparing data

# Load avian datasets, AB region raster, project and extract data for AB, 
load("D:/Avian data processed/data_package_2016-04-18.Rdata")	
load("D:/Avian data processed/offsets-v3_2016-04-18.Rdata")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSLCC <- as.data.frame(spTransform(SS, LCC))
ABSS <- SSLCC[SSLCC$JURS=="AB",]
rm(SSLCC)
gc()

offl <- data.table(melt(OFF))
names(offl) <- c("PKEY","SPECIES","logoffset")
offl$SPECIES <- as.character(offl$SPECIES)
offl$PKEY <-as.character(offl$PKEY)
rm(OFF)
gc()

eco <- raster("D:/CHID regional Alberta BRT/spatial/albertaeco1.tif") 
alberta <- raster("D:/CHID regional Alberta BRT/spatial/AlbertaLCC.tif")
alberta<-projectRaster(alberta,crs=LCC)
plot(alberta)

# Load Beaudoin layers (2011 and 2001), crop and mask for AB, save as rasters
b2011 <- list.files("D:/Beaudoin/2011/Processed/sppBiomass_Canada",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
abs2011 <- crop(bs2011,alberta)
abs2011 <- mask(abs2011,alberta)
writeRaster(abs2011, filename="D:/Beaudoin/2011/Processed/AB/abs2011_250m.grd", format="raster",overwrite=TRUE)
abs2011<-stack("D:/Beaudoin/2011/Processed/AB/abs2011_250m.grd")

b2001 <- list.files("D:/Beaudoin/2001/Processed/sppBiomass_Canada",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
abs2001 <- crop(bs2001,alberta)
abs2001<-mask(abs2001,alberta)
writeRaster(abs2001, filename="D:/Beaudoin/2001/Processed/AB/abs2001_250m.grd", format="raster",overwrite=TRUE)
abs2001<-stack("D:/Beaudoin/2001/Processed/AB/abs2001_250m.grd")

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters

#Beaudoin 2011
## sigma = 250m
fw250<-focalWeight(x=abs2011,d=250,type="Gauss")
abs2011_Gauss250<-stack(focal(abs2011[[1]],w=fw250,na.rm=TRUE))
names(abs2011_Gauss250)<-names(abs2011)[[1]]
for(i in 2:nlayers(abs2011)){
 abs2011_Gauss250<-addLayer(abs2011_Gauss250,focal(abs2011[[i]],w=fw250,na.rm=TRUE))
 names(abs2011_Gauss250)[i]<-names(abs2011)[[i]]
}
abs2011_Gauss250<-brick(abs2011_Gauss250)
writeRaster(abs2011_Gauss250, filename="D:/Beaudoin/2011/Processed/AB/abs2011_250_Gauss250m.grd", format="raster",overwrite=TRUE)
#abs2011_Gauss250<-brick("D:/Beaudoin/2011/Processed/AB/abs2011_250_Gauss250m.grd")

## sigma = 750m
fw750<-focalWeight(x=abs2011,d=750,type="Gauss")
abs2011_Gauss750<-brick(focal(abs2011[[1]],w=fw750,na.rm=TRUE))
names(abs2011_Gauss750)<-names(abs2011)[[1]]
for(i in 2:nlayers(abs2011)){
 abs2011_Gauss750<-addLayer(abs2011_Gauss750,focal(abs2011[[i]],w=fw750,na.rm=TRUE))
 names(abs2011_Gauss750)[i]<-names(abs2011)[[i]]
}
abs2011_Gauss750<-brick(abs2011_Gauss750)
writeRaster(abs2011_Gauss750, filename="D:/Beaudoin/2011/Processed/AB/abs2011_250_Gauss750m.grd", format="raster",overwrite=TRUE)
#abs2011_Gauss750<-brick("D:/Beaudoin/2011/Processed/AB/abs2011_250_Gauss750m.grd")


#Beaudoin 2001
## sigma = 250m
abs2001_Gauss250<-stack(focal(abs2001[[1]],w=fw250,na.rm=TRUE))
names(abs2001_Gauss250)<-names(abs2001)[[1]]
for(i in 2:nlayers(abs2001)){
 abs2001_Gauss250<-addLayer(abs2001_Gauss250,focal(abs2001[[i]],w=fw250,na.rm=TRUE))
 names(abs2001_Gauss250)[i]<-names(abs2001)[[i]]
}
abs2001_Gauss250<-brick(abs2001_Gauss250)
writeRaster(abs2001_Gauss250, filename="D:/Beaudoin/2001/Processed/AB/abs2001_250_Gauss250m.grd", format="raster",overwrite=TRUE)
#abs2001_Gauss250<-brick("D:/Beaudoin/2001/Processed/AB/abs2001_250_Gauss250m.grd")


# ## sigma = 750m
abs2001_Gauss750<-brick(focal(abs2001[[1]],w=fw750,na.rm=TRUE))
names(abs2001_Gauss750)<-names(abs2001)[[1]]
for(i in 2:nlayers(abs2001)){
 abs2001_Gauss750<-addLayer(abs2001_Gauss750,focal(abs2001[[i]],w=fw750,na.rm=TRUE))
 names(abs2001_Gauss750)[i]<-names(abs2001)[[i]]
}
abs2001_Gauss750<-brick(abs2001_Gauss750)
writeRaster(abs2001_Gauss750, filename="D:/Beaudoin/2001/Processed/AB/abs2001_250_Gauss750m.grd", format="raster",overwrite=TRUE)
#abs2001_Gauss750<-brick("D:/Beaudoin/2001/Processed/AB/abs2001_250_Gauss750m.grd")

# cti data
cti100<-raster("D:/ABTerrainNielsen/100-m/cti.asc")
cti100<-projectRaster(cti100,crs=LCC)
cti250<-resample(cti100,alberta)
writeRaster(cti250, filename="D:/ABTerrainNielsen/250-m/cti.asc", format="ascii",overwrite=TRUE)
#cti250<-raster("D:/ABTerrainNielsen/250-m/cti.asc")

cti250_Gauss250<-focal(cti250,w=fw250,na.rm=TRUE)
cti250_Gauss750<-focal(cti250,w=fw750,na.rm=TRUE)

writeRaster(cti250_Gauss250, filename="D:/ABTerrainNielsen/250-m/Processed/cti250_Gauss250m.asc", format="ascii",overwrite=TRUE)
writeRaster(cti250_Gauss750, filename="D:/ABTerrainNielsen/250-m/Processed/cti250_Gauss750m.asc", format="ascii",overwrite=TRUE)

#cti250_Gauss250<-raster("D:/ABTerrainNielsen/250-m/Processed/cti250_Gauss250m.asc")
#cti250_Gauss750<-raster("D:/ABTerrainNielsen/250-m/Processed/cti250_Gauss750m.asc")

rm(cti100)
gc()

# Human Footprint data- upload and resample to 250m resolution to match other layers, and attach to abs2011
HF<-list.files("D:/ABMIstuff/",pattern = "tif$")
setwd("D:/ABMIstuff")
for(i in 1:length(HF)){ 
  HFlayer<-raster(HF[i])
  HFlayer<-resample(HFlayer,alberta)
  abs2011 <- addLayer(abs2011, HFlayer)
  names(abs2011)[nlayers(abs2011)] <- gsub("_AB_1km","",names(HFlayer))
}


# water and "lake edge density"
wat <- raster("D:/wat2011_lcc1/wat2011_lcc1.tif")
watAB<- crop(wat,alberta)
watAB<- mask(watAB,abs2011[[1]])

watAB_Gauss250<-focal(watAB,w=fw250,na.rm=TRUE)
watAB_Gauss750<-focal(watAB,w=fw750,na.rm=TRUE)

# climate data- upload and resample to 250m resolution to match other layers, and attach to abs2011
#climateAW2010 <- list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="asc$")
#setwd("D:/ClimateAdaptWest/baseline19812010/")
#clim2010 <- stack(raster(climateAW2010[1]))
#for (i in 2:length(climateAW2010)) { clim2010 <- addLayer(clim2010, raster(climateAW2010[i]))}
#proj4string(clim2010)<-LCC
#aclim2010 <- crop(clim2010,abs2011)
#aclim2010<-resample(clim2010,abs2011)
#aclim2010<-mask(aclim2010,abs2011$LandCover_NonVeg_v1)
 
#writeRaster(aclim2010, filename="D:/ClimateAdaptWest/baseline19812010/Processed/AB/aclim2010.grd", format="raster",overwrite=TRUE)
#aclim2010<-stack("D:/ClimateAdaptWest/baseline19812010/Processed/AB/aclim2010.grd")

#for(i in 1:length(names(aclim2010))){ 
#  abs2011 <- addLayer(abs2011, aclim2010[[i]])
#  names(abs2011)[nlayers(abs2011)] <- names(aclim2010[[i]])
#}

# put together prediction rasterstack
## need to update names of layers with Gaussian filters first to differentiate them
for(i in 1:nlayers(abs2001_Gauss250)){
  names(abs2001_Gauss250)[i] <- paste(names(abs2001_Gauss250)[i],"_Gauss250",sep="")
  names(abs2001_Gauss750)[i] <- paste(names(abs2001_Gauss750)[i],"_Gauss750",sep="")
  names(abs2011_Gauss250)[i] <- paste(names(abs2011_Gauss250)[i],"_Gauss250",sep="")
  names(abs2011_Gauss750)[i] <- paste(names(abs2011_Gauss750)[i],"_Gauss750",sep="")
}

names(cti250) <- "cti250"
names(cti250_Gauss250) <- "cti250_Gauss250"
names(cti250_Gauss750) <- "cti250_Gauss750"

names(watAB) <- "watAB"
names(watAB_Gauss250) <- "watAB_Gauss250"
names(watAB_Gauss750) <- "watAB_Gauss750"

pred_abs_2011<-stack(abs2011,abs2011_Gauss250,abs2011_Gauss750,cti250,cti250_Gauss250,cti250_Gauss750,watAB, watAB_Gauss250, watAB_Gauss750,quick=TRUE)
writeRaster(pred_abs_2011, filename="D:/CHID regional Alberta BRT/prediction dataset/abs2011_250m.grd", format="raster",overwrite=TRUE)


# Extracting data from rasters for surveyed locations (ABSS)
## 2011
dat2011 <- cbind(ABSS, extract(abs2011,as.matrix(cbind(ABSS$X,ABSS$Y)))) # includes Beaudoin layers (250m, no Gaussian filter), HF 

for(i in 1:nlayers(abs2011_Gauss250)){
  dat2011 <- cbind(dat2011, extract(abs2011_Gauss250[[i]],as.matrix(cbind(ABSS$X,ABSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=250m
  names(dat2011)[ncol(dat2011)] <- names(abs2011_Gauss250)[i]
}

for(i in 1:nlayers(abs2011_Gauss750)){
  dat2011 <- cbind(dat2011, extract(abs2011_Gauss750[[i]],as.matrix(cbind(ABSS$X,ABSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=750m
  names(dat2011)[ncol(dat2011)] <- names(abs2011_Gauss750)[i]
}

dat2011<-cbind(dat2011,extract(cti250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes cti 250m resolution data
names(dat2011)[ncol(dat2011)] <- "cti250"

dat2011<-cbind(dat2011,extract(cti250_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes cti 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "cti250_Gauss250"

dat2011<-cbind(dat2011,extract(cti250_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes cti 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "cti250_Gauss750"

dat2011<-cbind(dat2011,extract(watAB,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data
names(dat2011)[ncol(dat2011)] <- "watAB"

dat2011<-cbind(dat2011,extract(watAB_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "watAB_Gauss250"

dat2011<-cbind(dat2011,extract(watAB_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "watAB_Gauss750"

dat2011 <-cbind(dat2011,extract(eco,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes ecoregions layer
names(dat2011)[ncol(dat2011)] <- "eco"

### set up weight matrix for SS, and calculate weight values for each row in dat2011
r2 <- abs2011[[1]]
samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1)
sampsum25 <- focal(samprast2011, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)

dat2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "sampsum25"
dat2011$wt <- 1/dat2011$sampsum25

dat2011$SS <- as.character(dat2011$SS)
dat2011$PCODE <- as.character(dat2011$PCODE)
dat2011<-dat2011[,-c(25:43)] # remove columns with climate data (from avian dataset, since using Adaptwest Climate data instead)

setwd("D:/CHID regional Alberta BRT/")
write.csv(dat2011,"ABdat2011.csv")

## 2001
dat2001 <- cbind(ABSS, extract(abs2001,as.matrix(cbind(ABSS$X,ABSS$Y)))) # includes Beaudoin layers (250m, no Gaussian filter), HF and climate

for(i in which(names(abs2011)=="CultivationAbandoned"):which(names(abs2011)=="TransmissionLine") ) { # copy human footprint  data from dat2011
  dat2001<-cbind(dat2001,extract(abs2011[[i]],as.matrix(cbind(ABSS$X,ABSS$Y))))
  names(dat2001)[ncol(dat2001)] <- names(abs2011)[[i]]
}

for(i in 1:nlayers(abs2001_Gauss250)){
  dat2001 <- cbind(dat2001, extract(abs2001_Gauss250[[i]],as.matrix(cbind(ABSS$X,ABSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=250m
  names(dat2001)[ncol(dat2001)] <- names(abs2001_Gauss250)[i]
}

for(i in 1:nlayers(abs2001_Gauss750)){
  dat2001 <- cbind(dat2001, extract(abs2001_Gauss750[[i]],as.matrix(cbind(ABSS$X,ABSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=750m
  names(dat2001)[ncol(dat2001)] <- names(abs2001_Gauss750)[i]
}

dat2001<-cbind(dat2001,extract(cti250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes cti 250m resolution data
names(dat2001)[ncol(dat2001)] <- "cti250"

dat2001<-cbind(dat2001,extract(cti250_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes cti 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "cti250_Gauss250"

dat2001<-cbind(dat2001,extract(cti250_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes cti 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "cti250_Gauss750"

dat2001<-cbind(dat2001,extract(watAB,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data
names(dat2001)[ncol(dat2001)] <- "watAB"

dat2001<-cbind(dat2001,extract(watAB_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "watAB_Gauss250"

dat2001<-cbind(dat2001,extract(watAB_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "watAB_Gauss750"

dat2001 <-cbind(dat2001,extract(eco,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes ecoregions layer
names(dat2001)[ncol(dat2001)] <- "eco"

### set up weight matrix for SS, and calculate weight values for each row in dat2001
samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25

dat2001$SS <- as.character(dat2001$SS)
dat2001$PCODE <- as.character(dat2001$PCODE)
dat2001<-dat2001[,-c(25:43)] # remove columns with climate data (from avian dataset, since using Adaptwest Climate data instead)
write.csv(dat2001,"ABdat2001.csv")


# Prepare point count data for each SS and aggregate for 2001 and 2011.


PC <- inner_join(PCTBL,PKEY[,1:8],by=c("PKEY","SS"))[,-1]
colnames(PC)[10]<-"PCODE"
PC <- inner_join(PC,SS@data[,c(2,5)],by="SS")
ABPC <- PC[PC$JURS=="AB",]

ABPC$SS <- as.character(ABPC$SS)
ABPC$PKEY <- as.character(ABPC$PKEY)
ABPC$PCODE <- as.character(ABPC$PCODE)
ABPC$SPECIES <- as.character(ABPC$SPECIES)
ABPC2001 <- ABPC[ABPC$YEAR < 2006,] #n=418337
ABPC2011 <- ABPC[ABPC$YEAR > 2005,] #n=400780
survey2001 <- aggregate(ABPC2001$ABUND, by=list("PKEY"=ABPC2001$PKEY,"SS"=ABPC2001$SS,"PCODE"=ABPC2001$PCODE), FUN=sum) #n=31838
survey2011 <- aggregate(ABPC2011$ABUND, by=list("PKEY"=ABPC2011$PKEY,"SS"=ABPC2011$SS,"PCODE"=ABPC2011$PCODE), FUN=sum) #n=26239

speclist<-levels(as.factor(offl$SPECIES))



rm(list=setdiff(ls(),c("pred_abs_2011","LCC","speclist","offl","ABPC2001","survey2001","dat2001","ABPC2011","survey2011","dat2011")))
gc()
save.image("D:/CHID regional Alberta BRT/AB_BRT_Rproject/data_pack.RData")

library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

# Preparing data


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
nalc <- raster("D:/CHID regional Alberta BRT/spatial/NA_LandCover_2005_LCC.img")
alberta <- raster("D:/CHID regional Alberta BRT/spatial/AlbertaLCC.tif")
alberta<-projectRaster(alberta,crs=LCC)
plot(alberta)

b2011 <- list.files("D:/Beaudoin/2011/",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))

abs2011 <- crop(bs2011,alberta)
abs2011 <- mask(abs2011,alberta)

writeRaster(abs2011, filename="abs2011_250m.grd", format="raster",overwrite=TRUE)

gc()
abs2011_1km <- aggregate(abs2011, fact=4, fun=mean)

ab1km<-aggregate(alberta, fact=4, fun=mean)
abs2011_1km<-mask(abs2011_1km,ab1km)

rm(ab1km)
gc()

r2 <- abs2011_1km[[1]]

ecor1km <- resample(eco, abs2011_1km)
abs2011_1km <- addLayer(abs2011_1km, ecor1km)
names(abs2011_1km)[nlayers(abs2011_1km)] <- "eco"

cti500<-raster("D:/ABTerrainNielsen/500-m/cti.asc")
cti500<-projectRaster(cti500,crs=LCC)
cti_1km<-aggregate(cti500,fact=4, fun=mean)
cti_1km<-resample(cti_1km,abs2011_1km)
gc()

abs2011_1km <- addLayer(abs2011_1km,cti_1km)
names(abs2011_1km)[nlayers(abs2011_1km)] <- "cti"


HF<-list.files("D:/ABMIstuff/",pattern = "tif$")
setwd("D:/ABMIstuff")
for(i in 1:length(HF)){ 
  HFlayer<-raster(HF[i])
  HFlayer<-resample(HFlayer,abs2011_1km)
  abs2011_1km <- addLayer(abs2011_1km, HFlayer)
  names(abs2011_1km)[nlayers(abs2011_1km)] <- gsub("_AB_1km","",names(HFlayer))
  }


climateAW2010 <- list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="asc$")
setwd("D:/ClimateAdaptWest/baseline19812010/")
clim2010 <- stack(raster(climateAW2010[1]))
for (i in 2:length(climateAW2010)) { clim2010 <- addLayer(clim2010, raster(climateAW2010[i]))}
proj4string(clim2010)<-LCC
aclim2010 <- crop(clim2010,abs2011_1km)
aclim2010<-resample(clim2010,abs2011_1km)
aclim2010<-mask(aclim2010,abs2011_1km$LandCover_NonVeg_v1)

for(i in 1:length(names(aclim2010))){ 
  abs2011_1km <- addLayer(abs2011_1km, aclim2010[[i]])
  names(abs2011_1km)[nlayers(abs2011_1km)] <- names(aclim2010[[i]])
}



b2001 <- list.files("D:/Beaudoin/2001/",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
abs2001 <- crop(bs2001,alberta)
abs2001<-mask(abs2001,alberta)

dat2011 <- cbind(ABSS, extract(abs2011,as.matrix(cbind(ABSS$X,ABSS$Y))))
dat2011 <-cbind(dat2011,extract(nalc,as.matrix(cbind(dat2011$X,dat2011$Y)))) 
names(dat2011)[ncol(dat2011)] <- "LCC"
dat2011 <-cbind(dat2011,extract(eco,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "eco"

gc()

####

cti100<-raster("D:/ABTerrainNielsen/100-m/cti.asc")
cti100<-projectRaster(cti100,crs=LCC)
dat2011<-cbind(dat2011,extract(cti100,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "cti"

for(i in which(names(abs2011_1km)=="CultivationAbandoned"):which(names(abs2011_1km)=="TD") ) {
  dat2011<-cbind(dat2011,extract(abs2011_1km[[i]],as.matrix(cbind(ABSS$X,ABSS$Y))))
  names(dat2011)[ncol(dat2011)] <- names(abs2011_1km)[[i]]
}



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

dat2001 <- cbind(ABSS, extract(abs2001,as.matrix(cbind(ABSS$X,ABSS$Y))))
dat2001 <-cbind(dat2001,extract(nalc,as.matrix(cbind(dat2001$X,dat2001$Y)))) 
names(dat2001)[ncol(dat2001)] <- "LCC"
dat2001 <-cbind(dat2001,extract(eco,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "eco"
dat2001<-cbind(dat2001,extract(cti100,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "cti"

for(i in which(names(abs2011_1km)=="CultivationAbandoned"):which(names(abs2011_1km)=="TD") ) {
  dat2001<-cbind(dat2001,extract(abs2011_1km[[i]],as.matrix(cbind(ABSS$X,ABSS$Y))))
  names(dat2001)[ncol(dat2001)] <- names(abs2011_1km)[[i]]
}



samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25

dat2001$SS <- as.character(dat2001$SS)
dat2001$PCODE <- as.character(dat2001$PCODE)
dat2001<-dat2001[,-c(25:43)] # remove columns with climate data (from avian dataset, since using Adaptwest Climate data instead)
write.csv(dat2001,"ABdat2001.csv")

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

writeRaster(abs2011_1km, filename="abs2011_1km.grd", format="raster",overwrite=TRUE)

save.image("D:/CHID regional Alberta BRT/data_pack.RData")

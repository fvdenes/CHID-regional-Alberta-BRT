library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

# Preparing data


load("M:/DataStuff/AvianData/Processed/data_package_2016-04-18.Rdata")	
load("M:/DataStuff/AvianData/Processed/offsets-v3_2016-04-18.Rdata")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSLCC <- as.data.frame(spTransform(SS, LCC))
ABSS <- SSLCC[SSLCC$JURS=="AB",]

offl <- data.table(melt(OFF))
names(offl) <- c("PKEY","SPECIES","logoffset")
offl$SPECIES <- as.character(offl$SPECIES)
offl$PKEY <-as.character(offl$PKEY)
rm(OFF)


eco <- raster("D:/CHID regional Alberta BRT/albertaeco1.tif") 
nalc <- raster("D:/CHID regional Alberta BRT/NA_LandCover_2005_LCC.img")
alberta <- raster("D:/CHID regional Alberta BRT/AlbertaLCC.tif") 

b2011 <- list.files("M:/DataStuff/SpatialData/Beaudoin/2011/",pattern="tif$")
setwd("M:/DataStuff/SpatialData/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
abs2011 <- crop(bs2011,alberta)
abs2011_1km <- aggregate(abs2011, fact=4, fun=mean)
r2 <- abs2011_1km[[1]]

ecor1km <- resample(eco, abs2011_1km)
abs2011_1km <- addLayer(abs2011_1km, ecor1km)
names(abs2011_1km)[nlayers(abs2011_1km)] <- "eco"
writeRaster(abs2011_1km,"AB2011rasters")

b2001 <- list.files("M:/DataStuff/SpatialData/Beaudoin/2001/",pattern="tif$")
setwd("M:/DataStuff/SpatialData/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
abs2001 <- crop(bs2001,alberta)

dat2011 <- cbind(ABSS, extract(abs2011,as.matrix(cbind(ABSS$X,ABSS$Y))))
dat2011 <-cbind(dat2011,extract(nalc,as.matrix(cbind(dat2011$X,dat2011$Y)))) 
names(dat2011)[ncol(dat2011)] <- "LCC"
dat2011 <-cbind(dat2011,extract(eco,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "eco"

samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1)
sampsum25 <- focal(samprast2011, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)

dat2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "sampsum25"
dat2011$wt <- 1/dat2011$sampsum25

dat2011$SS <- as.character(dat2011$SS)
dat2011$PCODE <- as.character(dat2011$PCODE)

setwd("D://CHID regional Alberta BRT/")
write.csv(dat2011,"ABdat2011.csv")

dat2001 <- cbind(ABSS, extract(abs2001,as.matrix(cbind(ABSS$X,ABSS$Y))))
dat2001 <-cbind(dat2001,extract(nalc,as.matrix(cbind(dat2001$X,dat2001$Y)))) 
names(dat2001)[ncol(dat2001)] <- "LCC"
dat2001 <-cbind(dat2001,extract(eco,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "eco"

samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25

dat2001$SS <- as.character(dat2001$SS)
dat2001$PCODE <- as.character(dat2001$PCODE)
write.csv(dat2001,"ABdat2001.csv")

PC <- inner_join(PCTBL,PKEY[,1:8],by=c("PKEY","SS","PCODE"))
PC <- inner_join(PC,SS@data[,c(2,5)],by="SS")
ABPC <- PC[PC$JURS=="AB",]

ABPC$SS <- as.character(ABPC$SS)
ABPC$PKEY <- as.character(ABPC$PKEY)
ABPC$PCODE <- as.character(ABPC$PCODE)
ABPC$SPECIES <- as.character(ABPC$SPECIES)
ABPC2001 <- ABPC[ABPC$YEAR < 2006,] #n=180998
ABPC2011 <- ABPC[ABPC$YEAR > 2005,] #n=215457
survey2001 <- aggregate(ABPC2001$ABUND, by=list("PKEY"=ABPC2001$PKEY,"SS"=ABPC2001$SS,"PCODE"=ABPC2001$PCODE), FUN=sum) #n=31838
survey2011 <- aggregate(ABPC2011$ABUND, by=list("PKEY"=ABPC2011$PKEY,"SS"=ABPC2011$SS,"PCODE"=ABPC2011$PCODE), FUN=sum) #n=26239

speclist<-levels(as.factor(offl$SPECIES))

save.image("D:/CHID regional Alberta BRT/data_pack.RData")
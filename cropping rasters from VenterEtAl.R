library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)




LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
alberta <- raster("D:/CHID regional Alberta BRT/CAWA_pred1km.tif") 
croplands<-raster("D:/CHID regional Alberta BRT/croplands2005_AB.tif")

croplands<-crop(croplands,alberta)
croplands<-mask(croplands,alberta)
croplands<-projectRaster(croplandsM,crs=LCC)
plot(croplands)
croplands<-crop(croplands,alberta)

writeRaster(croplands,"croplands2005AB",format="GTiff",overwrite=TRUE)


pasture<-raster("D:/CHID regional Alberta BRT/pasture2009AB.tif")
extent(pasture)

pasture<-mask(pasture,alberta)
writeRaster(pasture,"pasture2009AB",format="GTiff",overwrite=TRUE)

roads<-raster("D:/CHID regional Alberta BRT/roads_AB.tif")
plot(roads)
roads<-mask(roads,alberta)

writeRaster(roads,"roads_AB",format="GTiff",overwrite=TRUE)

library(raster)

b2001 <- list.files("D:/Beaudoin/2001/",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}

b2011 <- list.files("D:/Beaudoin/2011/",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}

biomass2001<-bs2001
biomass2011<-bs2011


start_time<-Sys.time()
for (i in 5:82){
  biomass2001[[i]]<-bs2001[[i]]/100*bs2001[[88]]
  biomass2011[[i]]<-bs2011[[i]]/100*bs2011[[88]]
}
end_time<-Sys.time()
end_time-start_time

start_time<-Sys.time()
for(i in 1:nlayers(biomass2001)){
  writeRaster(biomass2001[[i]], filename=paste("D:/Beaudoin/2001/Processed/sppBiomass_Canada_t_per_ha/",names(bs2001)[i],sep=""), format="GTiff",overwrite=TRUE)
  writeRaster(biomass2011[[i]], filename=paste("D:/Beaudoin/2011/Processed/sppBiomass_Canada_t_per_ha/",names(bs2011)[i],sep=""), format="GTiff",overwrite=TRUE)
}
end_time<-Sys.time()
end_time-start_time
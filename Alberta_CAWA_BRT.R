# R code to develop boosted regression tree models and generate current predictions
library(raster) #Reading, writing, manipulating, analyzing and modeling of gridded spatial data
library(dismo) #Species distribution modeling
library(gbm) #Generalized boosted regression models
library(sp)
library(rgdal)
library(viridis)


load("D:/CHID regional Alberta BRT/data_pack.RData")

head(ABDAT)
head(pred_data)

colnames(pred_data)[6]<-"HABTR"
pred_data$YR<-12


summary(ABDAT$CAWAoffset)
table(ABDAT$CAWAoffset)

which(colnames(ABDAT)=="HABTR")
which(colnames(ABDAT)=="HGT")
which(colnames(ABDAT)=="CTI")
which(colnames(ABDAT)=="MSP")
which(colnames(ABDAT)=="CMI")
which(colnames(ABDAT)=="TD")

which(colnames(ABDAT)=="CAWAoffset")

which(colnames(ABDAT)=="count")





set.seed(140918)
cawa.AB.tc3.lr0005<- gbm.step(data=ABDAT, gbm.x = c(11,17,29,31,33,34,35,36,37,41,43,44), gbm.y = 45, family = "poisson", offset=ABDAT$CAWAoffset, tree.complexity = 3,learning.rate = 0.001, bag.fraction = 0.5)
summary(cawa.AB.tc3.lr0005)
gbm.plot(cawa.AB.tc3.lr0005, write.title = F)
gbm.plot.fits(cawa.AB.tc3.lr0005)
find.int <- gbm.interactions(cawa.AB.tc3.lr0005)
find.int$rank.list
find.int$interactions

gbm.perspec(cawa.AB.tc3.lr0005, 1, 3)
gbm.perspec(cawa.AB.tc3.lr0005, 2, 8)


pred_data$predictions<- predict.gbm(cawa.AB.tc3.lr0005, pred_data, type="response", n.trees=cawa.AB.tc3.lr0005$n.trees)

sp_preds<-SpatialPixelsDataFrame(points=cbind(pred_data$POINT_X,pred_data$POINT_Y),proj4string =  CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"),data=pred_data, tolerance=0.9)

save.image("D:/CHID regional Alberta BRT/CAWA_Alberta_BRT.RData")

plot(sp_preds["predictions"], col=rev(magma(n=50,  end=1)),scale.shrink=2, zlim=c(0,0.2))

writeRaster(raster(sp_preds["predictions"]), filename="brt_CAWA_AB.asc",prj=TRUE, format="ascii",overwrite=TRUE)


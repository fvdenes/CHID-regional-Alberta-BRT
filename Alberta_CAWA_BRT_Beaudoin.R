library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

load("D:/CHID regional Alberta BRT/AB_BRT_Rproject/data_pack.RData")

# list of tree spp for AB:


# define list 
w <-"D://CHID regional Alberta BRT/BRT_outputs/"

#for (j in 1:length(speclist)) { # to run for all spp. Make sure to uncomment the bracket at the end
j<-which(speclist=="CAWA") # to run only for CAWA
  
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
  
  rm(list=setdiff(ls(),c("datcombo","pred_abs_2011","w","LCC","speclist","randomCV_brt","j")))
  gc()
  #Setting up list of variables for model
randomCV_brt(datcombo,nmodels=19,holdout = 0.3)

randomCV_brt(datcombo,nmodels=4,holdout = 0.3)

  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = "ABUND", 
                             gbm.x = c("LandCover_NonVeg_v1",
                                       "Species_Abie_Bal_v1",
                                       "Species_Abie_Las_v1",
                                       "Species_Betu_Pap_v1",
                                       "Species_Lari_Lar_v1",
                                       "Species_Pice_Gla_v1",
                                       "Species_Pice_Mar_v1",
                                       "Species_Pinu_Ban_v1",
                                       "Species_Pinu_Con_v1",
                                       "Species_Popu_Bal_v1",
                                       "Species_Popu_Tre_v1",
                                       "Species_Pseu_Men_v1",
                                       "Species_Thuj_Pli_v1",
                                       "Species_Tsug_Het_v1",
                                       "Structure_Biomass_TotalLiveAboveGround_v1",
                                       "Structure_Stand_Age_v1",
                                       "LandCover_NonVeg_v1_Gauss250",
                                       "Species_Abie_Bal_v1_Gauss250",
                                       "Species_Abie_Las_v1_Gauss250",
                                       "Species_Betu_Pap_v1_Gauss250",
                                       "Species_Lari_Lar_v1_Gauss250",
                                       "Species_Pice_Gla_v1_Gauss250",
                                       "Species_Pice_Mar_v1_Gauss250",
                                       "Species_Pinu_Ban_v1_Gauss250",
                                       "Species_Pinu_Con_v1_Gauss250",
                                       "Species_Popu_Bal_v1_Gauss250",
                                       "Species_Popu_Tre_v1_Gauss250",
                                       "Species_Pseu_Men_v1_Gauss250",
                                       "Species_Thuj_Pli_v1_Gauss250",
                                       "Species_Tsug_Het_v1_Gauss250",
                                       "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
                                       "Structure_Stand_Age_v1_Gauss250",
                                       "LandCover_NonVeg_v1_Gauss500",
                                       "Species_Abie_Bal_v1_Gauss500",
                                       "Species_Abie_Las_v1_Gauss500",
                                       "Species_Betu_Pap_v1_Gauss500",
                                       "Species_Lari_Lar_v1_Gauss500",
                                       "Species_Pice_Gla_v1_Gauss500",
                                       "Species_Pice_Mar_v1_Gauss500",
                                       "Species_Pinu_Ban_v1_Gauss500",
                                       "Species_Pinu_Con_v1_Gauss500",
                                       "Species_Popu_Bal_v1_Gauss500",
                                       "Species_Popu_Tre_v1_Gauss500",
                                       "Species_Pseu_Men_v1_Gauss500",
                                       "Species_Thuj_Pli_v1_Gauss500",
                                       "Species_Tsug_Het_v1_Gauss500",
                                       "Structure_Biomass_TotalLiveAboveGround_v1_Gauss500",
                                       "Structure_Stand_Age_v1_Gauss500",
                                       "cti250",
                                       "cti250_Gauss250",
                                       "cti250_Gauss500",
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
                                       "MAP",#in AB it might make more sense to use this instead of MSP because winters are so dry
                                       "FFP",
                                       "MWMT",
                                       "MCMT"), 
                             family = "poisson", tree.complexity = 3, learning.rate = 0.025, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "try-error") {
    save(brt1,file=paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","brtAB.R",sep=""))
    varimp <- as.data.frame(brt1$contributions)
    write.csv(varimp,file=paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","varimp.csv",sep=""))
    cvstats <- t(as.data.frame(brt1$cv.statistics))
    write.csv(cvstats,file=paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","cvstats.csv",sep=""))
    #find.int<-gbm.interactions(x1)
    #find.int$rank.list
    #find.int$interactions
    #gbm.perspec(x1,14,16)
    #pdf(paste(w,speclist[j],"_plot.pdf",sep=""))
    #gbm.plot(brt1,n.plots=32,smooth=TRUE, plot.layout = c(3,3),common.scale=T)
    #dev.off()
    pdf(paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","cawa_plot.var-scale.pdf",sep=""))
    gbm.plot(brt1,n.plots=66,smooth=TRUE, plot.layout = c(3,3),common.scale=F,write.title = F)
    dev.off()
    #pdf(paste(w,speclist[j],"_plot.fitted.pdf",sep=""))
    #gbm.plot.fits(brt1)
    #dev.off()
    rast <- predict(pred_abs_2011, brt1, type="response", n.trees=brt1$n.trees)
    writeRaster(rast, filename=paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","cawa_pred250m",sep=""), format="GTiff",overwrite=TRUE)
    png(paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","cawa_pred250m.png",sep=""))
    plot(rast, zlim=c(0,1))
    points(datcombo$X, datcombo$Y, cex=0.05)
    dev.off()
  }
  
  # reordering varimp according to variable name (without Gauss filer) and relative influence
  varimp$ordervar<-gsub("_Gauss500","",varimp$var)
  varimp$ordervar<-gsub("_Gauss250","",varimp$ordervar)
  varimp[order(varimp$ordervar,-varimp$rel.inf),]
#}

  ### "weeding out less influential covariates (among those with Gaussian filters and those with 0 influence)
  simp1<- gbm.simplify(brt1)
  
  
  
  
  
  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = "ABUND", 
                             gbm.x = c(#"LandCover_NonVeg_v1",
                                       #"Species_Abie_Bal_v1",
                                       #"Species_Abie_Las_v1",
                                       #"Species_Betu_Pap_v1",
                                       #"Species_Lari_Lar_v1",
                                       #"Species_Pice_Gla_v1",
                                       #"Species_Pice_Mar_v1",
                                       #"Species_Pinu_Ban_v1",
                                       #"Species_Pinu_Con_v1",
                                       #"Species_Popu_Bal_v1",
                                       #"Species_Popu_Tre_v1",
                                       #"Species_Pseu_Men_v1",
                                       #"Species_Thuj_Pli_v1",
                                       #"Species_Tsug_Het_v1",
                                       #"Structure_Biomass_TotalLiveAboveGround_v1",
                                       #"Structure_Stand_Age_v1",
                                       #"LandCover_NonVeg_v1_Gauss250",
                                       #"Species_Abie_Bal_v1_Gauss250",
                                       #"Species_Abie_Las_v1_Gauss250",
                                       #"Species_Betu_Pap_v1_Gauss250",
                                       "Species_Lari_Lar_v1_Gauss250",
                                       #"Species_Pice_Gla_v1_Gauss250",
                                       #"Species_Pice_Mar_v1_Gauss250",
                                       #"Species_Pinu_Ban_v1_Gauss250",
                                       #"Species_Pinu_Con_v1_Gauss250",
                                       #"Species_Popu_Bal_v1_Gauss250",
                                       "Species_Popu_Tre_v1_Gauss250",
                                       #"Species_Pseu_Men_v1_Gauss250",
                                       #"Species_Thuj_Pli_v1_Gauss250",
                                       #"Species_Tsug_Het_v1_Gauss250",
                                       "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
                                       "Structure_Stand_Age_v1_Gauss250",
                                       "LandCover_NonVeg_v1_Gauss500",
                                       "Species_Abie_Bal_v1_Gauss500",
                                       #"Species_Abie_Las_v1_Gauss500",
                                       "Species_Betu_Pap_v1_Gauss500",
                                       #"Species_Lari_Lar_v1_Gauss500",
                                       "Species_Pice_Gla_v1_Gauss500",
                                       "Species_Pice_Mar_v1_Gauss500",
                                       "Species_Pinu_Ban_v1_Gauss500",
                                       "Species_Pinu_Con_v1_Gauss500",
                                       "Species_Popu_Bal_v1_Gauss500",
                                       #"Species_Popu_Tre_v1_Gauss500",
                                       #"Species_Pseu_Men_v1_Gauss500",
                                       #"Species_Thuj_Pli_v1_Gauss500",
                                       #"Species_Tsug_Het_v1_Gauss500",
                                       #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss500",
                                       #"Structure_Stand_Age_v1_Gauss500",
                                       #"cti250",
                                       #"cti250_Gauss250",
                                       "cti250_Gauss500",
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
                                       "MAP",#in AB it might make more sense to use this instead of MSP because winters are so dry
                                       "FFP",
                                       "MWMT",
                                       "MCMT"), 
                             family = "poisson", tree.complexity = 3, learning.rate = 0.025, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "try-error") {
    save(brt1,file=paste(w,speclist[j],"brtAB.R",sep=""))
    varimp <- as.data.frame(brt1$contributions)
    write.csv(varimp,file=paste(speclist[j],"varimp.csv",sep=""))
    cvstats <- t(as.data.frame(brt1$cv.statistics))
    write.csv(cvstats,file=paste(w,speclist[j],"cvstats.csv",sep=""))
    #find.int<-gbm.interactions(x1)
    #find.int$rank.list
    #find.int$interactions
    #gbm.perspec(x1,14,16)
    pdf(paste(w,speclist[j],"_plot.pdf",sep=""))
    gbm.plot(brt1,n.plots=28,smooth=TRUE, plot.layout = c(3,3),common.scale=T)
    dev.off()
    pdf(paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","cawa_simp_plot.var-scale.pdf",sep=""))
    gbm.plot(brt1,n.plots=28,smooth=TRUE, plot.layout = c(3,3),common.scale=F,write.title = F)
    dev.off()
    pdf(paste(w,speclist[j],"_plot.fitted.pdf",sep=""))
    gbm.plot.fits(brt1)
    dev.off()
    rast2 <- predict(pred_abs_2011, brt1, type="response", n.trees=brt1$n.trees)
    writeRaster(rast2, filename=paste("D:/CHID regional Alberta BRT/BRT_outputs/fullmodel_Gaussianfilters/","cawa_pred250m_simplified",sep=""), format="GTiff",overwrite=TRUE)
    png(paste(speclist[j],"_pred1km.png",sep=""))
    plot(rast, zlim=c(0,1))
    points(datcombo$X, datcombo$Y, cex=0.05)
    dev.off()
  }
  
    
 surveypoints<-SpatialPointsDataFrame(coords=cbind(datcombo$X,datcombo$Y),data=datcombo,proj4string = LCC)
 writeOGR(surveypoints,paste(w,"surveypoints.shp",sep=""),"surveypoints", driver="ESRI Shapefile")  

 
 rm(surveypoints)
 gc()
 
 # Assessing sampling representativeness
 load("D:/CHID regional Alberta BRT/BRT_outputs/model1/CAWAbrtAB.R") # load example brt
 brt1$var.names
  mess1<-mess(abs2011_1km[[brt1$var.names]], extract(abs2011_1km[[brt1$var.names]],cbind(datcombo$X,datcombo$Y)))
 
  writeRaster(mess1, filename="mess_brt_preds.tif", format="GTiff",overwrite=TRUE)

 
 
 # mean and coefficient of variation among the 20 brt models
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
 
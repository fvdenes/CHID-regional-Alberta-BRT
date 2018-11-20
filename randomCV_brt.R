randomCV_brt<-function(data1,nmodels,holdout){
  
  for(i in 1:nmodels){
    
    # sample data
    dat<-sample_frac(data1,size=1-holdout,replace=T,weight=data1$wt)
    
    #recalculate weights within sample
    weightraster <- rasterize(cbind(dat$X,dat$Y), pred_abs_2011[[1]], field=1)
    sampsum25 <- focal(weightraster, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
    dat$wt <- 1/extract(sampsum25,as.matrix(cbind(dat$X,dat$Y)))
    
    #fit BRT 
    x1 <- try(brt1 <- gbm.step(dat, gbm.y = "ABUND", 
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
                               family = "poisson", tree.complexity = 3, learning.rate = 0.025, bag.fraction = 0.5, offset=dat$logoffset, site.weights=dat$wt))
    if (class(x1) != "try-error") {
      w <-"D://CHID regional Alberta BRT/BRT_outputs/"
      z<-paste(w,"model",(1:nmodels)[i],"/",sep="")
      if (file.exists(z)==FALSE){dir.create(z)} 
      
      save(brt1,file=paste(z,speclist[j],"brtAB.R",sep=""))
      varimp <- as.data.frame(brt1$contributions)
      write.csv(varimp,file=paste(z,speclist[j],"varimp.csv",sep=""))
      cvstats <- t(as.data.frame(brt1$cv.statistics))
      write.csv(cvstats,file=paste(z,speclist[j],"cvstats.csv",sep=""))
      pdf(paste(z,speclist[j],"_plot.pdf",sep=""))
      gbm.plot(brt1,n.plots=32,smooth=TRUE, plot.layout = c(3,3),common.scale=T)
      dev.off()
      pdf(paste(z,speclist[j],"_plot.var-scale.pdf",sep=""))
      gbm.plot(brt1,n.plots=32,smooth=TRUE, plot.layout = c(3,3),common.scale=F,write.title = F)
      dev.off()
      pred_abs_2011<-brick("D:/CHID regional Alberta BRT/prediction dataset/abs2011_250m.grd")
      rast <- predict(pred_abs_2011, brt1, type="response", n.trees=brt1$n.trees)
      writeRaster(rast, filename=paste(z,speclist[j],"_pred1km",sep=""), format="GTiff",overwrite=TRUE)
      png(paste(z,speclist[j],"_pred1km.png",sep=""))
      plot(rast, zlim=c(0,1))
      points(dat$X, dat$Y, cex=0.05)
      dev.off()
      
      surveypoints<-SpatialPointsDataFrame(coords=cbind(dat$X,dat$Y),data=dat,proj4string = LCC)
      writeOGR(surveypoints,dsn=paste(z,"surveypoints.shp",sep=""),layer="surveypoints", driver="ESRI Shapefile")  
      
      rm(dat,brt1,varimp,cvstats,surveypoints,rast,abs2011_1km)
      gc()
    }
  }
  
}


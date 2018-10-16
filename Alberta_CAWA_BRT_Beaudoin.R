library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

load(("D:/CHID regional Alberta BRT/data_pack.RData"))

# list of tree spp for AB:"Species_Abie_Bal_v1"    
#"Species_Abie_Las_v1"                  
#"Species_Alnu_Rub_v1"     
#"Species_Betu_Pap_v1"                 
#"Species_Lari_Lar_v1"                                  
#"Species_Pice_Gla_v1"               
#"Species_Pice_Mar_v1"            
#"Species_Pinu_Ban_v1"                       
#"Species_Pinu_Con_v1"                   
#"Species_Popu_Bal_v1"                     
#"Species_Popu_Tre_v1"                                
#"Species_Pseu_Men_v1"                     
#"Species_Thuj_Pli_v1"                      
#"Species_Tsug_Het_v1"  

# define list 
w <-"D://CHID regional Alberta BRT/"

#for (j in 1:length(speclist)) { # to run for all spp. Make sure to uncomment the bracket at the end
j<-which(speclist=="CAWA") # to run only for CAWA

  specoff <- offl[offl$SPECIES==as.character(speclist[j]),]
  
  specdat2001 <- ABPC2001[ABPC2001$SPECIES == as.character(speclist[j]),] #n=423 for CAWA
  dat1 <- right_join(specdat2001[,c(1:5)],survey2001[,1:3],by=c("SS","PCODE","PKEY")) #n=31864
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
  d2001 <- left_join(s2001, dat2001, by=c("SS","PCODE")) 
  
  specdat2011 <- ABPC2011[ABPC2011$SPECIES == as.character(speclist[j]),] #n=989 for CAWA
  dat1 <- right_join(specdat2011[,c(1:5)],survey2011[,1:3],by=c("SS","PCODE","PKEY")) #n=26314
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
  s2011 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
  d2011 <- left_join(s2011, dat2011, by=c("SS","PCODE")) 
  
  datcombo <- rbind(d2001,d2011)
  datcombo$eco <- as.factor(datcombo$eco)
  
  x1 <- try(brt1 <- gbm.step(datcombo, gbm.y = 5, gbm.x = c(60,70,74,89,97,98,103,104,111,114,118,126,129,141,142), family = "poisson", tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.5, offset=datcombo$logoffset, site.weights=datcombo$wt))
  if (class(x1) != "try-error") {
    save(brt1,file=paste(w,speclist[j],"brtAB.R",sep=""))
    varimp <- as.data.frame(brt1$contributions)
    write.csv(varimp,file=paste(speclist[j],"varimp.csv",sep=""))
    cvstats <- t(as.data.frame(brt1$cv.statistics))
    write.csv(cvstats,file=paste(w,speclist[j],"cvstats.csv",sep=""))
    pdf(paste(w,speclist[j],"_plot.pdf",sep=""))
    gbm.plot(brt1,n.plots=9,smooth=TRUE, plot.layout = c(3,3))
    dev.off()
    rast <- predict(abs2011_1km, brt1, type="response", n.trees=brt1$n.trees)
    writeRaster(rast, filename=paste(speclist[j],"_pred1km",sep=""), format="GTiff",overwrite=TRUE)
    png(paste(speclist[j],"_pred1km.png",sep=""))
    plot(rast, zlim=c(0,1))
    points(datcombo$X, datcombo$Y, cex=0.05)
    dev.off()
  }
  
#}
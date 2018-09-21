
load("D:/CHID subunit delineation/pack_2016-12-01.Rdata")

# Preparing count dataset ####
## Combine CAWA counts, offsets and covariates in a single dataset
count <- YY[,"CAWA"]
offset <- OFF[,"CAWA"]

DAT$count<-count
DAT$CAWAoffset<-offset
DAT$YR<-as.factor(DAT$YR)
str(DAT)

## subset dataset to Alberta

ABDAT<-subset(DAT, DAT$JURS=="AB")

head(ABDAT)
## filter dataset to balance oversampled areas

#### A function to randomly sample from quadrants defined by latitude and longitude intervals.
sample_quad<-function(data,N=20,Xintervals=40,Yintervals=40,seed=123,onlypres=F){
  set.seed(seed)
  if(onlypres==F){
    maxlongitude<-max(data$X)
    minlongitude<-min(data$X)
    Xwidth<-(maxlongitude-minlongitude)/Xintervals
    Xlimits<-minlongitude+c(0,Xwidth*1:Xintervals)
    
    maxlatitude<-max(data$Y)
    minlatitude<-min(data$Y)
    Ywidth<-(maxlatitude-minlatitude)/Yintervals
    Ylimits<-minlatitude+c(0,Ywidth*1:Yintervals)
    
    ls<-array(NA,dim=c(Xintervals,Yintervals,N))
    for(i in 1:Xintervals){
      data2<-data[which(data$X>Xlimits[i] & data$X<Xlimits[i+1]),]
      for(j in 1:Yintervals){
        data3<-data2[which(data2$Y>Ylimits[j] & data2$Y<Ylimits[j+1]),]
        if(nrow(data3)>0){
          if(nrow(data3)<N){
            ls[i,j,1:nrow(data3)]<-row.names(data3)
          }
          else{
            ls[i,j,]<-row.names(data3[sample(nrow(data3),size=N),]) 
          }
        }
      }
    } 
  }
  
  else{
    data2<-data[which(data$count>0),]
    maxlongitude<-max(data2$X)
    minlongitude<-min(data2$X)
    Xwidth<-(maxlongitude-minlongitude)/Xintervals
    Xlimits<-minlongitude+c(0,Xwidth*1:Xintervals)
    
    maxlatitude<-max(data2$Y)
    minlatitude<-min(data2$Y)
    Ywidth<-(maxlatitude-minlatitude)/Yintervals
    Ylimits<-minlatitude+c(0,Ywidth*1:Yintervals)
    
    ls<-array(NA,dim=c(Xintervals,Yintervals,N))
    for(i in 1:Xintervals){
      data3<-data2[which(data2$X>Xlimits[i] & data2$X<Xlimits[i+1]),]
      for(j in 1:Yintervals){
        data4<-data3[which(data3$Y>Ylimits[j] & data3$Y<Ylimits[j+1]),]
        if(nrow(data4)>0){
          if(nrow(data4)<N){
            ls[i,j,1:nrow(data4)]<-row.names(data4)
          }
          else{
            ls[i,j,]<-row.names(data4[sample(nrow(data4),size=N),]) 
          }
        }
      }
    } 
  }
  
  ls<-c(ls)
  ls<-ls[-which(is.na(ls))]
  out<-match(x=ls,table = row.names(data))
}

#### Another function to plot selected point counts
plot_sampled<-function(obj,index,onlypres=T,show="count"){
  thin_cawa_mm<-obj[index,]
  mmsp <- SpatialPointsDataFrame(coords=cbind(thin_cawa_mm$X,thin_cawa_mm$Y),proj4string =  CRS("+proj=longlat +ellps=WGS84"),data=thin_cawa_mm)
  mmsp<-spTransform(mmsp,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
  if(onlypres==T){
    spplot(mmsp[which(mmsp$count>0),],show,do.log = F,cuts=4,legendEntries=as.character(1:4),
           key.space=list(x=0.5,y=0.9,corner=c(0,1))
           ,sp.layout=list(basemap)
           )
  }
  else{
    spplot(mmsp,show,do.log = F,legendEntries=as.character(0:4),
           key.space=list(x=0.02,y=0.3,corner=c(0,1))
           ,sp.layout=list(basemap)
           )
  }
}
basemap<-readOGR("D:/CHID subunit delineation/province_state_lcc.shp")

sample<-sample_quad(ABDAT, N=10,Xintervals=50,Yintervals=100,onlypres = F, seed=1111)


plot_sampled(ABDAT,index=sample, onlypres = F)

#comparing with all points
plot_sampled(ABDAT,index=1:nrow(ABDAT), onlypres = F)

par(mfrow=c(1,1))

ABDAT<-ABDAT[sample,]

## preparing prediction dataset AB: ####

load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-3.3.2.Rdata")
dat1<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-3.4.5.Rdata")
dat2<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.1.1.Rdata")
dat3<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.4.1.Rdata")
dat4<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.4.2.Rdata")
dat5<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-6.2.4.Rdata")
dat6<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-9.2.1.Rdata")
dat7<-dat



dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7)

levels(dat$JURS)
pred_data<-dat[which(dat$JURS=="ALBERTA"),]

rm(list=ls()[! ls() %in% c("ABDAT","sample","pred_data")]) 

save.image("D:/CHID regional Alberta BRT/data_pack.RData")

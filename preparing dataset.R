
load("D:/CHID subunit delineation/pack_2016-12-01.Rdata")

# Combine CAWA counts, offsets and covariates in a single dataset
count <- YY[,"CAWA"]
offset <- OFF[,"CAWA"]

DAT$CAWAcount<-count
DAT$CAWAoffset<-offset

str(DAT)

# subset dataset to Alberta

ABDAT<-subset(DAT, DAT$JURS=="AB")


# preparing prediction dataset AB:

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
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-9.3.1.Rdata")
dat8<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-6.2.6.Rdata")
dat9<-dat


dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9)

levels(dat$JURS)
pred_data<-dat[which(dat$JURS=="ALBERTA"),]

rm(list=ls()[! ls() %in% c("ABDAT","pred_data")]) 

save.image("D:/CHID regional Alberta BRT/data_pack.RData")

rm(list=ls())
library(zoo)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
dat<-read.csv("ct2015_weatherstationdata.csv")
austin<-read.csv("ct2015_mabry.csv")
austin$id<-"austin"
austin$TMAX[austin$TMAX==-9999]<-NA
austin$PRCP[austin$PRCP==-9999]<-0
austin<-austin[,3:7]
dat<-dat[-grep("AUSTIN", dat$STATION_NAME),]
dat$TMAX[dat$TMAX==-9999]<-NA
dat$PRCP[dat$PRCP==-9999]<-0
dat$id<-"corpuschristi"
cc<-dat[,3:7]
aus<-austin
aus$precip.mm<-aus$PRCP/10
aus$tmax.c<-aus$TMAX/10
aus$tmin.c<-aus$TMIN/10
cc$precip.mm<-cc$PRCP/10
cc$tmax.c<-cc$TMAX/10
cc$tmin.c<-cc$TMIN/10

wfc.drought<-c(0.00,0.00,0.01,0.00,0.00,0.00,0.05,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.12,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.00,0.00,0.00,0.00,0.00,0.00)
wfc.tempmax<-c(90.50,90.30,86.30,91.10,92.90,96.00,91.10,94.50,93.80,94.60,98.30,88.20,91.50,91.60,96.10,95.50,96.40,96.40,91.80,91.70,98.50,98.90,96.00,98.80,98.90,98.70,99.50,97.90,99.50,100.70,101.50,102.40,107.10,107.50,99.00,92.10,96.30,98.60,97.40,101.40,99.20,100.40,90.40,100.60,102.60,101.60,103.10,105.40)
wfc.tempmin<-c(69.95,72.83,72.11,74.44,75.64,75.52,67.67,60.62,65.93,65.99,67.37,68.09,59.00,70.70,66.35,73.42,70.07,73.66,72.53,76.12,74.32,76.05,70.54,74.08,73.49,75.99,75.40,75.75,75.28,76.17,76.34,74.91,73.85,71.39,70.07,67.49,66.11,64.49,63.12,65.51,67.55,71.92,75.50,73.48,70.85,74.02,73.12,70.13)
wfc.tempmax<-(wfc.tempmax - 32) * 5/9
wfc.tempmin<-(wfc.tempmin - 32) * 5/9
wfc.drought<-wfc.drought*25.4
aus$month<-as.numeric(substring(aus$DATE,5,6))
cc$month<-as.numeric(substring(cc$DATE,5,6))
aus$day<-as.numeric(substring(aus$DATE,7,8))
cc$day<-as.numeric(substring(cc$DATE,7,8))
aus$year<-as.numeric(substring(aus$DATE,1,4))
cc$year<-as.numeric(substring(cc$DATE,1,4))
cc<-cc[cc$month==6 | cc$month==5 & cc$day>25 | cc$month==7 & cc$day>7,]
aus<-aus[aus$month==6 | aus$month==5 & aus$day>25 | aus$month==7 & aus$day>7,]
library(plyr)
aus.summary <- ddply(aus, c("year"), summarise,
               sum.precip    = sum(precip.mm),
               mean.tmax.c = mean(tmax.c),
               mean.tmin.c = mean(tmin.c))
cc.summary <- ddply(cc, c("year"), summarise,
               sum.precip    = sum(precip.mm),
               mean.tmax.c = mean(tmax.c),
               mean.tmin.c = mean(tmin.c))
cc.summary<-cc.summary[complete.cases(cc.summary),]
aus.summary<-aus.summary[complete.cases(aus.summary),]
par(mfrow=c(2,1))
aus.density.precip<-density(log10(aus.summary$sum.precip),na.rm=T)
plot(10^(aus.density.precip$x),aus.density.precip$y)
hist(aus.summary$sum.precip,breaks=1000)

cc.density.precip<-density(log10(cc.summary$sum.precip),na.rm=T)
plot(10^(cc.density.precip$x),cc.density.precip$y)
hist(cc.summary$sum.precip,breaks=100)

cc.density.precip<-density(sw.cc.precip,na.rm=T)
aus.density.tmax<-density(sw.aus.maxtemp,na.rm=T)
cc.density.tmax<-density(sw.cc.maxtemp,na.rm=T)


aus<-aus[aus$month %in% c(5,6,7),]
cc<-cc[-cc$month==5 & cc$day<25,]
aus<-aus[-aus$month==5 & aus$day<25,]
cc<-cc[-cc$month==7 & cc$day>7,]
aus<-aus[-aus$month==7 & aus$day<7,]
#sliding window analysis of drought intensity:
mean.drought.precip<-mean(wfc.drought)
mean.drought.maxtemp<-mean(wfc.tempmax)
sw.aus.precip<-rollapply(aus$precip.mm,width=40, by=7, function(x) mean(x, na.omit=T))
sw.aus.maxtemp<-rollapply(aus$tmax.c,width=40, by=7, function(x) mean(x, na.omit=T))
sw.cc.precip<-rollapply(cc$precip.mm,width=40, by=7, function(x) mean(x, na.omit=T))
sw.cc.maxtemp<-rollapply(cc$tmax.c,width=40, by=7, function(x) mean(x, na.omit=F))

pdf("ct2015_climatecomps.pdf")
par(mfrow=c(2,1))
plot(cc.density.precip, col="blue", xlim=c(0,10), main="mean precip (mm) May-July", bty="n")
lines(aus.density.precip, col="red")
abline(v=mean.drought.precip, col="black")
legend("topright",c("corpuschristi","austin","wfc drought"),col=c("blue","red","black"),lty=1)
plot(cc.density.tmax, col="blue", xlim=c(22,40), main="mean of max temp (c) May-July", bty="n")
lines(aus.density.tmax, col="red")
abline(v=mean.drought.maxtemp, col="black")
legend("topleft",c("corpuschristi","austin","wfc drought"),col=c("blue","red","black"),lty=1)
dev.off()

drought.beg<-which(aus$DATE=="20130527")
drought.end<-which(aus$DATE=="20130705")


wfc<-weather[weather$id=="austin",]
wfc[drought.beg:drought.end,"precip.mm"]<-wfc.drought[1:40]
wfc[drought.beg:drought.end,"tmax.c"]<-wfc.tempmax[1:40]
wfc[drought.beg:drought.end,"tmin.c"]<-wfc.tempmin[1:40]
hr<-hargreaves(Tmin=wfc$tmin.c, Tmax=wfc$tmax.c, lat=30.19)
sp<-spei(wfc$precip.mm-hr, scale=.1)
wfc$spei<-sp[[2]]
drought.intensity<-min(wfc$spei[drought.beg:drought.end])

hr.aus<-hargreaves(Tmin=aus$tmin.c, Tmax=aus$tmax.c, lat=30.19)
sp.aus<-spi(aus$precip.mm, scale=.1)
aus$spei<-sp.aus[[2]]
aus$spei[drought.beg:drought.end]

hr.cc<-hargreaves(Tmin=cc$tmin.c, Tmax=cc$tmax.c, lat=27.5)
sp.cc<-spi(cc$precip.mm, scale=.1)
cc$spei<-sp.cc[[2]]
cc$spei[drought.beg:drought.end]


#read in experimental design data
library(reshape)
harvest<-read.csv("ct2015_wfc_expdesigndata.csv")
harvest<-harvest[,c("PlantID_FINAL","Day_July","MD_LWP_July","PD_LWP_July","VWC_July")]
harvest<-harvest[-grep("FH", harvest$PlantID_FINAL),]
harvest<-harvest[complete.cases(harvest),]
harvest$predawn<-(-1*(harvest$PD_LWP_July/10))
harvest$midday<-(-1*(harvest$MD_LWP_July/10))
harvest<-harvest[,c("PlantID_FINAL","Day_July","predawn","midday")]
harvest<-harvest[harvest$PlantID_FINAL %in% c("FIL2","HAL2"),]
harvest.long<-melt(harvest, id.vars=c("PlantID_FINAL","Day_July"))
harvest.long$transformed.lwp<-log10(harvest.long$value*-1)
harvest.long$Day_July[harvest.long$Day_July==1]<-"Drought"
harvest.long$Day_July[harvest.long$Day_July==2]<-"Recovery"
cdata <- ddply(harvest.long, c("PlantID_FINAL","Day_July","variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)
cdata<-cdata[cdata$variable=="midday",]
harvest.long<-harvest.long[harvest.long$variable=="midday",]
cdata$measure<-paste(cdata$Day_July,cdata$variable,sep="_")
harvest.long$measure<-paste(harvest.long$Day_July,cdata$variable,sep="_")
limits <- with(cdata, aes(ymax = mean + se, ymin=mean - se))
cdata$geno<-as.character(cdata$PlantID_FINAL)
harvest.long$geno<-as.character(harvest.long$PlantID_FINAL)

ggplot(cdata, aes(x=Day_July,y=mean,col=geno,group=geno))+
  geom_line() + geom_errorbar(limits, width=0.2)+theme_classic()+
  scale_y_continuous(limits=c(-4,0))
  
harvest.midday<-harvest.long[harvest.long$variable=="midday",]
summary(aov(value~PlantID_FINAL+Day_July+PlantID_FINAL*Day_July, data=harvest.long))

weather<-read.csv("ct2015_wfc_weatherdata.csv")
weather<-weather[,c("Date","Day","Air.T.F.Max","Soil.T.F.Max","Total.Rain..in.","Daily.Rain..in.")]
colnames(weather)<-c("date","day","air.max.temp","soil.max.temp","total.precip","daily,precip")
weather<-weather[complete.cases(weather),]
weather$air.max.temp<-(weather$air.max.temp-32)*(5/9)
weather$soil.max.temp<-(weather$soil.max.temp-32)*(5/9)
weather$month<-substring(weather$date,1,1)
weather<-weather[weather$month %in% c(4,5,6,7),]
par(mfrow=c(3,1))
plot(weather$total.precip, type="l",bty="n", col="blue")
plot(weather$air.max.temp, type="l",bty="n", ylim=c(15,50),col="red")
plot(weather$soil.max.temp, type="l",bty="n", ylim=c(15,50), col="red", lty=2)

test<-harvest.long[harvest.long$Day_July=="Drought" & harvest.long$variable=="midday",]

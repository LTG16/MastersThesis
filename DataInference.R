library(openair)
library(imputeTS)
library(rgdal)
library(raster)
library(geosphere)
########This File contains the code for the Data Inference
### Downloading the UB stations
dataUB<-importKCL(site =c("BL0","CT3","EA7","IS6","KC1", "LW1","WA2","WM0"), year = 2008:2009,
                  pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
dataUB<-timeAverage(dataUB, period="day",type=c("site","code","site.type"))
timePlot(dataUB,type="site")
sites<-unique(dataUB$code)

df <- unique(dataUB[,c(2,6,7)])
df1 <- df
d<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/1000

#### Preparing the UB-data

polution.data <- dataUB[ which(dataUB$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB[ which(dataUB$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data <- as.data.frame(polution.data)

data<-matrix(nrow=233,ncol=8)

data[,1] <- na.interpolation(pol.data[351:583,1])
data[,2] <- na.interpolation(pol.data[351:583,2])
data[,3] <- na.interpolation(pol.data[351:583,3])
data[,4] <- na.interpolation(pol.data[351:583,4])
data[,5] <- na.interpolation(pol.data[351:583,5])
data[,6] <- na.interpolation(pol.data[351:583,6])
data[,7] <- na.interpolation(pol.data[351:583,7])
data[,8] <- na.interpolation(pol.data[351:583,8])
View(data)


##########Downloading the data for the assumed to be fixed stations and preparing it

dataUB2<-importKCL(site =c("KC1", "LW1","WA2","WM0"), year = 2008:2009,
                  pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
#timePlot(dataUB,type="site")
dataUB2<-timeAverage(dataUB2, period="day",type=c("site","code","site.type"))
timePlot(dataUB2,type="site")
sites<-unique(dataUB2$code)


polution.data <- dataUB2[ which(dataUB2$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB2[ which(dataUB2$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data <- as.data.frame(polution.data)


datastart<-matrix(nrow=233,ncol=4)
datastart[,1] <- na.interpolation(pol.data[351:583,1])
datastart[,2] <- na.interpolation(pol.data[351:583,2])
datastart[,3] <- na.interpolation(pol.data[351:583,3])
datastart[,4] <- na.interpolation(pol.data[351:583,4])
datastart<-log(datastart)


df3 <- unique(dataUB2[,c(2,6,7)])
df1 <- df3
d2<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000
df <- unique(dataUB[,c(2,6,7)])

##############Plotting the UB stations and the assumed to be fixed stations
library(rgdal)
library(raster)
par(mfrow=c(1,1),mar=c(0,0,0,0))
lnd <- readOGR(dsn = path.expand("~/Desktop/Masterarbeit/data"), layer="london_sport")
plot(lnd,col="gray95") # plot the lnd object (not shown)



crs(lnd)
coordinates(df) <- ~longitude + latitude
crs(df) <- CRS("+proj=longlat +datum=WGS84")
df <- spTransform(df, crs(lnd))
coordinates(df3) <- ~longitude + latitude
crs(df3) <- CRS("+proj=longlat +datum=WGS84")
df3 <- spTransform(df3, crs(lnd))


points(df$longitude, df$latitude, col = "red", pch =17,lwd=3,cex=1.2)

points(df3$longitude, df3$latitude, col = "green", pch =17,lwd=3,cex=1.2)
legend("bottomright",col=c("red","green"),legend=c("Non-fixed Stations","Fixed Stations")
       ,pch=c(17,17),bty="n")








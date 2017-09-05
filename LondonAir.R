#####Master thesis descriptive analysis of airquality data NOX 2015

library(openair)
library(imputeTS)

####Here we import all sites in London from 2015, which contain pollutant NOX 
#### Since there is no function which automatically downloades the data of the measurement
#### sites in London, they have to be typed in manually. 
datakclraw<-importKCL(site =c("BL0","CD1","CT3","CT6",
                           "EA6","EA7","GB6","GN0","GN2","GN3",
                           "GR4","GR5","GR7","GR8","GR9",
                           "HG1","HI0","HK6","IS2","IS6","KC1","KC2","KC3",
                           "KC4","KC5",
                           "RI1","RI2","TD0","TH2","TH4",
                           "MY1","WM0","WA2","WA7","LW1"), 
                      year = 2008:2009, pollutant = "nox", met = FALSE,units = "mass",
                   extra = FALSE, meta = TRUE)
#### Here the hourly data is converted to daily data
datakcl1<-timeAverage(datakclraw, period="day",type=c("site","code","site.type"))
datakcl<-as.data.frame(datakcl1)
##############################
##############################
############################## Plots for Data~Overview Chapter in Masterthesis

####timeplot shows the different timeseries of the sites
timePlot(datakcl, pollutant="nox", type="site", col="blue")

#### Now we plot all the sites into a map of london, to check where they are etc. 
library(rgdal)
library(raster)
par(mfrow=c(1,1),mar=c(0,0,0,0))
lnd <- readOGR(dsn = path.expand("~/Desktop/Masterarbeit/data"), layer="london_sport")
plot(lnd,col="gray95") # plot the lnd object (not shown)


df <- unique(datakcl[,c(2,6,7)])
df<-df[-c(6),]
df1 <- df
crs(lnd)

coordinates(df) <- ~longitude + latitude
crs(df) <- CRS("+proj=longlat +datum=WGS84")
df <- spTransform(df, crs(lnd))

points(df$longitude, df$latitude, col = "blue", pch =17,lwd=2)

#################################################################

### Now we want to calculate the distance in 10 km and then check, if the correlation 
# of NOX decreases by increasing distance of the different sites. 
library(geosphere)

d<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000

sites <- c("BL0","CD1","CT3","CT6",
           "EA6","EA7","GB6","GN0","GN2","GN3",
           "GR4","GR5","GR7","GR8","GR9",
           "HG1","HI0","HK6","IS2","IS6","KC1","KC2","KC3",
           "KC4","KC5",
           "RI1","RI2","TD0","TH2","TH4",
           "MY1","WM0","WA2","WA7","LW1")

polution.data <- datakcl[ which(datakcl$code==sites[1] ), 5]
for(i in 2:35){
  polution.data  <- cbind(polution.data,datakcl[ which(datakcl$code==sites[i] ), 5])
}
for(r in 1:35){
  colnames(polution.data)[r] <- sites[r]
}
pol.data <- as.data.frame(polution.data)
test <- corPlot(pol.data)

corr.vals <- test[2]

corr <- corr.vals[[1]]$cor
corr.x <-corr.vals[[1]]$x
corr.y <-corr.vals[[1]]$y

cx<-c()
cy<-c()

for(i in 1:length(corr)){
  cx <- cbind(cx,as.vector(corr.x[i]))
  cy <- cbind(cy,as.vector(corr.y[i]))
}
korrelation.matrix <- as.data.frame(cbind(as.vector(cx),as.vector(cy),corr))
korrelation.matrix <- korrelation.matrix[order(korrelation.matrix$V1,korrelation.matrix$V2),]

dist.vals <- as.vector(d)

sites <- rep(df1$code,35)

sites2<-c()
for(i in 1:35){
  sites2 <- c(sites2,rep(sites[i],35))
}

distanz <- as.data.frame(cbind(sites,sites2,dist.vals))
distanz <- distanz[order(distanz$sites,distanz$sites2),]

dist.sub <- distanz
korr.sub <-korrelation.matrix
x <- as.numeric(as.vector(dist.sub$dist.vals))
y <- as.numeric(as.vector(korr.sub$corr))

plot(x,y,xlab="Distance in 10 [km]",ylab="Correlation",col="darkblue",ylim=c(-0.2,1.2),xaxt="n")
axis(1, at=seq(0,4,by=0.5))
model <- lm(y~x)
abline(model,col = "red")
predict(model,1)
abline(v=c(0,0.5,1,1.5,2), h=c(-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2), col="gray", lty=3)
summary(model)

# since there is no correlation in dependency of distance we look at the different groups
#####################################################################################


### we first split the dataset into the different groups and look at the different timeplots of each group
X<-split(datakcl,datakcl$site.type)

X1<-X$Kerbside
timePlot(X1,pollutant="nox", type="site", col="blue")
X2<-X$Industrial
timePlot(X2,pollutant="nox", type="site", col="blue")
X3<-X$Roadside
timePlot(X3,pollutant="nox", type="site", col="blue")
X4<-X$Suburban
timePlot(X4,pollutant="nox", type="site", col="blue")
X5<-X$`Urban Background`
timePlot(X5,pollutant="nox", type="site", col="blue")


par(mfrow=c(1,2))
timePlot(X3,pollutant="nox", type="site", col="blue")
timePlot(X5,pollutant="nox", type="site", col="blue")



###################Now we look at the average correlation of the different site-types 
# Furthermore we transform the data, so that we can work with it better. 



####################   Kerbside
sites1 <- unique(X1$code)
polution.data1 <- X1[ which(X1$code==sites1[1] ), 5]
for(i in 2:length(sites1)){
  polution.data1  <- cbind(polution.data1,X1[which(X1$code==sites1[i] ), 5])
}
pol.data1 <- as.data.frame(polution.data1)
test1 <- corPlot(pol.data1)

corr.vals1 <- test1[2]

corr1 <- corr.vals1[[1]]$cor

mean(corr1[which(corr1!=1)])   ###0.3409132


####################   INDUSTRIAL   
sites2 <- unique(X2$code)
polution.data2 <- X2[ which(X2$code==sites2[1] ), 5]
for(i in 2:length(sites2)){
  polution.data2  <- cbind(polution.data2,X2[which(X2$code==sites2[i] ), 5])
}
pol.data2 <- as.data.frame(polution.data2)
test2 <- corPlot(pol.data2)

corr.vals2 <- test2[2]

corr2 <- corr.vals2[[1]]$cor

mean(corr2[which(corr2!=1)])   ###0.6882196







####################   Roadside
sites3 <- unique(X3$code)
polution.data3 <- X3[ which(X3$code==sites3[1] ), 5]
for(i in 2:length(sites3)){
  polution.data3  <- cbind(polution.data3,X3[which(X3$code==sites3[i] ), 5])
}
pol.data3 <- as.data.frame(polution.data3)
test3 <- corPlot(pol.data3)

corr.vals3 <- test3[2]

corr3 <- corr.vals3[[1]]$cor

mean(corr3[which(corr3!=1)])   #0.53353


for(r in 1:length(sites3)){
  colnames(polution.data3)[r] <- sites3[r]
}
test3 <- corPlot(pol.data3)

corr.vals3 <- test3[2]

corr3 <- corr.vals3[[1]]$cor
corr.x3 <-corr.vals3[[1]]$x
corr.y3 <-corr.vals3[[1]]$y

cx3<-c()
cy3<-c()

for(i in 1:length(corr3)){
  cx3 <- cbind(cx3,as.vector(corr.x3[i]))
  cy3 <- cbind(cy3,as.vector(corr.y3[i]))
}
korrelation.matrix3 <- as.data.frame(cbind(as.vector(cx3),as.vector(cy3),corr3))
korrelation.matrix3 <- korrelation.matrix3[order(korrelation.matrix3$V1,korrelation.matrix3$V2),]
df <- unique(X3[,c(2,6,7)])
df<-df[-3,]
df1 <- df
d1<-distm(df1[,c(3,2)], df1[,c(3,2)], fun = distHaversine)/10000
#as.matrix(dtemp<-dist(df1[,c(3,2)]))

dist.vals3 <- as.vector(d1)

sites3 <- rep(df1$code,length(unique(X3$code)))
sites33<-c()
for(i in 1:length(unique(X3$code))){
  sites33 <- c(sites33,rep(sites3[i],length(unique(X3$code))))
}


distanz3 <- as.data.frame(cbind(sites3,sites33,dist.vals3))
distanz3 <- distanz3[order(distanz3$sites3,distanz3$sites33),]

dist.sub <- distanz3

korr.sub <-korrelation.matrix3
x1 <- as.numeric(as.vector(dist.sub$dist.vals3))
y1 <- as.numeric(as.vector(korr.sub$corr3))

plot(x1,y1,xlab="Distance in 10 [km]",ylab="Correlation",col="darkblue",ylim=c(-0.2,1.2),xaxt="n",
     main="Roadside")
axis(1, at=seq(0,3,by=0.5))
model <- lm(y1~x1)
abline(model,col = "red")
predict(model,1)
abline(v=c(0,0.5,1,1.5,2,2.5), h=c(-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2), col="gray", lty=3)
summary(model)



####################   SUBURBAN
sites4 <- unique(X4$code)
polution.data4 <- X4[ which(X4$code==sites4[1] ), 5]
for(i in 2:length(sites4)){
  polution.data4  <- cbind(polution.data4,X4[which(X4$code==sites4[i] ), 5])
}
pol.data4 <- as.data.frame(polution.data4)
test4 <- corPlot(pol.data4)

corr.vals4 <- test4[2]

corr4 <- corr.vals4[[1]]$cor

mean(corr4[which(corr4!=1)])    ###0.6372444


####################   URBAN BACKGROUND
sites5 <- unique(X5$code)
polution.data5 <- X5[ which(X5$code==sites5[1] ), 5]
for(i in 2:length(sites5)){
  polution.data5  <- cbind(polution.data5,X5[which(X5$code==sites5[i] ), 5])
}
pol.data5 <- as.data.frame(polution.data5)
test5 <- corPlot(pol.data5)

corr.vals5 <- test5[2]

corr5 <- corr.vals5[[1]]$cor

mean(corr5[which(corr5!=1)])   ### 0.6996025

for(r in 1:length(sites5)){
  colnames(polution.data5)[r] <- sites5[r]
}
test5 <- corPlot(pol.data5)

corr.vals5 <- test5[2]

corr5 <- corr.vals5[[1]]$cor
corr.x5 <-corr.vals5[[1]]$x
corr.y5 <-corr.vals5[[1]]$y

cx5<-c()
cy5<-c()

for(i in 1:length(corr5)){
  cx5 <- cbind(cx5,as.vector(corr.x5[i]))
  cy5 <- cbind(cy5,as.vector(corr.y5[i]))
}
korrelation.matrix5 <- as.data.frame(cbind(as.vector(cx5),as.vector(cy5),corr5))
korrelation.matrix5 <- korrelation.matrix5[order(korrelation.matrix5$V1,korrelation.matrix5$V2),]
df <- unique(X5[,c(2,6,7)])
df1 <- df
d1 <- distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000
d1<-as.matrix(d1)

dist.vals5 <- as.vector(d1)

sites5 <- rep(df1$code,length(unique(X5$code)))
sites6<-c()
for(i in 1:length(unique(X5$code))){
  sites6 <- c(sites6,rep(sites5[i],length(unique(X5$code))))
}


distanz5 <- as.data.frame(cbind(sites5,sites6,dist.vals5))
distanz5 <- distanz5[order(distanz5$sites5,distanz5$sites6),]

dist.sub <- distanz5

korr.sub <-korrelation.matrix5
x <- as.numeric(as.vector(dist.sub$dist.vals5))
y <- as.numeric(as.vector(korr.sub$corr5))
plot(x,y,xlab="Distance in 10 [km]",ylab="Correlation",col="darkblue",ylim=c(0.2,1.2),xaxt="n",main=
       "Urban Background")
axis(1, at=seq(0,2.5,by=0.5))
model <- lm(y~x)
abline(model,col = "red")
predict(model,1)
abline(v=c(0,0.5,1,1.5,2,2.5), h=c(0.2,0.4,0.6,0.8,1.0,1.2), col="gray", lty=3)

summary(model)

########################

n <- ncol(pol.data)
B <- diag(n)
for(i in 1:n){
  for(k in 1:n){
    if(B[i,k]!=1){
      if((sites[i]%in% sites1 & sites[k]%in% sites1 )| 
         (sites[i]%in% sites2 & sites[k]%in% sites2 )|
         (sites[i]%in% sites3 & sites[k]%in% sites3 )|
         (sites[i]%in% sites4 & sites[k]%in% sites4 )|
         (sites[i]%in% sites5 & sites[k]%in% sites5 )){
        if(sites[k] %in% sites1){
          B[i,k] <- mean(corr1)
        }
        if(sites[k] %in% sites2){
          B[i,k] <- mean(corr2)
        }
        if(sites[k] %in% sites3){
          B[i,k] <- mean(corr3)
        }
        if(sites[k] %in% sites4){
          B[i,k] <- mean(corr4)
        }
        if(sites[k] %in% sites5){
          B[i,k] <- mean(corr5)
        }
      }
    }
  }
}



########################################## Mean for table
polmean11<-colMeans(pol.data1,na.rm=TRUE)
polmean22<-colMeans(pol.data2,na.rm=TRUE)
polmean33<-colMeans(pol.data3,na.rm=TRUE)
polmean44<-colMeans(pol.data4,na.rm=TRUE)
polmean55<-colMeans(pol.data5,na.rm=TRUE)


mean(polmean11)
mean(polmean22)
mean(polmean33)
mean(polmean44)
mean(polmean55)

##########################################  Variance for table



polvar11<-apply(pol.data1,2,FUN=var,na.rm=TRUE)
polvar22<-apply(pol.data2,2,FUN=var,na.rm=TRUE)
polvar33<-apply(pol.data3,2,FUN=var,na.rm=TRUE)
polvar44<-apply(pol.data4,2,FUN=var,na.rm=TRUE)
polvar55<-apply(pol.data5,2,FUN=var,na.rm=TRUE)


mean(polvar11)
mean(polvar22)
mean(polvar33)
mean(polvar44)
mean(polvar55)


##########################################################


################Average N0x concentration of the different types over the time
# This plot is not part of the thesis


plot(polmean1,type="l",xlab="Day",ylab="NOx",ylim=c(0,700),xaxt="n")
axis(1,at=seq(0,750,by=50))
lines(polmean2,col="blue")
lines(polmean3,col="red")
lines(polmean4,col="green")
lines(polmean5,col="purple")
abline(v=c(0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750), h=c(0,100,200,300,400,500,600,700), col="gray", lty=3)

legend("topright",c("Kerbside","Industrial","Roadside","Suburban","Urban Background"),
       lty = c(1, 1,1,1,1),col=c("black","blue","red","green","purple"),cex=0.55)
####################################################

library(zoo)

############################################
############################################
###########################  Urban Background Plots and Data
library(openair)
### We download the UB Dataset and transform the data so that we can use it easier
dataUB<-importKCL(site =c("BL0","CT3","EA7","IS6","KC1", "LW1","WA2","WM0"), year = 2008:2009,
                  pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
dataUB<-timeAverage(dataUB, period="day",type=c("site","code","site.type"))
timePlot(dataUB,type="site")
sites<-unique(dataUB$code)
polution.data <- dataUB[ which(dataUB$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB[ which(dataUB$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data <- as.data.frame(polution.data)

####Checking how many NAs there are, and which we need to cut out or interpolate
length(which(is.na(pol.data[,1]))) 
length(which(is.na(pol.data[,2]))) 
length(which(is.na(pol.data[,3]))) 
length(which(is.na(pol.data[,4]))) 
length(which(is.na(pol.data[,5]))) 
length(which(is.na(pol.data[,6]))) 
length(which(is.na(pol.data[,7]))) 
length(which(is.na(pol.data[,8]))) 

########################   Interpolation of 1,2,7,8
na.interpolation(pol.data[,1])
na.interpolation(pol.data[,2])
na.interpolation(pol.data[,7])
na.interpolation(pol.data[,8])

## Plot that shows missing values
par(mfrow=c(2,4))
plotNA.distribution(pol.data[,1],colPoints="",colBackgroundMV ="red",col="darkblue",main="BLO",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,2],colPoints="",colBackgroundMV ="red",col="darkblue",main="CT3",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,3],colPoints="",colBackgroundMV ="red",col="darkblue",main="EA7",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,4],colPoints="",colBackgroundMV ="red",col="darkblue",main="IS6",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,5],colPoints="",colBackgroundMV ="red",col="darkblue",main="KC1",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,6],colPoints="",colBackgroundMV ="red",col="darkblue",main="LW1",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,7],colPoints="",colBackgroundMV ="red",col="darkblue",main="WA2",xlab="Day",ylab="N0x Concentration")
plotNA.distribution(pol.data[,8],colPoints="",colBackgroundMV ="red",col="darkblue",main="WM0",xlab="Day",ylab="N0x Concentration")

#### use either 351:534 or 351:583 if interpolation of 5 values in a row, is ok 
which(is.na(pol.data[221:632,3]))
which(is.na(pol.data[306:583,4]))  
which(is.na(pol.data[310:583,5]))
which(is.na(pol.data[351:632,6]))

# Step 1: Perform imputation for x using na.interpolation
ts.imp1 <- na.interpolation(pol.data[351:583,1])
ts.imp2 <- na.interpolation(pol.data[351:583,2])
ts.imp3 <- na.interpolation(pol.data[351:583,3])
ts.imp4 <- na.interpolation(pol.data[351:583,4])
ts.imp5 <- na.interpolation(pol.data[351:583,5])
ts.imp6 <- na.interpolation(pol.data[351:583,6])
ts.imp7 <- na.interpolation(pol.data[351:583,7])
ts.imp8 <- na.interpolation(pol.data[351:583,8])


# Step 2: Visualize the imputed values in the time series
plotNA.imputations(pol.data[351:583,1], ts.imp1)
plotNA.imputations(pol.data[351:583,2], ts.imp2)
plotNA.imputations(pol.data[351:583,3], ts.imp3)
plotNA.imputations(pol.data[351:583,4], ts.imp4)
plotNA.imputations(pol.data[351:583,5], ts.imp5)
plotNA.imputations(pol.data[351:583,6], ts.imp6)
plotNA.imputations(pol.data[351:583,7], ts.imp7)
plotNA.imputations(pol.data[351:583,8], ts.imp8)




####################ACF Plots
par(mfrow=c(2,4))
acf(pol.data[,1],na.action=na.pass,main="BLO")
acf(pol.data[,2],na.action=na.pass,main="CT3")
acf(pol.data[,3],na.action=na.pass,main="EA7")
acf(pol.data[,4],na.action=na.pass,main="IS6")
acf(pol.data[,5],na.action=na.pass,main="KC1")
acf(pol.data[,6],na.action=na.pass,main="LW1")
acf(pol.data[,7],na.action=na.pass,main="WA2")
acf(pol.data[,8],na.action=na.pass,main="WM0")


###################Histogramm of Log data
par(mfrow=c(2,4))
hist(log(ts.imp1),xlab="log(NOx)",main="BLO",col="grey")
hist(log(ts.imp2),xlab="log(NOx)",main="CT3",col="grey")
hist(log(ts.imp3),xlab="log(NOx)",main="EA7",col="grey")
hist(log(ts.imp4),xlab="log(NOx)",main="IS6",col="grey")
hist(log(ts.imp5),xlab="log(NOx)",main="KC1",col="grey")
hist(log(ts.imp6),xlab="log(NOx)",main="LW1",col="grey")
hist(log(ts.imp7),xlab="log(NOx)",main="WA2",col="grey")
hist(log(ts.imp8),xlab="log(NOx)",main="WMO",col="grey")

###################Histogramm of raw data

par(mfrow=c(2,4))
hist(ts.imp1,xlab="NOx",main="BLO",col="grey")
hist(ts.imp2,xlab="NOx",main="CT3",col="grey")
hist(ts.imp3,xlab="NOx",main="EA7",col="grey")
hist(ts.imp4,xlab="NOx",main="IS6",col="grey")
hist(ts.imp5,xlab="NOx",main="KC1",col="grey")
hist(ts.imp6,xlab="NOx",main="LW1",col="grey")
hist(ts.imp7,xlab="NOx",main="WA2",col="grey")
hist(ts.imp8,xlab="NOx",main="WMO",col="grey")


###################################################################

######Plot that shows all stations in London and the Urban Background Stations

par(mfrow=c(1,1),mar=c(0,0,0,0))
plot(lnd,col="gray95") # plot the lnd object (not shown)


df <- unique(datakcl[,c(2,6,7)])
df <- df[c(-6),]
crs(lnd)

coordinates(df) <- ~longitude + latitude
crs(df) <- CRS("+proj=longlat +datum=WGS84")
df <- spTransform(df, crs(lnd))

points(df$longitude, df$latitude, col = "blue", pch =17,lwd=2)

df <- unique(dataUB[,c(2,6,7)])
df1 <- df
crs(lnd)

coordinates(df) <- ~longitude + latitude
crs(df) <- CRS("+proj=longlat +datum=WGS84")
df <- spTransform(df, crs(lnd))

points(df$longitude, df$latitude, col = "red", pch =17,lwd=2)
legend("bottomright",col=c("blue","red"),legend=c("Stations","Urban Background Stations")
       ,pch=c(17,17),bty="n",cex=0.8)

#########################################################################################





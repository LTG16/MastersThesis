### In this file we exten the F matrix for the four different non-fixed stations
## We do this here for the Kalman-Filter estimates
F<-matrix(c(0.4675871,0.3672238,-0.9784721,-0.4740237,0.07150804,0.72246269,
            0.28619770,0.16485973,0.3680380,1.7403075,-0.1698203,-0.1421445,
            0.00529898,-1.89673653,1.81604593,1.39362019),4,4)


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
df3 <- unique(dataUB2[,c(2,6,7)])
df1 <- df3
d2<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000


as.vector(F)
as.vector(d2)
plot(x=as.vector(d2),y=as.vector(F))
model<-lm(as.vector(F)~as.vector(d2))
abline(model)
summary(model)
#### Plot for the model
plot(x=as.vector(d2),y=as.vector(F),xlab="Distance in 10 [km]",ylab="F-value",col="darkblue")
model<-lm(as.vector(F)~as.vector(d2))
abline(model,col = "red")
abline(v=c(0,0.5,1,1.5), h=c(-2,-1,0,1), col="gray", lty=3)
summary(model)





####Extending F for BL0

##############################################################
d3
F3<-matrix(0,nrow=5,ncol=5)
F3[2:5,2:5] <- F
F3[1,1] <- 0.3961-0.2286*d3[1,1]
F3[1,2] <- 0.3961-0.2286*d3[1,2]
F3[2,1] <- 0.3961-0.2286*d3[2,1]
F3[1,3] <- 0.3961-0.2286*d3[1,3]
F3[3,1] <- 0.3961-0.2286*d3[3,1]
F3[1,4] <-0.3961-0.2286*d3[1,4]
F3[4,1] <- 0.3961-0.2286*d3[4,1]
F3[1,5] <- 0.3961-0.2286*d3[1,5]
F3[5,1] <- 0.3961-0.2286*d3[5,1]
##############################################################
####Extending F for CT3

##############################################################
d4
F4<-matrix(0,nrow=5,ncol=5)
F4[2:5,2:5] <- F
F4[1,1] <- 0.3961-0.2286*d4[1,1]
F4[1,2] <- 0.3961-0.2286*d4[1,2]
F4[2,1] <- 0.3961-0.2286*d4[2,1]
F4[1,3] <- 0.3961-0.2286*d4[1,3]
F4[3,1] <- 0.3961-0.2286*d4[3,1]
F4[1,4] <- 0.3961-0.2286*d4[1,4]
F4[4,1] <- 0.3961-0.2286*d4[4,1]
F4[1,5] <- 0.3961-0.2286*d4[1,5]
F4[5,1] <- 0.3961-0.2286*d4[5,1]
##############################################################
####Extending F for EA7

##############################################################
d5
F5<-matrix(0,nrow=5,ncol=5)
F5[2:5,2:5] <- F
F5[1,1] <- 0.3961-0.2286*d5[1,1]
F5[1,2] <- 0.3961-0.2286*d5[1,2]
F5[2,1] <- 0.3961-0.2286*d5[2,1]
F5[1,3] <- 0.3961-0.2286*d5[1,3]
F5[3,1] <- 0.3961-0.2286*d5[3,1]
F5[1,4] <- 0.3961-0.2286*d5[1,4]
F5[4,1] <- 0.3961-0.2286*d5[4,1]
F5[1,5] <- 0.3961-0.2286*d5[1,5]
F5[5,1] <- 0.3961-0.2286*d5[5,1]
##############################################################
####Extending F for IS6

##############################################################
d6
F6<-matrix(0,nrow=5,ncol=5)
F6[2:5,2:5] <- F
F6[1,1] <- 0.3961-0.2286*d6[1,1]
F6[1,2] <- 0.3961-0.2286*d6[1,2]
F6[2,1] <- 0.3961-0.2286*d6[2,1]
F6[1,3] <- 0.3961-0.2286*d6[1,3]
F6[3,1] <- 0.3961-0.2286*d6[3,1]
F6[1,4] <-0.3961-0.2286*d6[1,4]
F6[4,1] <- 0.3961-0.2286*d6[4,1]
F6[1,5] <- 0.3961-0.2286*d6[1,5]
F6[5,1] <- 0.3961-0.2286*d6[5,1]
##############################################################












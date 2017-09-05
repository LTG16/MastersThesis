library(MASS)
library(mvtnorm)
library(openair)
library(imputeTS)
library(geosphere)##########

##### In this code we use the F matrices, which we obtained in "ExtendingFKF"

# We start with the kalman filter function to obtain S_T etc. 
kalmanOptimalDesignPred<-function(y,F,H,R,sigma_eta,d2,phi,sites,l){
  N <- dim(y)[1]
  n <- length(sites)
  mu0 <- colMeans(y)
  Q<-sigma_eta*exp(-d2*phi^2)
  Sigma0 <- Q
  
  m_n <- P_n <- K_n <- m_nn1 <- P_nn1 <- lik_m<- lik_s<-Y_n<-NULL
  loglikeli<-0
  m_nn1 <- matrix(rep(0,n),nrow=1)
  P_nn1 <- Q
  for(i in 1:N){
    
    m_n <- rbind(m_n,t(F%*%m_nn1[i,]))
    P_n <- rbind(P_n,F%*%P_nn1[1:(n)+(i-1)*n,1:n]%*%t(F)+Q)
    K_n <- rbind(K_n,P_n[1:(n)+(i-1)*n,1:n]%*%t(H)%*%
                   solve(H%*%P_n[1:(n)+(i-1)*n,1:n]%*%t(H)+R))
    lik_m <- rbind(lik_m,t(H%*%m_n[i,]))
    lik_s <- rbind(lik_s,H%*%P_n[1:(n)+(i-1)*n,1:n]%*%t(H)+R)
    m_nn1 <- rbind(m_nn1,m_n[i,]+t(K_n[1:(n)+(i-1)*n,1:n]%*%(y[i,]-H%*%m_n[i,])))
    P_nn1 <- rbind(P_nn1, (diag(1,n)- K_n[1:(n)+(i-1)*n,1:n]%*%H)%*%P_n[1:(n)+(i-1)*n,1:n])
  }
  ######### Now the prediction starts from N to l=10
  
  m_nn1 <- m_nn1[-1,]
  P_nn1 <- P_nn1[-c(1:n),]
  lik_spred<-list()
  for(j in N+1:l){
    
    y <- rbind(y,m_nn1[j-1,])
    m_n <- rbind(m_n,t(F%*%m_nn1[j-1,]))
    P_n <- rbind(P_n,F%*%P_nn1[1:(n)+(j-2)*n,1:n]%*%t(F)+Q)
    K_n <- rbind(K_n,P_n[1:(n)+(j-2)*n,1:n]%*%t(H)%*%
                   solve(H%*%P_n[1:(n)+(j-2)*n,1:n]%*%t(H)+R))
    lik_m <- rbind(lik_m,t(H%*%m_n[j-1,]))
    #lik_mpred <- t(H%*%m_n[j-1,])
    #lik_s <- rbind(lik_s,H%*%P_n[1:(n)+(j-1)*n,1:n]%*%t(H)+R)
    lik_spred[[j]]<-H%*%P_n[1:(n)+(j-1)*n,1:n]%*%t(H)+R
    m_nn1 <- rbind(m_nn1,m_n[j,]+t(K_n[1:(n)+(j-2)*n,1:n]%*%(y[j,]-H%*%m_n[j,])))
    P_nn1 <- rbind(P_nn1, (diag(1,n)- K_n[1:(n)+(j-1)*n,1:n]%*%H)%*%P_n[1:(n)+(j-1)*n,1:n])
  }
  
  #######
  ret <- list(lik_spred=lik_spred,lik_s=lik_s)
  return(ret)
}


##########we start by initialising the datasets
library(MASS)
library(mvtnorm)
##################################################################################
dataUB3<-importKCL(site =c("BL0","KC1", "LW1","WA2","WM0"), year = 2008:2009,
                   pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
#timePlot(dataUB,type="site")
dataUB3<-timeAverage(dataUB3, period="day",type=c("site","code","site.type"))
timePlot(dataUB3,type="site")
sites<-unique(dataUB3$code)



polution.data <- dataUB3[ which(dataUB3$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB3[ which(dataUB3$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data3 <- as.data.frame(polution.data)
data3<-matrix(nrow=233,ncol=5)

data3[,1] <- na.interpolation(pol.data3[351:583,1])
data3[,2] <- na.interpolation(pol.data3[351:583,2])
data3[,3] <- na.interpolation(pol.data3[351:583,3])
data3[,4] <- na.interpolation(pol.data3[351:583,4])
data3[,5] <- na.interpolation(pol.data3[351:583,5])
df3 <- unique(dataUB3[,c(2,6,7)])
df1 <- df3
d3<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000

#####################################################################################
##################################################################################
dataUB4<-importKCL(site =c("CT3","KC1", "LW1","WA2","WM0"), year = 2008:2009,
                   pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
#timePlot(dataUB,type="site")
dataUB4<-timeAverage(dataUB4, period="day",type=c("site","code","site.type"))
timePlot(dataUB4,type="site")
sites<-unique(dataUB4$code)



polution.data <- dataUB4[ which(dataUB2$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB4[ which(dataUB4$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data4 <- as.data.frame(polution.data)
data4<-matrix(nrow=233,ncol=5)

data4[,1] <- na.interpolation(pol.data4[351:583,1])
data4[,2] <- na.interpolation(pol.data4[351:583,2])
data4[,3] <- na.interpolation(pol.data4[351:583,3])
data4[,4] <- na.interpolation(pol.data4[351:583,4])
data4[,5] <- na.interpolation(pol.data4[351:583,5])
df3 <- unique(dataUB4[,c(2,6,7)])
df1 <- df3
d4<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000

#####################################################################################
##################################################################################
dataUB5<-importKCL(site =c("EA7","KC1", "LW1","WA2","WM0"), year = 2008:2009,
                   pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
#timePlot(dataUB,type="site")
dataUB5<-timeAverage(dataUB5, period="day",type=c("site","code","site.type"))
timePlot(dataUB5,type="site")
sites<-unique(dataUB5$code)



polution.data <- dataUB5[ which(dataUB5$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB5[ which(dataUB5$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data5 <- as.data.frame(polution.data)
data5<-matrix(nrow=233,ncol=5)

data5[,1] <- na.interpolation(pol.data5[351:583,1])
data5[,2] <- na.interpolation(pol.data5[351:583,2])
data5[,3] <- na.interpolation(pol.data5[351:583,3])
data5[,4] <- na.interpolation(pol.data5[351:583,4])
data5[,5] <- na.interpolation(pol.data5[351:583,5])

df3 <- unique(dataUB5[,c(2,6,7)])
df1 <- df3
d5<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000

#####################################################################################
##################################################################################
dataUB6<-importKCL(site =c("IS6","KC1","LW1","WA2","WM0"), year = 2008:2009,
                   pollutant = "nox", met = FALSE,units = "mass", extra = FALSE, meta = TRUE)
#timePlot(dataUB,type="site")
dataUB6<-timeAverage(dataUB6, period="day",type=c("site","code","site.type"))
timePlot(dataUB6,type="site")
sites<-unique(dataUB6$code)



polution.data <- dataUB6[ which(dataUB6$code==sites[1] ), 5]
for(i in 2:length(sites)){
  polution.data  <- cbind(polution.data,dataUB6[ which(dataUB6$code==sites[i] ), 5])
}
for(r in 1:length(sites)){
  colnames(polution.data)[r] <- sites[r]
}
pol.data6 <- as.data.frame(polution.data)
data6<-matrix(nrow=233,ncol=5)

data6[,1] <- na.interpolation(pol.data6[351:583,1])
data6[,2] <- na.interpolation(pol.data6[351:583,2])
data6[,3] <- na.interpolation(pol.data6[351:583,3])
data6[,4] <- na.interpolation(pol.data6[351:583,4])
data6[,5] <- na.interpolation(pol.data6[351:583,5])

df3 <- unique(dataUB6[,c(2,6,7)])
df1 <- df3
d6<-distm(df1[,c(3,2)],df1[,c(3,2)],fun = distHaversine)/10000

#####################################################################################
y3<-log(data3)

#### here we initialise the fixed values, for F we need to run the code "ExtendingFKF" before this one
H<-diag(length(sites))
phi3<-0.2621347
sigma_eta3<-0.28995
R3<-0.02582537*diag(length(sites))
#####################################################################################
y4<-log(data4)

#####################################################################################
y5<-log(data5)

#####################################################################################
y6<-log(data6)


#############################################################################################




#######################################################################################
# Optimal design for adding BL0

X3<-kalmanOptimalDesignPred(y=y3,F=F3,H=H,R=R3,sigma_eta=sigma_eta3,d2=d3,phi=phi3,sites,l=10)

l<-10
N<-dim(y3)[1]
n<-dim(y3)[2]
opt<-0
opt2<-0
for(j in N+1:l){
  opt2<- (-(1/2)*n)-((1/2)*log((2*pi)^n*det(X3$lik_spred[[j]])))
  opt<-opt+opt2
}
#############2. Optimal Design kriterium

det(X3$lik_spred[[243]])

#######################################################################################
#######################################################################################
# Optimal design for adding CT3

X4<-kalmanOptimalDesignPred(y=y4,F=F4,H=H,R=R3,sigma_eta=sigma_eta3,d2=d4,phi=phi3,sites,l=10)

l<-10
N<-dim(y4)[1]
n<-dim(y4)[2]
opt<-0
opt2<-0
for(j in N+1:l){
  opt2<- (-(1/2)*n)-((1/2)*log((2*pi)^n*det(X4$lik_spred[[j]])))
  opt<-opt+opt2
}

#############2. Optimal Design kriterium

det(X4$lik_spred[[243]])

#######################################################################################
#######################################################################################
# Optimal design for adding EA7

X5<-kalmanOptimalDesignPred(y=y5,F=F5,H=H,R=R3,sigma_eta=sigma_eta3,d2=d5,phi=phi3,sites,l=10)

l<-10
N<-dim(y3)[1]
n<-dim(y3)[2]
opt<-0
opt2<-0
for(j in N+1:l){
  opt2<- (-(1/2)*n)-((1/2)*log((2*pi)^n*det(X5$lik_spred[[j]])))
  opt<-opt+opt2
}
#############2. Optimal Design kriterium

det(X5$lik_spred[[243]])


#######################################################################################
#######################################################################################
# Optimal design for adding IS6

X6<-kalmanOptimalDesignPred(y=y6,F=F6,H=H,R=R3,sigma_eta=sigma_eta3,d2=d6,phi=phi3,sites,l=10)

l<-10
N<-dim(y3)[1]
n<-dim(y3)[2]
opt<-0
opt2<-0
for(j in N+1:l){
  opt2<- (-(1/2)*n)-((1/2)*log((2*pi)^n*det(X6$lik_spred[[j]])))
  opt<-opt+opt2
}
#############2. Optimal Design kriterium

det(X6$lik_spred[[243]])

#######################################################################################



plot(X$m_nn1[501:511,1],type="l")




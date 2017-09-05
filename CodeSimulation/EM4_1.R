##Necessary Functions

#######################################################################################
#EM Q spatial and F restricted for simulated Data
#######################################################################################

#start value calculation functions
calc_C0 <- function(y, mu){
  T <- dim(y)[1]
  C  <- matrix(0,dim(y)[2],dim(y)[2])
  for(t in 1:T){
    C <- C + (y[t,]-mu)%*%t(y[t,]-mu)
  }
  return(C/T)
}

calc_C1 <- function(y,mu){
  T <- dim(y)[1]
  C  <- matrix(0,dim(y)[2],dim(y)[2])
  for(t in 1:(T-1)){
    C <- C + (y[t+1,]-mu)%*%t(y[t,]-mu)
  }
  return(C/T)
}
#######################################################################################
#stopping criteria function

calc.differences <- function(Est.old,Est){
  F.diff     <- abs(Est.old[[1]] - Est[[1]])
  H.diff     <- abs(Est.old[[2]] - Est[[2]])
  R.diff     <- abs(Est.old[[3]] - Est[[3]])
  sigma_eta.diff <- abs(Est.old[[4]] - Est[[4]])
  phi.diff   <- abs(Est.old[[5]] - Est[[5]])
  
  eps <- 0.0001
  
  diff.matrix<- matrix(eps,ncol(Est[[1]]),nrow(Est[[1]]))
  
  if(all(F.diff <= diff.matrix) && all(H.diff <= diff.matrix) && all(R.diff <= diff.matrix) && 
     all(sigma_eta.diff <= diff.matrix) && all(phi.diff <= diff.matrix)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
#######################################################################################


s00 <- function(x,P){
  N<- dim(x)[1]
  n<- dim(x)[2]
  s <- matrix(0,n,n)
  for(t in 1:N){
    s <- s + as.numeric(x[t,])%*%t(as.numeric(x[t,])) + P[1:n+(t-1)*n,1:n]
  }
  return(s)
}

###############################################################################################################
s10 <- function(x,P){
  n<- dim(x)[2]
  N <- dim(x)[1]
  s <- matrix(0,n,n)
  for(t in 2:N){
    s <- s+as.numeric(x[t,])%*%t(as.numeric(x[t-1,])) + P[1:n+(t-2)*n,1:n]  # hier nochmal schauen
  }
  r <- s
  return(r)
}
###############################################################################################################
s11 <- function(x,P){
  n<- dim(x)[2]
  N <- dim(x)[1]
  s <- matrix(0,n,n)
  for(t in 2:N){
    s <- s+as.numeric(x[t,])%*%t(as.numeric(x[t,])) + P[1:n+(t-1)*n,1:n]
  }
  return(s)
}
####################################################################################################################
#Function to estimate R
sum4 <- function(x,y,H,P.new){
  r4 <- 0
  n<- dim(x)[2]
  N <- dim(x)[1]
  for(i in 1:N){
    
    a <- as.matrix(as.numeric((y[i,]-H%*%as.numeric(x[i,]))))
    b <- H%*%P.new[1:n+(i-1)*n,1:n]%*%t(H)
    #r4 <- r4 + sum(diag(a%*%t(a)+b))
    r4 <- r4 + sum(diag((a%*%t(a)+b)))
  }
  r4
}

###################   ESTIMATE THE PARAMETERS OF THE KALMAN FILTER ##################################################

estimate <- function(x,y,P.new,P.smoothed,phi,sigma_eta,d2){
  #estimate <- function(x,y,P.new,P.smoothed){
  
  N <- dim(x)[1]
  m <- dim(x)[2]
  Q<-sigma_eta*exp(-phi^2*d2)
  ################## CALCULATE HELPSUMS
  S00 <- s00(s$m_s_n,s$P_s_n)
  S10 <- s10(s$m_s_n,s$P_s_nn1)
  S11 <- s11(s$m_s_n,s$P_s_n)
  ################## ESTIMATE F,Q AND  H (which is fixed in this case)
  restrf <-function(x){
    -2*sum(diag(solve(Q)%*%(x*diag(m))%*%t(S10)))+sum(diag(solve(Q)%*%(x*diag(m))%*%S00%*%t(x*diag(m))))
  }
  F1<-optimize(restrf,interval=c(0,1),maximum=FALSE)$minimum
  F<-F1*diag(m)
  H   <- diag(1,m,m)
  
  ################## ESTIMATE SIGMA_ETA AND PHI
  
  A   <- S11-S10%*%t(F)-F%*%t(S10) + F%*%S00%*%t(F)
  parameters <- update.phi(phi,d2,A,sigma_eta,N=N,m=m)
  sigma_eta <-parameters[1]
  phi   <- parameters[2] 
  
  sigma_eps <- 1/(N*m)*sum4(x=s$m_s_n,y=y,H=H,P.new=s$P_s_n)
  R<-sigma_eps*diag(m)
  
  ret <- list(F,H,R,sigma_eta,phi)
  names(ret) <- c("F","H","R","sigma_eta","phi")
  ret
}

##################################################################
abl1 <-function(d2,phi) {
  -2*phi*d2*exp(-phi^2*d2)
}
abl2<- function(d2,phi){
  -2*d2*exp(-phi^2*d2)+4*phi^2*d2^2*exp(-phi^2*d2)
}
update.phi <- function(phi,d2,A,sigma_eta,N,m){
  
  # first and second derivative of the spatial matrix
  C         <- exp(-phi^2*d2)
  C.first   <-  abl1(d2,phi)
  C.second  <- abl2(d2,phi)
  C.solve <- solve(C)
  
  g.zero <- function(x){
    N*log(det(sigma_eta*exp(-x^2*d2)))+sum(diag(solve(sigma_eta*exp(-x^2*d2))%*%A))
  }
  g.sigma <- function(x){
    N*m*log(x)+N*log(det(C))+(1/x)*sum(diag(C.solve%*%A))
  }
  sigma_eta <- (1/(N*m))*sum(diag(C.solve%*%A))
  phi <- optimize(g.zero,interval=c(0,10),maximum=FALSE)$minimum
  
  
  r <- c(sigma_eta,phi)
  return(r) 
}



#########################################################################

#KF Function
kalman<-function(y,F,H,R,sigma_eta,d2,phi,sites){
  N <- dim(y)[1]
  n <- length(sites)
  mu0 <- colMeans(y,na.rm=T)
  Q<-sigma_eta*exp(-d2*phi^2)
  Sigma0 <- Q
  
  m_n <- P_n <- K_n <- m_nn1 <- P_nn1 <- NULL
  m_nn1 <- matrix(rep(0,n),nrow=1)
  P_nn1 <- Q
  for(i in 1:N){
    
    m_n <- rbind(m_n,t(F%*%m_nn1[i,]))
    P_n <- rbind(P_n,F%*%P_nn1[1:(n)+(i-1)*n,1:n]%*%t(F)+Q)
    K_n <- rbind(K_n,P_n[1:(n)+(i-1)*n,1:n]%*%t(H)%*%
                   solve(H%*%P_n[1:(n)+(i-1)*n,1:n]%*%t(H)+R))
    m_nn1 <- rbind(m_nn1,m_n[i,]+t(K_n[1:(n)+(i-1)*n,1:n]%*%(y[i,]-H%*%m_n[i,])))
    P_nn1 <- rbind(P_nn1, (diag(1,n)- K_n[1:(n)+(i-1)*n,1:n]%*%H)%*%P_n[1:(n)+(i-1)*n,1:n])
  }
  m_nn1 <- m_nn1[-1,]
  P_nn1 <- P_nn1[-c(1:n),]
  ret <- list(m_n=m_n,m_nn1=m_nn1,P_n=P_n,P_nn1=P_nn1,K_n=K_n)
  return(ret)
}

smoothing <- function(y,m_n,m_nn1,P_n,P_nn1,F,H,K_n){
  N <- dim(y)[1]
  n <- dim(y)[2]
  
  m_s_n <- matrix(m_nn1[N,],nrow=1)
  P_s_n <- P_nn1[1:(n)+(N-1)*n,1:n]
  P_s_nn1 <- (diag(1,n)-K_n[1:(n)+(N-1)*n,1:n]%*%H)%*%F%*%P_nn1[1:(n)+(N-2)*n,1:n]
  G_n <- NULL
  for(i in (N-1):1){
    G_n <- rbind(P_nn1[1:(n)+(i-1)*n,1:n]%*%t(F)%*%solve(P_n[1:(n)+(i)*n,1:n]),G_n)
    m_s_n <- rbind(m_nn1[i,]+t(G_n[1:n,1:n]%*%(m_s_n[1,]-m_n[i+1,])),m_s_n)
    P_s_n <- rbind(P_nn1[1:(n)+(i-1)*n,1:n]+G_n[1:n,1:n]%*%(P_s_n[1:n,1:n]-
                                                              P_n[1:(n)+(i-1)*n,1:n])%*%t(G_n[1:n,1:n]),P_s_n)
  }
  
  for(i in (N-2):1){
    P_s_nn1 <- rbind(P_nn1[1:(n)+(i)*n,1:n]%*%t(G_n[1:(n)+(i-1)*n,1:n])+
                       G_n[1:(n)+(i)*n,1:n]%*%(P_s_nn1[1:n,1:n]-F%*%P_nn1[1:(n)+(i)*n,1:n])
                     %*%t(G_n[1:(n)+(i)*n,1:n]),P_s_nn1)
  }
  
  r <- list(m_s_n=m_s_n,P_s_n=P_s_n,P_s_nn1=P_s_nn1)
  return(r)
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#initialise dataset
library(MASS)
set.seed(190493)
sites<-c(1,1,1,1)
coords<-matrix(runif(4*length(sites),0,1),ncol=4)
d2<-as.matrix(dist(coords))
phi<-0.5
sigma_eta<-0.8
Q<-sigma_eta*exp(-phi^2*d2)
H<-diag(length(sites))
F<-0.5*diag(4)
sigma2<-0.5
R<-sigma2*diag(length(sites))
y<-matrix(0,200,4)
x<-matrix(0,200,4)
x[1,]<-c(0,0,0,0)
for (i in 2:200){
  x[i,]<-F%*%x[i-1,]+mvrnorm(n=1,c(0,0,0,0),Q)
  y[i,]<-H%*%x[i,]+mvrnorm(n=1,c(0,0,0,0),R)
}

############# STARTING VALUES ###############################################################################

mu <- colMeans(y,na.rm=TRUE)

C_0 <- calc_C0(y,mu)
C_1 <- calc_C1(y,mu)

# Nugget effect
sigma_R <- 0.25

H1 = diag(length(sites))
F1 = C_1%*%solve(C_0)
Q1 = C_0 - C_1%*%solve(C_0)%*%C_1
R1 = sigma_R * diag(length(sites))
phi_start<-0.5 #select starting values
sigma_eta_start<-0.5 #select starting values

##################   INITIAL KALMAN FILTERING STEP WITH ESTIMATED VALUES   ##########################
###########################
P      <- list()
k      <- 1
l      <- list()




l <- kalman(y,F=F1,H=H1,R=R1,sigma_eta=sigma_eta_start,phi=phi_start,sites,d2=d2)


s <- smoothing(y=y,m_n=l$m_n,m_nn1=l$m_nn1,P_n=l$P_n,P_nn1=l$P_nn1,F=F,H=H,K_n=l$K_n)

x           <- s$m_s_n
P           <- l$P_nn1
P.smoothed  <- s$P_s_nn1
P.new       <- s$P_s_n
x.kalman    <- l$m_nn1




num_iter        <- 1
converged       <- FALSE

F.old           <- matrix(0,length(sites),length(sites))
H.old           <- matrix(0,length(sites),length(sites))
R.old           <- matrix(0,length(sites),length(sites))
sigma_eta.old       <- 0
phi.old         <- 0
Est.old         <- list(F.old,H.old,R.old,sigma_eta.old,phi.old)

#####################   EM - ALGORITHM       ##############################################

pb <- txtProgressBar(min = 0, max = 10000, style = 3)
phi<-0.5  #to start with
sigma<-0.5 #to start with



#Start model calibration

while(converged == FALSE && num_iter < 10000){
  
  setTxtProgressBar(pb, num_iter)
  
  
  ########
  Est <- estimate(x=s$m_s_n,y=y,P.new=s$P_s_n,P.smoothed=s$P_s_nn1,phi=phi,sigma_eta = sigma_eta,d2)
  
  ########
  
  l <- kalman(y,F=Est[[1]],H=Est[[2]],R=Est[[3]],sigma_eta=Est[[4]],phi=Est[[5]],sites=sites,d2=d2)
  s <- smoothing(y=y,m_n=l$m_n,m_nn1=l$m_nn1,P_n=l$P_n,P_nn1=l$P_nn1,F=Est[[1]],H=Est[[2]],K_n=l$K_n)
  
  P.new       <- s$P_s_n
  x.kalman    <- l$m_nn1
  x           <- s$m_s_n
  
 
  P           <- l$P_nn1
  P.smoothed  <- s$P_s_n
  
  diff <- calc.differences(Est.old,Est)
  if(diff == FALSE){
    Est.old[[1]]  <- Est[[1]]
    Est.old[[2]]  <- Est[[2]]
    Est.old[[3]]  <- Est[[3]]
    Est.old[[4]]  <- Est[[4]]
    Est.old[[5]]  <- Est[[5]]
  }
  
  phi <- Est[[5]]
  sigma_eta <- Est[[4]]
  
  
  if(diff == TRUE){
    converged <- TRUE
  }
  num_iter <- num_iter+1
}

close(pb) #closes progress bar       num_iter < 10000


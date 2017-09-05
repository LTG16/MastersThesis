###########################
##########################
### Gibbs sampler for parameters. Q is restricted and F, which is called M, is not

library(tmvtnorm)
library(MASS)
library(MCMCpack)
library(geosphere)
library(truncnorm)
######initialise dataset
set.seed(190493)
TMax<-200
sites<-c(1,1,1,1)
coords<-matrix(runif(2*length(sites),0,1),ncol=2)
d2<-as.matrix(dist(coords))
phi<-0.5
sigma_eta<-0.8
Q<-sigma_eta*exp(-d2*phi^2)
H<-diag(length(sites))
M<-0.5*diag(length(sites))
sigma2<-0.5
R<-sigma2*diag(length(sites))
z<-matrix(0,TMax+1,length(sites))
y<-matrix(0,TMax+1,length(sites))
y[1,]<-c(0,0,0,0)
for (i in 2:TMax+1){
  y[i,]<-M%*%y[i-1,]+mvrnorm(n=1,c(0,0,0,0),Q)
  z[i,]<-H%*%y[i,]+mvrnorm(n=1,c(0,0,0,0),R)
}


#############

### select hyperparameters
mu_0<-colMeans(z)
sigma2_0<-0.5
phi_0<-0.0001
Sigma_0<-sigma2_0*exp(-phi_0^2*d2)
a<-2
b<-1
vQ<-length(sites)
CQ<-cov(z)
mu_m<-0
Sigma_m<-1

phi_a<-2
phi_b<-1


#####Select initial values

Y_init <- t(matrix(colMeans(z,na.rm = TRUE),length(sites),dim(z)[1]-1))
M_init <- 0.1*diag(length(sites))
sigma2_init <- 0.9
sigma2eta_init <- 1.2
phi_init<-0.5

####Daten vorbereitung und initialise
num_iter        <- 2
num_total<-10000
gibbs_list <- list()
gibbs_list$Yt <-list()
length(gibbs_list$Yt) <- num_total
for( i in 1:num_total){
  gibbs_list$Yt[[i]]<-matrix(nrow=TMax+1,ncol=length(sites))
}
gibbs_list$sigma <- vector(length=num_total)
gibbs_list$sigma_eta <- vector(length=num_total)
gibbs_list$phi <- vector(length=num_total)
gibbs_list$Q <- list()
length(gibbs_list$Q)<-num_total
for( i in 1:num_total){
  gibbs_list$Q[[i]]<-matrix(nrow=length(sites),ncol=length(sites))
}
gibbs_list$M <- list()
length(gibbs_list$M)<-num_total
for( i in 1:num_total){
  gibbs_list$M[[i]]<-matrix(nrow=length(sites),ncol=length(sites))
}
##################
##################   Initialise


gibbs_list$Yt[[1]][1,] <- mvrnorm(n=1,mu_0,Sigma_0)

gibbs_list$Yt[[1]][-1,] <- Y_init

gibbs_list$sigma[1] <- sigma2_init

gibbs_list$Q[[1]] <- sigma2eta_init*exp(-phi_init^2*d2)

gibbs_list$phi[1] <-phi_init

gibbs_list$sigma_eta[1] <-sigma2eta_init

gibbs_list$M[[1]] <- M_init

###calculations before sampling




####While loop
pb <- txtProgressBar(min = 0, max = num_total, style = 3)
system.time(
  while( num_iter < num_total){
    
    setTxtProgressBar(pb, num_iter)
    
    #sample Y0
    Sigma_0 <- sigma2_0*exp(-gibbs_list$phi[num_iter-1]^2*d2)
    calc_Sigma <- solve(Sigma_0)
    calc_0 <- solve(Sigma_0)%*%mu_0
    
    temp_MQ <- t(gibbs_list$M[[num_iter-1]])%*%solve(gibbs_list$Q[[num_iter-1]])
    temp_Lambda<-solve(temp_MQ%*%gibbs_list$M[[num_iter-1]]+calc_Sigma)
    
    gibbs_list$Yt[[num_iter]][1,] <- mvrnorm(n=1,
                                             mu=(temp_Lambda%*%(temp_MQ%*%gibbs_list$Yt[[num_iter-1]][2,]+calc_0)),
                                             Sigma= temp_Lambda)
    #sample Yt
    temp_Q <- solve(gibbs_list$Q[[num_iter-1]])
    temp_Lambda<-solve(t(H)%*%H/gibbs_list$sigma[[num_iter-1]]+temp_MQ%*%gibbs_list$M[[num_iter-1]]+temp_Q)
    for (j in 2:TMax){
      gibbs_list$Yt[[num_iter]][j,] <- mvrnorm(n=1,
                                               mu=(temp_Lambda%*%(t(H)%*%z[j,]/gibbs_list$sigma[[num_iter-1]]+temp_MQ%*%
                                                                    gibbs_list$Yt[[num_iter-1]][j+1,]+t(temp_MQ)%*%
                                                                    gibbs_list$Yt[[num_iter]][j-1,])),
                                               Sigma=temp_Lambda)
    }
    
    #sample YT
    temp_Lambda <- solve(t(H)%*%H/gibbs_list$sigma[[num_iter-1]]+temp_Q)
    
    gibbs_list$Yt[[num_iter]][TMax+1,] <- mvrnorm(n=1,
                                                  mu=(temp_Lambda%*%(t(H)%*%z[j,]/gibbs_list$sigma[[num_iter-1]]
                                                                     +t(temp_MQ)%*%gibbs_list$Yt[[num_iter]][j-1,])),
                                                  Sigma=temp_Lambda)
    
    
    #sample R (sigma2)
    
    gibbs_list$sigma[num_iter]<-rinvgamma(n=1,shape=((TMax*length(sites))/2+a),scale=(b+0.5*sum((as.vector(z[-1])-
                                                                                                   as.vector(gibbs_list$Yt[[num_iter]][-1]))^2)))
    
    #sample Sigma_eta
    
    expsolve <- solve(exp(-gibbs_list$phi[num_iter-1]^2*d2))
    term1<-t(gibbs_list$Yt[[num_iter]][-1,])-(gibbs_list$M[[num_iter-1]]%*%t(gibbs_list$Yt[[num_iter]][-(TMax+1),]))
    gibbs_list$sigma_eta[num_iter]<-rinvgamma(n=1,shape=((TMax*length(sites))/2+a),scale=b+0.5*t(as.vector(term1))%*%(diag(TMax)%x%expsolve)%*%as.vector(term1))
    
    
    
    #sample Phi 
    
    
    ####log function for mcmc()
    
    
    phi_func<-function(phi_est){
      
      dgamma(phi_est,phi_a,scale=phi_b,log=TRUE)-TMax/2*log(det(exp(-phi_est^2*d2)))-1/(2*gibbs_list$sigma_eta[num_iter])*
        t(as.vector(term1))%*%(diag(TMax)%x%solve(exp(-phi_est^2*d2)))%*%as.vector(term1)-0.5*
        log(det(exp(-phi_est^2*d2)))-1/(2*sigma2_0)*t(gibbs_list$Yt[[num_iter]][1,]-mu_0)%*%
        solve(exp(-phi_est^2*d2))%*%(gibbs_list$Yt[[num_iter]][1,]-mu_0)
    }
    #metrophastings for phi
    gibbs_list$phi[num_iter]<-mean(MCMCmetrop1R(fun=phi_func,theta.init=gibbs_list$phi[num_iter-1],mcmc=1000,tune=1,thin=1))[1]
    
    
    #sample Q
    
    gibbs_list$Q[[num_iter]]<-gibbs_list$sigma_eta[num_iter]*exp(-gibbs_list$phi[num_iter]^2*d2)
    
    #sample M
    Q_snake<-diag(TMax)%x%solve(gibbs_list$Q[[num_iter]])
    
    Solvedelt <- t(as.vector(t(gibbs_list$Yt[[num_iter]][-(TMax+1),])))%*%Q_snake%*%
      as.vector(t(gibbs_list$Yt[[num_iter]][-(TMax+1),]))+1/Sigma_m
    
    temp_b <- t(as.vector(t(gibbs_list$Yt[[num_iter]][-(TMax+1),])))%*%Q_snake%*%
      as.vector(t(gibbs_list$Yt[[num_iter]][-1,]))+mu_m/Sigma_m
    
    
    
    gibbs_list$M[[num_iter]]<-rnorm(n=1,mean=solve(Solvedelt)*temp_b,sd=solve(Solvedelt))*diag(dim(z)[2])
    
    
    
    num_iter <- num_iter+1
  })

close(pb)


#############################

###Burnin of 2000

X<-apply(simplify2array(gibbs_list$M[2001:(num_iter-1)]),1:2,mean)
X

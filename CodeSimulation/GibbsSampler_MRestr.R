###########################
##########################
### Gibbs sampler for parameters. Q is not restricted so no spatial covariance function is used
# However F which is called M is restricted
library(MASS)
library(MCMCpack)
######initialise dataset
TMax<-200
sites<-c(1,1,1,1)
Q<-0.7*diag(4)+matrix(0.2,4,4)
H<-diag(length(sites))
M<-0.5*diag(4)
sigma2<-0.5
R<-sigma2*diag(length(sites))
z<-matrix(0,TMax+1,length(sites))
y<-matrix(0,TMax+1,length(sites))
y[1,]<-c(0,0,0,0)
for (i in 2:TMax+1){
  y[i,]<-M%*%y[i-1,]+mvrnorm(n=1,c(0,0,0,0),Q)
  z[i,]<-H%*%y[i,]+mvrnorm(n=1,c(0,0,0,0),R)
}


### select hyperparameters
mu_0<-colMeans(y)
Sigma_0<-diag(length(sites))
a<-2
b<-1
vQ<-length(sites)
CQ<-cov(y)
mu_m<-0
Sigma_m<-1

#####Select initial values

Y_init <- t(matrix(colMeans(y),4,TMax))
M_init <- 0.1*diag(length(sites))
sigma2_init <- 0.9
Q_init <- 0.9*diag(length(sites))

####Daten vorbereitung und initialise
num_iter        <- 2
#num_total <- num_burnin+1000
num_total<-10000
gibbs_list <- list()
gibbs_list$Yt <-list()
length(gibbs_list$Yt) <- num_total
for( i in 1:num_total){
  gibbs_list$Yt[[i]]<-matrix(nrow=TMax+1,ncol=length(sites))
}
gibbs_list$sigma <- vector(length=num_total)
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

gibbs_list$Q[[1]] <- Q_init

gibbs_list$M[[1]] <- M_init

###calculations before sampling

calc_Sigma <- solve(Sigma_0)
calc_0 <- solve(Sigma_0)%*%mu_0
solveSigma_m<-solve(Sigma_m)

####While loop
pb <- txtProgressBar(min = 0, max = 10000, style = 3)
system.time(
  while( num_iter < 10000){
    
    setTxtProgressBar(pb, num_iter)
    
    #sample Y0
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
    
    
    #sample R
    
    gibbs_list$sigma[num_iter]<-rinvgamma(n=1,shape=((TMax*length(sites))/2+a),scale=(b+0.5*sum((as.vector(z[-1])-
                                                                                                   as.vector(gibbs_list$Yt[[num_iter]][-1]))^2)))
    
    #sample Q
    term1<-t(gibbs_list$Yt[[num_iter]][-1,])-(gibbs_list$M[[num_iter-1]]%*%t(gibbs_list$Yt[[num_iter]][-(TMax+1),]))
    gibbs_list$Q[[num_iter]]<-riwish(v=vQ+TMax,S=(term1%*%t(term1)+vQ*CQ))
    
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



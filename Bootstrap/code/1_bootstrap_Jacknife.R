## Build a bias, standard deviation, and confidence interval estimator for the mean
## based on the bootstrap (use 10000 = nboot) and the jackknife.

Jackknife_confidence_interval <- function( v1, statfunc=sd, alpha = 0.05 ) {
  
  n1 <- length(v1)
  jackvec <- NULL
  mu0 <- statfunc(v1) 
  for(i in 1:n1)
  {
    mua<-statfunc(v1[-i])
    jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua)
    
  }
  jackbias<-mean(jackvec)-mu0 
  jacksd<-sd(jackvec)

  NLB<-mean(v1)-(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1) 
  NUB<-mean(v1)+(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1)
  
  list(mu0=mu0,jackestimate = mean(jackvec), jackbias=jackbias,jacksd=jacksd, normal.confidence.interval=c(NLB,NUB)) 
}

my.bootstrapci.ml<-function(vec0, nboot=10000, alpha=0.1){
  #extract sample size, mean and standard deviation from the original data
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))
  
  # create a vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  bootvec_estimate <- NULL
  jacckbiasvec<-NULL
  jacksdvec<-NULL
  
  #create the bootstrap distribution using a for loop
  for( i in 1:nboot){
    vecb<-sample(vec0,replace=T) #create mean and standard deviation to studentize
    meanb<-mean(vecb)
    sdb<-sqrt(var(vecb))
    #note since resampling full vector we can use n0 for sample size of vecb
    bootvec<-c(bootvec,(meanb-mean0)/(sdb/sqrt(n0))) 
    bootbiasvec<-c(bootbiasvec, meanb-mean0 ) 
    bootvec_estimate <- c(bootvec_estimate, meanb)
    
  }
  hist(bootvec)
  bootsd <- sd(bootvec)

  #Calculate lower and upper quantile of the bootstrap distribution
  bootbias<- mean(bootbiasvec)
  
  lq<-quantile(bootvec,alpha/2) 
  uq<-quantile(bootvec,1-alpha/2)
  
  #ADD the other two confidence intervals.
  #incorporate into the bootstrap confidence interval (what algebra supports this?) and output result 
  LB<-mean0-(sd0/sqrt(n0))*uq
  UB<-mean0-(sd0/sqrt(n0))*lq
 
  #since I have the mean and standard deviation calculate the normal confidence interval here as well 
  NLB<-mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n0))*qnorm(1-alpha/2)
  
  list( bootestimate = mean(bootvec_estimate), bootbias = bootbias, bootsd = bootsd, bootstrap.confidence.interval=c(LB,UB),normal.confidence.interval=c(NLB,NUB)) 
}

v1 = rnorm(500, mean = 10, sd = 2)
Jackknife_confidence_interval(v1,mean)
my.bootstrapci.ml(v1, nboot = 1000,alpha=0.05)

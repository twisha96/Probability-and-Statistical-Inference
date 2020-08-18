## Compare the coverage rates for the bootstrap confidence interval,
##the jackknife normal approximation confidence interval and the 
##central limit theorem based confidence interval. 
##For sample sizes 10, 30, and 100 alpha=0.05 (95% confidence)


library(tidyverse)
library(dplyr)

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
  
  #list(mu0=mu0,jackbias=jackbias,jacksd=jacksd, mean0 = mean0, jack_mean = jack_mean) 
  list(mu0=mu0,jackbias=jackbias,jacksd=jacksd, normal.confidence.interval=c(NLB,NUB)) 
  
}

my.bootstrapci.ml<-function(vec0, nboot=10000, alpha=0.1){
  #extract sample size, mean and standard deviation from the original data
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))
  
  # create a vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  
  #create the bootstrap distribution using a for loop
  for( i in 1:nboot){
    vecb<-sample(vec0,replace=T) #create mean and standard deviation to studentize
    meanb<-mean(vecb)
    sdb<-sqrt(var(vecb))
    #note since resampling full vector we can use n0 for sample size of vecb
    bootvec<-c(bootvec,(meanb-mean0)/(sdb/sqrt(n0))) 
    bootbiasvec<-c(bootbiasvec,meanb-mean0)
  }
 
  bootsd <- sd(bootvec)
  #Calculate lower and upper quantile of the bootstrap distribution
  bootbias<-mean(bootbiasvec)
  
  lq<-quantile(bootvec,alpha/2) 
  uq<-quantile(bootvec,1-alpha/2)
  
  #ADD the other two confidence intervals.
  #incorporate into the bootstrap confidence interval (what algebra supports this?) and output result 
  LB<-mean0-(sd0/sqrt(n0))*uq
  UB<-mean0-(sd0/sqrt(n0))*lq

  NLB<-mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n0))*qnorm(1-alpha/2)
  
  #since I have the mean and standard deviation calculate the normal confidence interval here as well 
  TLB <- mean0-(sd0/sqrt(n0))*qt(1-alpha/2, df = n0) 
  TUB<- mean0+(sd0/sqrt(n0))*qt(1-alpha/2, df = n0)
  
  list( bootbias = bootbias, bootsd = bootsd, bootstrap.confidence.interval=c(LB,UB),normal.confidence.interval=c(NLB,NUB), t.confidence.interval = c(TLB, TUB) ) 
}

Sim.func<-function(mu.val=3,n=30,nsim=1000)
{
  #create coverage indicator vectors for bootstrap and normal 
  cvec.boot<-NULL
  cvec.norm<-NULL
  cvec.t <- NULL
  
  #calculate real mean
  mulnorm<-(exp(mu.val+1/2)) #taking standard normal
  
  #run simulation
  
  for(i in 1:nsim){
    if((i/10)==floor(i/10)){ 
      print(i)
      #let me know computer hasnt died
    }
    
    #sample the simulation vector 
    vec.sample<-rlnorm(n,mu.val)
    
    #bootstrap it
    boot.list<-my.bootstrapci.ml(vec.sample) 
    boot.conf<-boot.list$bootstrap.confidence.interval 
    norm.conf<-boot.list$normal.confidence.interval 
    t.conf <- boot.list$t.confidence.interval 
    
    ## bootstap CI
    cvec.boot<-c(cvec.boot, (boot.conf[1]<mulnorm) * (boot.conf[2]>mulnorm) )  
    ## no need to call Jackknife since normal of both will be the same
    cvec.norm<-c(cvec.norm,(norm.conf[1]<mulnorm)*(norm.conf[2]>mulnorm)) 
    ##t
    cvec.t<-c(cvec.t,(t.conf[1]<mulnorm)*(t.conf[2]>mulnorm)) 
    
  }
  #calculate and output coverage probability estimates 
  list(boot.coverage=(sum(cvec.boot)/nsim), norm.coverage=(sum(cvec.norm)/nsim), t.coverage=(sum(cvec.t)/nsim)) 
}


## Question - 2 ##
Sim.func(n=30, nsim=1000)

## Question - 3 ##
df <-  data.frame(n_sample = numeric(), boot.cov = numeric(), norm.cov = numeric(), t.cov = numeric())
sample_size = c(10,30,100)
nsim = 1000
for (n_sample in sample_size)
{
  coverage.list <- Sim.func(mu.val=3,n=n_sample,nsim= nsim)
  boot.cov<-coverage.list$boot.coverage
  norm.cov<-coverage.list$norm.coverage
  t.cov <- coverage.list$t.coverage
  
  l <- list(n_sample = n_sample, boot.cov = boot.cov, norm.cov = norm.cov, t.cov=t.cov )
  df <- rbind(df, l)
  print(df)
}


df %>%
  gather(key = "variable", value = "value", -n_sample) %>%
  ggplot(aes(n_sample, value, col = variable)) +
  geom_line( position=position_dodge(0.5)) +
  geom_point(position=position_dodge(0.5), size=2) +
  scale_colour_manual(name="Confidence Interval", labels = c("Bootstrap","CLT based", "Jackknife Normal \n Approximation (t distibution)"), values =  c( "#ee5253", "#2e86de", "#10ac84")) +
  labs(y = "Coverage Rate",
       x = "Sample size", 
       subtitle = paste0('#simulations = ',nsim) ) + 
  theme_bw()





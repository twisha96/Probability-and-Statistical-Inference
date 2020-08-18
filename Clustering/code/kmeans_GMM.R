library(MASS)
library(ggplot2)
library(tidyverse)
library(clusterGeneration)  #for the multivariate function

parameter.random.initialization <- function(x, k, num.features){
  
  mu.init = matrix( nrow = k, ncol = num.features)
  
  for( i in c(1:k)){
    for (j in c(1:num.features)){
      mu.init[i,j] = runif(n = 1, max = max(x[,j]), min = min(x[,j]))
    }
  }
 
  # weight initialization
  wt.init = c(0.5,0.1,0.4)
  
  # cov initialization 
  cov.init = array(dim = c(num.features, num.features, k))
  for (i in (1:k)){
    cov.init[,,i] = genPositiveDefMat("eigen", dim = num.features)$Sigma
  }
  list(mu = mu.init, cov = cov.init, wt = wt.init)
}

e.step <- function(x,k,N,p.list){
  
  respon <- matrix(nrow = N, ncol = k)
  for (i in c(1:N)){
    for(j in c(1:k)){
      respon[i,j] = (p.list$wt[j]) * ( dmvn( as.numeric(x[i,]) , t(p.list$mu[j,]), p.list$cov[,,j]) )
    }
  }
  respon = respon/rowSums(respon)
  return(respon)
}

m.step <- function(x,respon, k, N, num.features){
  
  ## mu update ##  
  mu.updated = matrix(nrow = k, ncol = num.features)
  for(i in c(1:k)){
    mu.updated[i, ] = colSums(x*respon[,i])/sum(respon[,i])
  }
  
  ## wt update ##
  wt.updated = apply(respon,2,sum)/10

  ##cov update ##
  cov.updated = array(NA, dim =c(num.features,num.features,k))
  for(i in c(1:k)){
    diff = sweep(x, 2, c(mu.updated[i,1], mu.updated[i,2]), "-")
    cov.updated[,,i] = t(diff)  %*%  as.matrix(respon[,i] * diff )  / sum(respon[,i])
  }
  
  list(mu = mu.updated, cov = cov.updated,  wt = wt.updated)
  
}

# log likelihood 
log.likelihood <- function(x, k, N, p.list){
  inner.sum = vector()
  
  for (i in c(1:N)){
    temp = 0
    for(j in c(1:k)){
      temp = ( (p.list$wt[j]) * ( dmvn(as.numeric(x[i,]), t(p.list$mu[j,]), p.list$cov[,,j]) ) ) + temp
    }
    inner.sum = c(inner.sum, temp)
  }
  return(sum(log(inner.sum)))
}

plot.clusters <- function(x,respon){
  temp = respon
  temp = cbind(temp, apply(temp,1,which.max))
  temp1 = cbind(x, temp[,4]) %>% rename("cluster" = "temp[, 4]")
  print(temp1 %>%  
    ggplot(aes(temp1[,2],temp1[,4], color =as.factor(cluster))) + 
    geom_point() +
    scale_color_manual(
                      values=c("yellow", "red", "green"))+
    theme_bw() 
    )
}


#### Main function
em <- function(x, k=3){
  
  # no of features 
  num.features = ncol(x)
  # no. of data rows 
  N = nrow(x)
  
  log.likelihood.prev = 0
  Delta = 5
  # Initialization 
  p.list = parameter.random.initialization(x, k, num.features)  # parameter list
  n.steps = 0
  delta.list = c(1)
  
  while(Delta>0.1)
  {
   
    print(paste("----", n.steps, "----"))
    n.steps = n.steps + 1
    
    #E-STEP - recalculates responsibility matrix based on new cluster assignments.
    respon = e.step(x, k, N, p.list)
    
    #M-step - p.list contains the updated mean and variance of the clusters
    p.list = m.step(x,respon, k, N, num.features)
    
    ## log-likelihood ##
    log.likelihhod.new = log.likelihood(x,k,N, p.list)
    
    Delta = abs(log.likelihhod.new - log.likelihood.prev)
    delta.list = c(delta.list, Delta)
    log.likelihood.prev = log.likelihhod.new
    print(paste("Delta", Delta))
    
  }

  final.cluster = respon
  final.cluster = cbind(final.cluster, apply(respon,1,which.max))
  final.cluster = cbind(x, final.cluster[,4]) %>% rename("cluster" = "final.cluster[, 4]")
  
  with(x, pairs(x, col=final.cluster$cluster) )
  return(list(data_with_clusters = final.cluster, respon = respon, delta.list = delta.list, clusters.parameters = p.list))

}


# data
df <- read.csv("kmeans.csv")
#no of clusters 
k = 3
# running em algo
out.list = em(df, k) 
# extracting output
respon = out.list$respon
# data_with_clusters = out.list$data_with_clusters


## colors according to respon values ##
# rbPal <- colorRampPalette(c('red','green', 'blue'))
# data_with_clusters$Col <- rbPal(5)[as.numeric(cut(respon[,1],breaks = 5))]
# plot(data_with_clusters$b,data_with_clusters$d,col = data_with_clusters$Col, pch = 19)


# plotting Delta
t <- out.list$delta.list
df_temp <- data.frame(col_name = unlist(t[3:length(t)]))
df_temp %>% ggplot(aes(c(3:length(t)),col_name)) + geom_line() + xlab("n_steps") + ylab("Delta")  +theme_bw()






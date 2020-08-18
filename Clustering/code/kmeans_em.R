library(MASS)
library(ggplot2)
library(tidyverse)
library(matlib) #for Ginv(matrix) function

# assigning each data point to a cluster randomly at first 
clusters.init <- function(data, k){
  
  n = nrow(data)
  ## assign each point to random clusters
  clusters = rep(1:k, length.out = n, replace = TRUE)
  data = cbind(data, clusters)
  
  return(data)
}

# recalculating the centroid/means and the covariance of all the clusters, based on the points assigned to them
clusters.parameters <- function(data, k){
 
  num.features = ncol(data)-1
  #taking mean feature wise
  mean = aggregate(data, by = list(data$clusters), mean, na.rm = TRUE)
  cov.list <- lapply( sort(unique(data$clusters)), function(x) cov(data[data$clusters==x,-(num.features+1)],use="na.or.complete"))
 
  return(list(mean = mean[-1], cov.list = cov.list ))
}


#trial approach
assign.cluster <- function(data, clusters.list){
  
  n = nrow(data)
  num.features = ncol(data)-1
  
  for(i in c(1:n)){
    min_distance = +Inf
    for(j in  c(1:k) ) {
     
      cov = matrix(unlist(clusters.list$cov.list[j]), nrow = num.features, byrow = T)
      mean = clusters.list$mean[j,-(num.features+1)]
      diff = data[i,-ncol(data)] - mean
      temp  = as.matrix(diff) %*% (cov)
      distance = temp %*% as.matrix(t(diff))
       
      if(distance < min_distance){
          min_distance = distance
          data[i,]$clusters = j
        }
    }
  }
  
  return(data)
  
}

euc.dist <- function(x1, x2){
  return(sqrt(sum((x1 - x2) ^ 2)))
}
  

#re-assign cluster to each data point based on the points euclidean distance from the cluster centroids
assign.cluster.eucliean <- function(data, clusters.list){
  
  n = nrow(data)
  num.features = ncol(data)-1
  
  for(i in c(1:n)){
    min_distance = +Inf
    for( j in  unique(data$clusters) ){
        mean = clusters.list$mean[j,-(num.features+1)]
        distance = euc.dist(data[i,-ncol(data)], mean )
        if(distance < min_distance){
          min_distance = distance
          #assign the cluster with min distance from the data point
          data[i,]$clusters = j
         }
      }
    }
  return(data)
}

## 2D plot any two fearutres of the data
plot.clusters <- function( data){
  print(data %>%  
          ggplot(aes(data[,1],data[,2], color = as.factor(clusters))) + 
          geom_point() +
          scale_color_manual(
            values=c("red", "blue", "green")))
}

#############################################################
##functions to calculate DUNN INDEX (for convergence) #######
#############################################################
intracluster_dist <- function(cluster_mean, cluster_data){
  n = nrow(cluster_data)
  max = 0
  for (i in c(1:n)){
    distance = euc.dist(cluster_mean, cluster_data[i,])
    if(distance>max){
      max = distance
    }
  }
  return(max)
}

get_max_intracluster_distance <- function(data, mean){
  num.features = ncol(data)-1
  for (i in c(1:k)){
    t = lapply(sort(unique(data$clusters)), function(x) intracluster_dist(mean[x,], data[data$clusters==x, -(num.features+1)]))
  }
  
  return(max(unlist(t)))
}

get_min_intercluster_distance <- function(mean, k){
  
  
  min = +Inf
  for (i in c(1:(k-1))){
    for (j in c((i+1):k)){
      distance = euc.dist(mean[i,], mean[j,])
      #print(distance)
      if(distance<min){
        min = distance
      }
    }
  }
  return(min)
}

dunn_index <- function(data, mean, k){
  num.features = ncol(data)
  
  min_intercluster_distance = get_min_intercluster_distance(mean[-num.features], k)
  #print(min_intercluster_distance)
  max_intracluster_distance = get_max_intracluster_distance(data, mean[-num.features])
  #print(max_intracluster_distance)
  return(min_intercluster_distance/ max_intracluster_distance)
}

############################################
############## K MEANS #####################
############################################

k.means <- function(data, k){
  data = clusters.init(data, k)
  data.init.clusters = data
  plot.clusters(data)
  num.features = ncol(data_with_clusters)

  ##convergence of dunn index 
  prev_dunn_index = 0
  delta = 1
  n_steps = 1
  
  # stop when dunn index becomes constant
  while(delta!=0){
    print(paste("----", n_steps, "----"))
    #clusters.list contains the cluster parameters of mean and covariance for each cluster 
    clusters.list = clusters.parameters(data, k)
    # new column of 'clusters' gets appended 
    data = assign.cluster.eucliean(data,clusters.list )
    plot.clusters(data)
    
    new_dunn_index = dunn_index(data, clusters.list$mean[-(num.features)], k)
    print(new_dunn_index)
    delta = abs(new_dunn_index - prev_dunn_index)
    prev_dunn_index = new_dunn_index
    n_steps = n_steps + 1
  }
  
  return(list(data = data, clusters.list = clusters.list))
}
  
## Call the functions here 
df <- read.csv("kmeans.csv")
cols = c(1:5)
#rows = 1:100
rows = 1:nrow(df)
data = df[rows,cols]

#no of clusters 
k = 3

out.list = k.means(data, k)

data_with_clusters = out.list$data
clusters.parameters = out.list$clusters.list

print(clusters.parameters)

# Visualization #
num.features = ncol(data_with_clusters)
with( data_with_clusters[-(num.features)], pairs(data_with_clusters[-(num.features)], col=data_with_clusters$clusters) )

# Density of each feature #
print(data_with_clusters %>% 
  group_by(clusters) %>% 
  ggplot(aes(x = c,  fill = as.factor(clusters))) +
  geom_density( alpha = 0.35) +
  theme_bw())




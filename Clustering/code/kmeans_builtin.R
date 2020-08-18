mydata <- read.csv("kmeans.csv")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 3) # 5 cluster solution

# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)

# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         , lines=0)

with(mydata, pairs(mydata[-6], col=c(1:5)[fit$cluster])) 

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
#plotcluster(mydata, fit$cluster)

#fviz_cluster(final, data = df)
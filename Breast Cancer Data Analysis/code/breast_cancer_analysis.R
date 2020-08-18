# Analysis of breast cancer dataset #

# # for function ggcorr
if (!require("GGally")) install.packages("GGally")
if (!require("heatmaply")) install.packages("heatmaply")
if (!require("devtools"))  install.packages("devtools")
if (!require("tree"))  install.packages("tree")
if (!require("randomForest"))  install.packages("randomForest")
if (!require("naivebayes"))  install.packages("naivebayes")
if (!require("e1071"))  install.packages("e1071")

library(devtools)
install_github("vqv/ggbiplot")

library(ggplot2)
library(ggbiplot)
library(GGally)
library(dplyr)
# library(tidyverse)
library(heatmaply)
library(stats)
library(ggbiplot)
library(tree)
library(randomForest)
library(naivebayes)
library(e1071)


Data <- function(df){
  print("First few rows of the data")
  print(head(df))
  cat("Number of Benign(B) and Malingnant(M)")
  print(count(df, as.character(df$diagnosis)))
}

correlation <- function(df){
  h <- heatmaply_cor(
    cor(df),
    scale = "none",
    # reforder = F
    Rowv = F,
    Colv = F,
    revC = F
    
  )
  print(h)
}


distribution <- function(df){
  p = ggpairs(df, aes(color=diagnosis, alpha=0.75), lower=list(continuous="smooth"))+
    theme_bw()+
    labs(title="Cancer Mean") +
    theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=6))
  print(p)
}



pca <- function(df, df_class){
  df.pca <- prcomp(df, center = TRUE, scale. = TRUE)
  print(summary(df.pca))
  p = ggbiplot(df.pca, ellipse=TRUE, groups= as.character(df_class), choices = c(1,2), alpha = 0.5)
  print(p)
}

######################################################
################## Classification ##################
#####################################################


################## Logistic regression ##################
logistic_regression <- function(df, df_class){
 
  # changing B amd Ms to 0s and 1s
  df_class[df_class=="B"] <- 0
  df_class[df_class =="M"] <- 1
  df$diagnosis = as.numeric(t(df_class))
  train <- df[1:512,]
  test <- df[513:569,]
  
  model <- glm(diagnosis ~.,family=binomial(link='logit'), data=train)
  print(summary(model))
  
  # prediction
  fitted.results <- predict(model,newdata=subset(test, ,select=c(2:length(test))),type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  t <- table(fitted.results, test$diagnosis)
  cat("\nConfusion matrix for test prediction")
  print(t)
  misClasificError <- mean(fitted.results != test$diagnosis)
  print(paste('Accuracy',1-misClasificError))
  
}

################## Decision Tree ##################

decision_tree <- function(df, df_class){
  # need to chance B amd Ms to 0s and 1s
  df_class[df_class=="B"] <- 0
  df_class[df_class =="M"] <- 1
  df$diagnosis = as.numeric(t(df_class))
  train <- df[1:512,]
  test <- df[513:569,]
  
  dum0tree <- tree(diagnosis~radius_mean+texture_mean+smoothness_mean+compactness_mean+concavity_mean+concave.points_mean+symmetry_mean,data=train)
  pred2tree<-predict(dum0tree,newdata=test)
  pred2tree <- ifelse(pred2tree>0.5,1,0)
  t = table(pred2tree, test$diagnosis)
  cat("\nConfusion matrix for test prediction")
  print(t)
  misClasificError <- mean(pred2tree != test$diagnosis)
  print(paste('Accuracy',1-misClasificError))
  
  ## visualising tree ##
  par(mfrow=c(1,1))
  plot(dum0tree)
  text(dum0tree)
  
}


################## Random Forest ##################
random_forest <- function(df){
  train <- df[1:512,]
  test <- df[513:569,]
  model1 <- randomForest(train$diagnosis ~ ., data = train, importance = TRUE)
  model1
  predTest <- predict(model1, test, type = "class")
  # Checking classification accuracy
  t = table(predTest, test$diagnosis)
  cat("\nConfusion matrix for test prediction")
  print(t)
  missclassified <- sum(predTest != test$diagnosis)
  Accuracy <- 1 - (missclassified/nrow(test))
  print(paste('Accuracy', Accuracy))
}


################## Navie Bayes ##################
naive_bayes <- function(df){
  train <- df[1:512,]
  test <- df[513:569,]
  
  model_nb <- naiveBayes(train$diagnosis ~ ., data = train)
  predTest <- predict(model_nb, test, type = "class")
  t = table(predTest, test$diagnosis)
  cat("\nConfusion matrix for test prediction")
  print(t)
  missclassified <- sum(predTest != test$diagnosis)
  Accuracy <- 1 - (missclassified/nrow(test))
  print(paste('Accuracy', Accuracy))
  
}



wrapper <- function(input){
  
  rdf <- read.csv("data.csv")
  # removing last column of NAs
  df <- rdf[-ncol(rdf)]
  df_diagnosis = df[2:length(df)]
  df_mean = df[3:12]
  df_mean_diagnosis = df[2:12]
  df_class <- data.frame(lapply(df$diagnosis, as.character), stringsAsFactors=FALSE)
  
  if(input == "1")
    Data(df)
  
  else if(input == "2"){
    # of all 30 features
    correlation(df[-(1:2)])
    # of the means of the features
    correlation(df_mean)
  }
  
  else if(input == "3"){
    distribution(df_mean_diagnosis)
  }
  
  else if(input =="4"){
    cat("\n ********* PCA on all features *********\n")
    pca(df[c(-1,-2)], df_class)
    cat("\n ********* PCA on the means of the features *********\n")
    pca(df_mean, df_class)
  }
    
  else if(input == "5"){
    cat("\n  ********* Logstic Regression *********")
    # logistic_regression(df[c(-1)], df_class)
    logistic_regression(df_mean_diagnosis, df_class)
  }
  
  else if(input == "6"){
    cat("\n  ********* Naive Bayes *********")
    naive_bayes(df[c(-1)])
  }
  
  else if(input == "7"){
    cat("\n  ********* Decision Tree *********")
    decision_tree(df_mean_diagnosis, df_class)
  }
  
  else if(input == "8"){
    cat("\n  ********* Random Forest *********")
    random_forest(df_mean_diagnosis)
  }
}



repeat{
  
  cat("\n\n ###########################################################################
1. Data \n2. Data Analysis:Correlation \n3. Data Analysis: Distribution \n4. Data Analysis:PCA \n5. Classification:Logistic Regression \n6. Classification:Naive Bayes \n7. Classification:Decision Tree  \n8. Classification:Random Forest 
###########################################################################")
  input <- readline(prompt="Enter the number to see the results (write 'quit' to exit): ")
  wrapper(input)
  
  if(input =="quit")
    break;
}
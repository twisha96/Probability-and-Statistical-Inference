library(MASS)
library(tidyverse)

rdf <- read.csv("category5.csv")
rdf <- rdf[ -c(1),]

rdf <- rdf %>%
  rename("frequency" = "Frequency",
         "decade" = "Decade")

### DATA - Myqqplot ###
probability_vector <- rdf %>%
  group_by(frequency) %>%
  summarise(count = n()) %>%
  mutate(cum_sum = cumsum(count)) %>%
  mutate( probability = cum_sum/nrow(rdf)) 

rdf_probability <- full_join(rdf, probability_vector) %>%
                   select(-c(count))
### calculate q plot ###
lambda <- mean(rdf[["frequency"]])
#prob <- pull(probability_vector, probability) 
prob <- pull(rdf_probability, probability) 
q_prob <- qpois( prob, lambda)
rdf_probability["qprob"] = q_prob

#### Q-Q plot #####

plot( c(0,sort(freq)),c(0, q_prob) , xlim = c(0,15), ylim = c(0,15), xlab="data", ylab="theoretical quantile", main = "Poisson distribution for Category 5 hurricanes") 
lines ( c(0,sort(freq)), c(0,sort(freq)))




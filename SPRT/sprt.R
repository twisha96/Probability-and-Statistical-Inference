# Test for a Bernoulli random variable p<=.45 vs p>=.55 using an SPRT ##

#for ggplot
library(tidyverse)

#####################################
#### STRP Function for Bernoulli ####
#####################################

# alpha0 is Type1 error 
# alpha1 is Type2 error 
# p<=p1: NULL Hypothesis
# p>=p1: Alternate Hypothesis

strp_bernoulli <- function(alpha0 = 0.01, alpha1 = 0.01, p1 = 0.45, p2 = 0.55, bern_p = 0.3){

  S = 0
  log_likelihood = 0
  ## to keep track of number of steps required for convergance ##
  n_converge = 0  
  
  # calculating threshhold for stopping #
  A = log(alpha1/(1 - alpha0))   
  B = log((1 - alpha1)/alpha0)
  hypo_accepted = -1
  
  while(TRUE){
    
    n_converge = n_converge + 1
    
    # generating bernoulli RV with p = bern_p 
    data_point = rbinom(1, 1, bern_p)   
    
    # log-likelihood ratio
    log_likelihood = (data_point*p2 + (1-data_point)*(1-p2)) - (data_point*p1 + (1-data_point)*(1-p1))   
    
    # cumulative sum of the log-likelihood ratio
    S = S + log_likelihood          
    
    # Stopping Rule #
    if(S>=B){
      #Accept H1
      hypo_accepted = 1              
      break
    }
    if(S<=A){
      #Accept H0
      hypo_accepted = 0             
      break
    }
    #print(paste("data_point", data_point, "log_likelihood", log_likelihood, "S", S))
  }
  
  return(list(n_converge = n_converge, hypo_accepted = hypo_accepted))
}

#####################################
####    Simulation function      ####
#####################################

simulate_strp_bernoulli <- function(bern_p = 0.3, nsim = 10){      
  
  sum_n = 0
  H0_count = 0
  H1_count = 0
  
  ## Averaging ovet the STPR funtion   
  for(i in c(1:nsim)){            
    # calling the SRTP funtion which generates bernoulli variables with p = bern_p #
    strp_result = strp_bernoulli(bern_p=bern_p)    
    
    ## if H0 is accepted
    if(strp_result$hypo_accepted == 0){                  
      H0_count = H0_count + 1
    }
    ## if H1 is true
    if(strp_result$hypo_accepted == 1){                   
      H1_count = H1_count + 1
    }
    sum_n = sum_n + strp_result$n_converge
  }
  avg_steps = sum_n/nsim
  # print(paste("bern_p=", bern_p, ", avg_n=", avg_steps))
  return(list(avg_steps=avg_steps, H0_count=H0_count, H1_count=H1_count))
}



######### Testing for particular values ###########
simulate_strp_bernoulli(0.5, 100)



###### Running the simulation function for a range of p values ####
bern_p_list = seq(0, 1, by=0.05)

# data frame to store the result 
df<-data.frame(bern_p=numeric(), avg_steps_to_converge=numeric(), H0_count=numeric(), H1_count=numeric(), final_hypothesis_accepted=character())

#iterating over bern_p_list
for(p in bern_p_list){
  print(p)
  sim_result = simulate_strp_bernoulli(bern_p = p)

  if(sim_result$H0_count < sim_result$H1_count)
    final_H = "H1"
  if(sim_result$H0_count > sim_result$H1_count ){
    final_H = "H0"
  }
  if(sim_result$H1_count == sim_result$H0_count){
    final_H = "H0/H1"
  }
  
  new_data = data.frame(bern_p=p, avg_steps_to_converge=sim_result$avg_steps, H0_count=sim_result$H0_count, H1_count=sim_result$H1_count,  final_hypothesis_accepted=final_H, stringsAsFactors = FALSE)
  df<-rbind(df, new_data)
}

View(df)

## plotting the results ##
df %>% 
  ggplot(aes(bern_p, avg_steps_to_converge) ) +
  geom_line(aes(color = final_hypothesis_accepted), alpha = 0.5)  +
  geom_point(aes(color = final_hypothesis_accepted)) +
  #xlim(0,1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45)) + 
  scale_x_continuous(limits = c(0,1), breaks = bern_p_list) +
  labs(y = "#Steps for convergence",
        x = "bernoulli p") +
  scale_colour_manual(name="Hypothesis Accepted", labels = c("H0", "H0/H1", "H1"), values =  c( "#ee5253",  "#10ac84", "#2e86de")) 





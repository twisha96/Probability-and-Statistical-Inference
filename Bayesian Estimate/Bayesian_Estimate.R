# Baysian Estimate for distributions using conjugate priors

library(actuar) #for rpareto function

bayes_estimate_wrapper <- function(distribution, n = 10000){
  
  cat("\n")
  print(paste("-------", distribution, "-------"))
  
  if(distribution == "Binomial"){
    p = 0.4
    sample <- rbinom(n = 1000, size = 1, p)
   
    # Conjugate prior - Beta(alpha, beta)
    prior_alpha = 1
    prior_beta = 1
    
    # Posterior - Beta(alpha + n*mean(X), beta + n - n*mean(X))
    posterior_alpha = prior_alpha + sum(sample)
    posterior_beta = prior_beta + length(sample) - sum(sample)
    
    # Bayes estimate of p - Mean of posterior distribution
    estimate = posterior_alpha/ (posterior_alpha + posterior_beta)
    
    print("Paramters of posterior beta distribution: ")
    print(c(posterior_alpha, posterior_beta))
    print(paste("p = ", p))
    print(paste("Bayes estimate = ", estimate))
    
    # Dentisy 
    posterior_sample <- rbeta(n, posterior_alpha, posterior_beta)
    plot(density(posterior_sample))
    
  }
  
  else if(distribution == "Poisson"){
    lambda = 5
    sample <- rpois(n = 1000, lambda = lambda )
    
    # Conjugate prior - Gamma(alpha, beta)
    prior_alpha = 1
    prior_beta = 1
    
    # Posterior - Gamma(alpha + n*mean(X), beta + n)
    posterior_alpha = prior_alpha + sum(sample)
    posterior_beta = prior_beta + length(sample) 
    
    # Bayes estimate - Mean of posterior distribution
    estimate = posterior_alpha/ posterior_beta

    print("Paramters of posterior gamma distribution: ")
    print(c(posterior_alpha, posterior_beta))
    print(paste("lambda = ", lambda))
    print(paste("Bayes estimate = ", estimate))
    
    # Dentisy 
    posterior_sample <- rgamma(n, posterior_alpha, posterior_beta)
    plot(density(posterior_sample))
    
  }
  
  else if(distribution == "Exponential"){
    rate = 5
    sample <- rexp(n = 1000, rate = rate )
    
    # Conjugate prior - Gamma(alpha, beta)
    prior_alpha = 1
    prior_beta = 1
    
    # Posterior -  Gamma(alpha + n, beta + n*mean(X))
    posterior_alpha = prior_alpha + length(sample) 
    posterior_beta = prior_beta + sum(sample) 
  
    # Bayes estimate - Mean of posterior distribution
    estimate = posterior_alpha/ posterior_beta
    
    print("Paramters of posterior gamma distribution: ")
    print(c(posterior_alpha, posterior_beta))
    print(paste("rate = ", rate))
    print(paste("Bayes estimate = ", estimate))
    
    # Dentisy 
    posterior_sample <- rgamma(n, posterior_alpha, posterior_beta)
    plot(density(posterior_sample))
    
  }
  
  else if(distribution == "Geometric"){
    p = 0.7
    sample <- rgeom(n = 1000, p )
    
    # Conjugate prior - Beta(alpha, beta) 
    prior_alpha = 1
    prior_beta = 1
    
    # Posterior - Beta(alpha + n, beta + n*mean(X))
    posterior_alpha = prior_alpha + length(sample) 
    posterior_beta = prior_beta + sum(sample) 
    
    # Bayes estimate - Mean of posterior distribution
    estimate = posterior_alpha/ (posterior_alpha + posterior_beta)
    
    print("Paramters of posterior beta distribution: ")
    print(c(posterior_alpha, posterior_beta))
    print(paste("p = ", p))
    print(paste("Bayes estimate = ", estimate))
    
    # Dentisy 
    posterior_sample <- rbeta(n, posterior_alpha, posterior_beta)
    plot(density(posterior_sample))
    
  }
  
  
  else if(distribution == "Uniform"){
   
    sample <- runif(n = 1000, min = 0, max = 10)
    
    # Conjugate prior - Pareto
    prior_v0 = 1
    prior_k = 1
    
    # Posterior
    posterior_v0 = max(c(prior_v0, sample))
    posterior_k = prior_k + length(sample) 
    
  
    print("Paramters of posterior pareto distribution: ")
    print(c(posterior_v0, posterior_k))
  
    # Dentisy 
    posterior_sample <- rpareto(n, posterior_v0, posterior_k)
    plot(density(posterior_sample))
  }
  
  else if(distribution == "Normal"){
    sample <- rnorm(n = 1000, mean = 10, sd = 20)
    
    # Assuming alpha and beta for the prior distribution to be 1
    r <- 1
    tau <- 5
    mu <- 4
    prior_alpha <- 1
    prior_beta <- 2
    
    # Getting the posterior distribution parameters
    M_conditional_distribution_mu <- (tau*mu + length(sample)*mean(sample))/(tau + length(sample))
    M_conditional_distribution_precision <- (tau + length(sample))*r
    print("The parameters of the conditional posterior normal distribution of M when R=r is:")
    print(c(M_conditional_distribution_mu, M_conditional_distribution_precision))
    
    R_marginal_distribution_alpha <- prior_alpha + length(sample)/2
    R_marginal_distribution_beta <- prior_beta + 1/2*(sum((sample - mean(sample))**2)) + tau*length(sample)*((mean(sample) - mu)**2)/2*(tau + length(sample))
    print("The parameters of the marginal posterior gamma distribution of R is:")
    print(c(R_marginal_distribution_alpha, R_marginal_distribution_beta))
    
    # Generate the distibutions
    conditional_joint_distribution_of_M <- rnorm(n, mean = M_conditional_distribution_mu, 1/sqrt(M_conditional_distribution_precision))
    marginal_joint_distribution_of_R <- rgamma(n, R_marginal_distribution_alpha, R_marginal_distribution_beta)
    
    plot(density(conditional_joint_distribution_of_M))
    plot(density(marginal_joint_distribution_of_R))
  }
  
  
}

#bayes_estimate_wrapper("Binomial")


################################################################
#### Main wrapper function to run bayesian estimate for all distributions ####
################################################################
main <- function(){
  
  distributions <- c("Binomial", "Geometric", "Poisson", "Uniform", "Normal", "Exponential")
  
  n = 100
  for (distribution in distributions){
    bayes_estimate_wrapper(distribution)
  }
  
}


main()



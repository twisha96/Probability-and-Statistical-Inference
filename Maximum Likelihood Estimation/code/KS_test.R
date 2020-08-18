# Parametric Bootstrap using KS Test

ks_test_func <- function(distribution, nboot = 1000, input_data,theta_hat){
  
  mle_name = get(paste( distribution,"_mle", sep = ""))
  n <- length(input_data)
  
  if(distribution == "poisson"){
    # Initializing the D0 value
    q_hat <- qpois(c(1:n)/(n+1),theta_hat)
    D0 <- ks.test(input_data, q_hat, exact = NULL)$statistic
    Dvec<-NULL
    
    # Bootsrapping nboot times
    for(i in 1:nboot){
      x_star <- rpois(n, theta_hat)
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qpois(c(1:n)/(n+1), theta_hat_star)
      D_star <- ks.test(x_star, q_hat_star, exact = NULL)$statistic
      Dvec <- c(Dvec, D_star)
    }
    # Calculating the p-value considering all the ks_test results above D0
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "normal"){
    q_hat <- qnorm(c(1:n)/(n+1),mean = theta_hat[1], sd = theta_hat[2])
    
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rnorm(n,mean = theta_hat[1], sd =theta_hat[2])
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qnorm(c(1:n)/(n+1),mean = theta_hat_star[1], sd =theta_hat_star[2])
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "uniform"){
    q_hat <- qunif(c(1:n)/(n+1), theta_hat[1], theta_hat[2])
    
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- runif(n, theta_hat[1], theta_hat[2])
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qunif(c(1:n)/(n+1), theta_hat_star[1], theta_hat_star[2])
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "gamma"){
    q_hat <- qgamma(c(1:n)/(n+1), shape = theta_hat[1],  scale = theta_hat[2])
    
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rgamma(n, shape = theta_hat[1], scale = theta_hat[2])
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qgamma(c(1:n)/(n+1), shape = theta_hat_star[1], scale = theta_hat_star[2])
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "beta"){
    q_hat <- qbeta(c(1:n)/(n+1),shape1 = theta_hat[1], shape2 = theta_hat[2])
    
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rbeta(n, shape1 =  theta_hat[1],shape2 =  theta_hat[2])
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qbeta(c(1:n)/(n+1), shape1 =  theta_hat_star[1], shape2 = theta_hat_star[2])
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "exponential"){
    q_hat <- qexp(c(1:n)/(n+1),theta_hat)
    
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rexp(n, theta_hat)
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qexp(c(1:n)/(n+1), theta_hat_star)
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "geometric"){
    q_hat <- qgeom(c(1:n)/(n+1),prob = theta_hat)
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rgeom(n, theta_hat)
      theta_hat_star <- mle_name(x_star)
      q_hat_star <- qgeom(c(1:n)/(n+1), theta_hat_star)
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "binomial"){
    q_hat <- qbinom(c(1:n)/(n+1), n, prob = theta_hat)
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rbinom(n,n, theta_hat)
      theta_hat_star <- mle_name(x_star,n)
      
      q_hat_star <- qbinom(c(1:n)/(n+1), size=n, theta_hat_star)
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
  else if(distribution == "bernoulli"){
    q_hat <- qbinom(c(1:n)/(n+1), 1, prob = theta_hat)
    D0 <- ks.test(input_data, q_hat)$statistic
    Dvec<-NULL
    
    for(i in 1:nboot){
      x_star <- rbinom(n,1, theta_hat)
      theta_hat_star <- mle_name(x_star)
      
      q_hat_star <- qbinom(c(1:n)/(n+1), size=1, theta_hat_star)
      D_star <- ks.test(x_star, q_hat_star)$statistic
      Dvec <- c(Dvec, D_star)
    }
    #print(Dvec)
    p_value <- sum(Dvec > D0)/nboot
    return(p_value)
  }
}

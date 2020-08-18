## Testing the MLE estimated parameter by sending data generated from a particular distribution
## Also, determining the p_value by doing parametric bootstrap of MLE using K-test.

#source('~/Documents/581_Probability_Statistics/MLE/MLE.R')
#source('~/Documents/581_Probability_Statistics/MLE/KS_test.R')

mle_test_function <-function(distribution, num_samples = 100)
{
   
  
  if(distribution=="bernoulli")
  {
    p = 0.6
    print("-------Bernoulli-------")
    print(paste("Population parameter: ", p))
    input_data <- rbinom(num_samples,1,p) 
    p_mle <- bernoulli_mle(input_data)
    print(paste("Estimated parameter: ", p_mle))

    return(p_mle)
  }
  else if(distribution=="binomial")
  {
    n = 10
    p = 0.6
    print("-------Binomial-------")
    print(paste("Population parameter: ", p))
    input_data <- rbinom(num_samples,n,p) 
    p_mle <- binomial_mle(input_data,n)
    print(paste("Estimated parameter: ", p_mle))
    
    return(p_mle)
  }
  else if(distribution=="geometric")
  {
    p = 0.6
    print("-------Geometric-------")
    print(paste("Population parameter: ", p))
    input_data <- rgeom(num_samples,p) 
    p_mle <- geometric_mle(input_data)
    print(paste("Estimated parameter: ", p_mle))

    # Doing parametric bootstrap of MLE using ks test
    p_value <- ks_test_func(distribution, input_data = input_data, theta_hat = p_mle)
    print(paste("p-value: ", p_value))

    return(p_mle)
  }
  else if(distribution=="poisson")
  {
    lamda = 0.02
    print("-------Poisson-------")
    print(paste("Population parameter: ", lamda))
    
    input_data <- rpois(num_samples,lamda) 
    lamda_mle <- poisson_mle(input_data)
    print(paste("Estimated parameter: ", lamda_mle))
    
    return(lamda_mle)
  }
  
  else if(distribution=="uniform")
  {
    a = 0
    b = 5
    print("-------Uniform-------")
    print(paste("Population parameter: ", a, b))
   
    input_data <- runif(num_samples,a,b) 
    theta_mle <- uniform_mle(input_data)
    print(paste("Estimated parameter: ", theta_mle))

    # Doing parametric bootstrap of MLE using ks test
    p_value <- ks_test_func(distribution, input_data = input_data, theta_hat =theta_mle)
    print(paste("p-value: ", p_value))

    return(theta_mle)
  }
  
  else if(distribution=="normal")
  {
    mean = 0
    variance = 1
    print("-------Normal-------")
    print(paste("Population parameter: ", "mean = ", mean,"variance = ",variance))
    
    input_data <- rnorm(num_samples,mean,variance) # mean = 0, variance = 1
    theta_mle <- normal_mle(input_data)
    print(paste("Estimated parameter: ", theta_mle))
    
    # Doing parametric bootstrap of MLE using ks test
    p_value <- ks_test_func(distribution, input_data = input_data, theta_hat =theta_mle)
    print(paste("p-value: ", p_value))
    
    return(theta_mle)
  }
  
  else if(distribution=="exponential")
  {
    beta = 0.2
    print("-------Exponential-------")
    print(paste("Population parameter: ", beta))
    
    input_data <- rexp(num_samples,beta) # beta = 6
    beta_mle <- exponential_mle(input_data)
    print(paste("Estimated parameter: ", beta_mle))
    
    # Doing parametric bootstrap of MLE using ks test
    p_value <- ks_test_func(distribution, input_data = input_data, theta_hat =beta_mle)
    print(paste("p-value: ", p_value))
    
    return(beta_mle)
  }
  else if(distribution=="gamma")
  {
    alpha = 10 
    beta = 0.4
    print("-------Gamma-------")
    print(paste("Population parameter: ", alpha,beta))
    
    input_data <- rgamma(num_samples,alpha,scale = beta) # alpha = 2, beta = 0.4
    theta_mle <- gamma_mle(input_data)
    print(paste("Estimated parameter: ", theta_mle))
    
    # Doing parametric bootstrap of MLE using ks test
    p_value <- ks_test_func(distribution, input_data = input_data, theta_hat =theta_mle)
    print(paste("p-value: ", p_value))
    
    return(theta_mle)
  }
  else if(distribution=="beta")
  {
    alpha = 10
    beta = 5
    print("-------Beta-------")
    print(paste("Population parameter: ", alpha,beta))
    
    input_data <- rbeta(num_samples, alpha, beta) # alpha = 3, beta = 5
    theta_mle <- beta_mle(input_data)
    print(paste("Estimated parameter: ", theta_mle))
    
    # Doing parametric bootstrap of MLE using ks test
    p_value <- ks_test_func(distribution, input_data = input_data, theta_hat =theta_mle)
    print(paste("p-value: ", p_value))
    
    return(theta_mle)
  }
  else if(distribution=="Chi-square")
  {
    v = 4
    print("-------Chi-square-------")
    print(paste("Population parameter: ", v))
    
    input_data <- rchisq(num_samples,v) # degree of freedom (v = 3)
    v_mle <- chi_sqaure_mle(input_data)
    print(paste("Estimated parameter: ", v_mle))
   
    return(v_mle)
  }
  
  else
  {
    return("Wrong distribution")
  }
  
}


num_samples = 100  #Number of samples generated for each distribution
### Test the MLE parameters and the KS_test of each distribution by running this code:
mle_test_function("bernoulli")
mle_test_function("binomial")
mle_test_function("geometric")
mle_test_function("exponential")
mle_test_function("poisson")
mle_test_function("normal")
mle_test_function("uniform")
mle_test_function("gamma")
mle_test_function("beta")
mle_test_function("Chi-square")

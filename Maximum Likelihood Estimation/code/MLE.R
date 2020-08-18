### Maximum Likelihood estimation for different distributions ### 

### These R file contains all the functions of the distribution. Each distribution takes the 
### generated data as an input and returns the estimated parameters of the distribution.


# Bernoulli distribution: estimate (probability p)
bernoulli_mle <- function(input_data)
{
  p_mle <- mean(input_data)
  return(p_mle)
}

# Binomial distribution: estimate(probability p) for a given n (assume that n is known to us)
binomial_mle <-function(input_data,n)
{
  m <- length(input_data)
  p_mle <- mean(input_data)/n
  return(p_mle)
}


# Poisson distribution: estimate(lamda)
poisson_mle <- function(input_data)
{
  lamda_mle <- mean(input_data)
  return(lamda_mle)
}

# Geometric distribution: estimate(probaility p)
geometric_mle <- function(input_data)
{
  p_mle <- (1.0/(mean(input_data)+1))
  return(p_mle)
}

# Exponential distribution: estimate(beta)
exponential_mle <- function(input_data)
{
  beta_mle <- 1/mean(input_data)
  return(beta_mle)
}

# Normal distribution : estimate(mean,variance)
normal_mle <-function(input_data)
{
  n <- length(input_data)
  mean_mle <- mean(input_data)
  variance_mle <- sum((input_data-mean_mle)^2)/n
  return(c(mean_mle,variance_mle))
}

# Uniform distribution : estimate(a,b)
uniform_mle <-function(input_data)
{
  a_mle <- min(input_data)
  b_mle <- max(input_data)
  return(c(a_mle,b_mle))
}

# Gamma distribution : estimate(alpha,beta)

### Taking the log-liklihood of gamma function and taking its derivative w.r.t alpha yields:
### ln(alpha) - digamma(alpha) - log(sample_mean) + mean(log(input_data)) = 0

### Since there is no closed form solution for alpha, I have used the following approximation to estimate alpha_mle:
### ln(alpha) - digamma(alpha) ~ (1/(2*alpha))*(1 + 1/6*alpha +1)
### s = log(mean(input_data)) - mean(log(input_data))
### Thus, alpha_mle ~ ((3-s) + sqrt((s-3)**2 + 24*s))/12*s

### beta_mle is simply estimated from alpha_mle as: sample_mean/alpha_mle

gamma_mle <-function(input_data)
{
  s <- log(mean(input_data)) - mean(log(input_data))
  alpha_mle <- (3-s)/(12*s) + sqrt((s-3)^2 + 24*s)/(12*s)
  beta_mle <- mean(input_data)/alpha_mle
  return(c(alpha_mle,beta_mle))
}

# Beta distribution : estimate(a,b)

## Since the MLE derivatives for a and b are not closed form equations, 
## we use the Newton-Raphsen method to iteratively find the MLE estimations for a and b.

beta_mle <- function(input_data)
{
  
  # Initializing alpha and beta parameters using Method of moment estimates
  input_data_mean <- mean(input_data)
  input_data_variance <- (sum(input_data * input_data))/length(input_data)
  a_mle <- ((input_data_mean ^ 2) - (input_data_mean * input_data_variance))/(input_data_variance - (input_data_mean ^ 2))
  b_mle <- (a_mle * (1 - input_data_mean))/(input_data_mean)
  
  final_val <- c(a_mle, b_mle)
  
  # Running the optimisation step for 100 iterations
  for(index in 1:100){
    g1 <- digamma(a_mle) - digamma(a_mle + b_mle) - (sum(log(input_data)))/length(input_data)
    g2 <- digamma(b_mle) - digamma(a_mle + b_mle) - (sum(log(1 - input_data))/length(input_data))
    g <- c(g1, g2)
    
    # Calculating G matrix of second derivatives:
    G1_val <- trigamma(a_mle) - trigamma(a_mle + b_mle)
    G2_val <- -trigamma(a_mle + b_mle)
    G3_val <- trigamma(b_mle) - trigamma(a_mle + b_mle)
    G <- matrix(c(G1_val, G2_val, G2_val, G3_val), nrow = 2, ncol = 2, byrow = TRUE)
    G_inverse <- solve(G)
    
    # Final values for the iteration: Theta_mle(i+1) = Theta_mle(i) - G_inverse*g
    
    final_val <- final_val - t(G_inverse %*% g)
    a_mle <- final_val[1]
    b_mle <- final_val[2]
  }
  return(c(a_mle,b_mle))
}

# Chi-squared distribution: estimate(v)
chi_sqaure_mle <- function(input_data)
{
  ## Here I have used the property of chi-square being a special case of Gamma distribution.
  ## A gamma distribution (alpha,beta) with alpha = v/2 and beta = 1/2, 
  ## is a chi-squared RV with v degrees of freedom.
  
  result_mle <- gamma_mle(input_data)
  v_mle <- 2* result_mle[1]
  # print(result_mle)
  # print(v_mle)
  return(v_mle)
}
  





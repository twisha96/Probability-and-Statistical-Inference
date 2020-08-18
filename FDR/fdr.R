fdr <- function(p_values, Q){
  
  sorted_p_values = sort(p_values) # Sort the p-values 
  m = length(sorted_p_values)      # Number of p-values
  
  # Hypothesis 1
  hypothesis_1 = Q*c(1:m)/m
  
  # Hypothesis 2 - If not independent
  d = m * sum(1/c(1:m))
  hypothesis_2 = (Q * c(1:m))/d
  
  # the p-vlaues less than the hypothesis line are considered to be interesting
  is_interesting = (sorted_p_values < hypothesis_1)
  
  # compute the index below which the sorted p-values are less than the line
  index = max(which(is_interesting == "TRUE"))
  print(paste0("Threshold Index: ", index))
  
  #################################################################################
  # NOTE: The index value we get depends on the given input vector
  #
  # Here, the input vector has 100 values in the range of 1e-6 and rest are large.
  # Thus the index turns out to some value near 100.
  #
  # Suppose, the input vector was something with 250 smaller values.
  # Example: v1 <- c(1e-5*runif(250),runif(750))
  # The index would be approximately 250.
  #################################################################################
  
  p_star = sorted_p_values[index]
  
  hypothesis = c(1:m)
  plot(hypothesis, sorted_p_values, col="black")
  lines(c(1:index), sorted_p_values[c(1:index)], col = "red", type="o")
  lines(hypothesis, hypothesis_1, type="l", col="black")
  
  # list index of hypothesis which are interesting in the original unsorted list of p values
  temp = (p_values < p_star)
  fd = which(temp == "TRUE")
  print("Index of intersting p values in the unsorted list:")
  print(fd)
  
  # false rejection rate
  frr = -log(p_star)/length(fd)
  print(paste0("False Rejection Rate: ", frr))
  
  #################################################################################
  # NOTE: FRR depends on the Q value given as an input to the function
  # Here, the Q value was 0.05 hence the frr is approximately 0.05
  #################################################################################
}

# Test Vector of p-values
v1 <- c(1e-5*runif(100),runif(900))
Q = 0.05
fdr(v1, Q)

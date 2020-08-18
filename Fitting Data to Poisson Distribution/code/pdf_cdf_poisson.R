library(MASS)
library(tidyverse)

pdf <- dpois(c(0:20), lambda)
plot(pdf, main = "PDF of poisson distribution", xlab = "data", ylab = "Probability")

df["pdf"] = pdf

df <- df %>% mutate(cdf = cumsum(pdf))
plot(df[["cdf"]], main = "CDF of poisson distribution", xlab = "data", ylab = "Cumulative Probability" )


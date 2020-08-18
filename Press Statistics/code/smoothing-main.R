#### Main fil to run the smoothing and bootstrap functions, defined in the R file "smoothing_functions.R" ####

# read the data
NOAA1 <- read.csv("NOAA+GISS.csv")

###### Output of all the smoothing ####
plot(NOAA1[,3], NOAA1[,2], xlab="temperature rise", ylab="rate of billion dollar weather disasters")

dum <- bin.mean(NOAA1[,3], NOAA1[,2], 6)
dum <- gauss.mean(NOAA1[,3],NOAA1[,2],.063)$df
gauss.reg(NOAA1[,3],NOAA1[,2],.078,do.plot=T)$df
gauss.mean.trunc(NOAA1[,3],NOAA1[,2],.063,20,do.plot=T)$df
gauss.reg.trunc(NOAA1[,3],NOAA1[,2],.08,17,do.plot=T)$df
# using the inbuilf funtion for Lowess and Smoothing spline smoother
lines(lowess(NOAA1[,3],NOAA1[,2]),col=7)
lines(smooth.spline(NOAA1[,3],NOAA1[,2]),col=8)


#### boostraping the press vector and calculate the confidence interval ####
bootstrap_press <- bootstrap(NOAA1[,3],NOAA1[,2])
#print(bootstrap_press$pressvec)
plot(density(bootstrap_press$pressvec), main = "PRESS vector bootstrap")
print("With confidence interval: ")
print(bootstrap_press$CI)
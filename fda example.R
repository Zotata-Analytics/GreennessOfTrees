# Load necessary libraries
library(fda)

# Load the data
data <- t(as.matrix(read.csv("greenness_data.csv", header = TRUE)))

############################
## Visualise the raw data ##
############################

daytime <- seq(4, 364, 8)
matplot(daytime, data, type='l', 
        xlab='Day', ylab='Greenness of Tree', 
        main='Greenness of Trees for 200 Locations\n over 1 year at 8-day intervals')

#############################################################
## Conduct Generalised Cross Validation with a             ##
## saturated B-spline basis and standard roughness penalty ##
#############################################################

# use standard roughness penalty
day8 = seq(4,364,8)
dayrng = c(4,364)
knots <- day8
norder <- 4
nbasis <- length(knots) + norder - 2

# basis functions
bbasis <- create.bspline.basis(dayrng, nbasis, norder, knots)

curv.Lfd <- int2Lfd(2)
lambda <- 1e1
curv.fdPar <- fdPar(bbasis,curv.Lfd,lambda)

lambdas <- 10^seq(-4,4,by=0.5)
mean.gcv <- rep(0,length(lambdas))

# GCV Process
for(ilam in 1:length(lambdas)){
  
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[ilam]
  
  dataSmoothi = smooth.basis(daytime,data,curv.fdPari)
  
  mean.gcv[ilam] = mean(dataSmoothi$gcv)
}

# GCV Curve
plot(lambdas,mean.gcv,type='b',log='x')

best <- which.min(mean.gcv)
lambdabest <- lambdas[best] # 3.16

curv.fdPar$lambda <- lambdabest
dataSmooth <- smooth.basis(daytime,data,curv.fdPar)

# plot smooth functions
par(mfrow=c(2,2))
plotfit.fd(data,daytime,dataSmooth$fd,index=1:4)
dataSmooth$df

# plot final composite visualisation
par(mfrow=c(1,1))
plot(dataSmooth)

####################################################################
## Conduct Generalised Cross Validation with a saturated B-spline ##
## basis and harmonic acceleration roughness penalty              ##
####################################################################

## use harmonic acceleration penalty instead
day8 = seq(4,364,8)
dayrng = c(4,364)
knots <- day8
norder <- 5
nbasis <- length(knots) + norder - 2
bbasis <- create.bspline.basis(dayrng, nbasis, norder, knots)

period <- 365
omega <- (2*pi/period)^2
harm.Lfd <- vec2Lfd(c(0, omega, 0), dayrng)
lambda <- 1e1
harm.fdPar <- fdPar(bbasis, harm.Lfd, lambda)

lambdas <- 10^seq(-4, 4, by=0.5)
mean.gcv <- rep(0, length(lambdas))

# GCV Process
for (ilam in 1:length(lambdas)) {
  harm.fdPari <- harm.fdPar
  harm.fdPari$lambda <- lambdas[ilam]
  
  dataSmoothi <- smooth.basis(daytime, data, harm.fdPari)
  
  mean.gcv[ilam] <- mean(dataSmoothi$gcv)
}

# GCV Curve
plot(lambdas, mean.gcv, type='b', log='x')

best <- which.min(mean.gcv)
lambdabest <- lambdas[best] # 316

harm.fdPar$lambda <- lambdabest
dataSmooth <- smooth.basis(daytime, data, harm.fdPar)

# plot smooth functions
par(mfrow=c(2,2))
plotfit.fd(data, daytime, dataSmooth$fd, index=1:4)
dataSmooth$df

# plot final composite visualisation
par(mfrow=c(1,1))
plot(dataSmooth)

############################################
## Obtain the first and second derivative ##
############################################

# Plot the graphs of the first and second derivatives of the curves
par(mfrow=c(1,1))
plot(deriv.fd(dataSmooth$fd))
plot(deriv.fd(dataSmooth$fd,2))

#####################################################################
## Perform an unpenalised Principal Component Analysis of the data ##
#####################################################################

datafd <- dataSmooth$fd
dataPCA <- pca.fd(datafd,nharm=6)

par(mfrow=c(1,1))
plot(dataPCA$varprop,type='b') # need at least 3 harmonics to explain at least 80% of the variance
cumsum(dataPCA$varprop)

par(mfrow=c(1,1))
plot(dataPCA$harmonics[1:3])

par(mfrow=c(2,2))
plot(dataPCA,harm=1:3)

###############################################################################
## Use GCV to find the best smoothing parameter for the principal components ##
###############################################################################

## Example of smoothed PCs with a lambda of 1e4
pca.fdPar <- fdPar(bbasis,curv.Lfd,1e4)
dataPCAsmooth <- pca.fd(datafd,nharm=6,harmfdPar=pca.fdPar)

par(mfrow=c(1,1))
plot(dataPCAsmooth$harmonics[1:3])

par(mfrow=c(2,2))
plot(dataPCAsmooth)

#############################

# use standard roughness penalty and GCV to find optimal lambda
dayrng = c(4, 364)
daytime = seq(4, 364, 8)
knots <- seq(4, 364, 8)
norder <- 4
nbasis <- length(knots) + norder - 2
bbasis <- create.bspline.basis(dayrng, nbasis, norder, knots)

lambdas <- 10^seq(-4, 4, by=0.5)
mean.gcv <- rep(0, length(lambdas))

# GCV Process
for (ilam in 1:length(lambdas)) {
  fdParObj <- fdPar(bbasis, Lfdobj=int2Lfd(2), lambda=lambdas[ilam])
  dataSmoothi <- smooth.basis(daytime, data, fdParObj)
  mean.gcv[ilam] <- mean(dataSmoothi$gcv)
}

# GCV plot
plot(lambdas, mean.gcv, type='b', log='x', xlab='Lambda', ylab='GCV Score',
     main='GCV Scores for Different Lambda Values')

best.lambda.index <- which.min(mean.gcv)
optimal.lambda <- lambdas[best.lambda.index] # 3.16

optimal.fdPar <- fdPar(bbasis, Lfdobj=int2Lfd(2), lambda=optimal.lambda)
dataSmooth <- smooth.basis(daytime, data, optimal.fdPar)
datafdSmooth <- dataSmooth$fd

# Conduct PCA
dataPCA <- pca.fd(datafdSmooth, nharm=6)

# variance explained
par(mfrow=c(1,1))
plot(dataPCA$varprop, type='b', xlab='Principal Component', ylab='Variance Explained',
     main='Variance Explained by Principal Components')

# plot first three PCs
par(mfrow=c(1,1))
plot(dataPCA$harmonics[1:3])

#####################################################################
## Use Leave One Out Cross Validation (LOOCV) to                   ##
## find the best smoothing parameter for the prinicipal components ##
#####################################################################

# Decided to try LOOCV instead of GCV to see if it would give a 
# different optimal lambda value. The optimal lambda is very different 
# hence I will just include both results.

day8 = seq(4,364,8)
dayrng = c(4,364)
knots <- day8
norder <- 5
nbasis <- length(knots) + norder - 2
bbasis <- create.bspline.basis(dayrng, nbasis, norder, knots)

# LOOCV Function for PCA
performLOOCV <- function(data, day8, basis, Lfd, lambdas) {
  n <- ncol(data)
  cv_errors <- numeric(length(lambdas))
  
  for (i in 1:length(lambdas)) {
    lambda <- lambdas[i]
    fdParObj <- fdPar(basis, Lfd, lambda)
    smooth_fd <- smooth.basis(daytime, data, fdParObj)
    data_fd <- smooth_fd$fd
    errors <- numeric(n)
    
    for (j in 1:n) {
      data_loo <- data[,-j]
      smooth_fd_loo <- smooth.basis(daytime, data_loo, fdParObj)
      data_fd_loo <- smooth_fd_loo$fd
      
      pca_loo <- pca.fd(data_fd_loo, nharm=3)
      
      predict_fd <- eval.fd(daytime, pca_loo$harmonics)
      
      actual <- data[,j]
      predicted <- predict_fd[,1]
      errors[j] <- sqrt(mean((actual - predicted)^2))
    }
    
    cv_errors[i] <- mean(errors)
  }
  
  return(list(lambda=lambdas, cv_errors=cv_errors))
}

# Apply the LOOCV function to find optimal lambda
lambdas <- 10^seq(-4, 4, by=0.5)
cv_results <- performLOOCV(data, day8, bbasis, int2Lfd(2), lambdas)
optimal_lambda <- lambdas[which.min(cv_results$cv_errors)]
print(paste("Optimal lambda:", optimal_lambda)) # 1000

# Plot CV errors
plot(lambdas, cv_results$cv_errors, type='b', log='x', xlab='Lambda', ylab='CV Error')
title('LOOCV Errors across Lambda Values')

# Smoothed PCA with optimal lambda
optimal_fdPar <- fdPar(bbasis, int2Lfd(2), optimal_lambda)
dataSmooth <- smooth.basis(daytime, data, optimal_fdPar)
data_fd <- dataSmooth$fd
dataPCA_optimal <- pca.fd(data_fd, nharm=3)

# visualise results
par(mfrow=c(2,2))
plot(dataPCA_optimal, harm=1:3)

# plot first three PCs
par(mfrow=c(1,1))
plot(dataPCA_optimal$harmonics[1:3])

##############################################################################
## Checking orthogonality of unpenalised and penalised principal components ##
##############################################################################

## Unpenalised PCA ##
       
# Create a matrix to hold the inner products
inprod_unpen <- matrix(nrow = 4, ncol = 4)

# Calculate the inner product for each pair of principal components
for(i in 1:4){
  for(j in 1:4){
    inprod_unpen[i,j] <- inprod(dataPCA$harmonics[i], dataPCA$harmonics[j])
  }
}

# Print the matrix of inner products for unpenalised components
print("Inner products of unpenalised principal components:")
print(inprod_unpen)

## Penalised PCA ##

## Standard Roughness Penalty##
penalisedDatafd <- curv.fdPar

# penalised PCA on the smoothed data
penalisedPCA_spline <- pca.fd(penalisedDatafd, nharm=4)

# Create a matrix to store the inner products
inprod_pen_spline <- matrix(nrow = 4, ncol = 4)

# Calculate the inner product for each pair of penalised principal components
for(i in 1:4){
  for(j in 1:4){
    inprod_pen_spline[i,j] <- inprod(penalisedPCA_spline$harmonics[i], penalisedPCA_spline$harmonics[j])
  }
}

# Print the matrix of inner products for penalised components (Spline)
print("Inner products of penalised principal components (Spline):")
print(inprod_pen_spline)

## Harmonic Accel Penalty ##
penalisedDatafd <- harm.fdPar

# Penalised PCA on the smoothed data
penalisedPCA_harm <- pca.fd(penalisedDatafd, nharm=4)

# Create a matrix to store the inner products
inprod_pen_harm <- matrix(nrow = 4, ncol = 4)

# Calculate the inner product for each pair of penalised principal components
for(i in 1:4){
  for(j in 1:4){
    inprod_pen_harm[i,j] <- inprod(penalisedPCA_harm$harmonics[i], penalisedPCA_harm$harmonics[j])
  }
}

# Print the matrix of inner products for penalised components (Harmonic Acceleration)
print("Inner products of penalised principal components (Harmonic Acceleration):")
print(inprod_pen_harm)

# Using the inprod function from the fda library numerically show that the un-penalised
# principal components are orthogonal, but the penalised principal components might not be.
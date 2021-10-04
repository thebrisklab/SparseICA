# Other packages may be required: 
# install.packages("steadyICA")
# install.packages("neuRosim")
# install.packages("ProDenICA")

source("./Code/RCode/refinedICA.R")
source('./Code/RCode/jngcaFunctions.R')

# You can directly use the provided datasets
# xmat = as.matrix(read.csv('./Data/X.csv'))
# smat = as.matrix(read.csv('./Data/S.csv'))

simData = SimFMRI123(var.inactive = 0, snr=0.2) 
xmat = simData$whitened.X
smat = simData$S
  
# note: here the data matrix x is prewhitened:
cov(xmat)[1:10,1:10]

# examine the first three time points:
par(mfrow=c(1,3))
image(matrix(xmat[,1],33))
image(matrix(xmat[,2],33))
image(matrix(xmat[,3],33))

# examine the true components:
# here, there are 3 components:
dim(smat)
par(mfrow=c(1,3))
image(matrix(smat[,1],33,33))
image(matrix(smat[,2],33))
image(matrix(smat[,3],33))

# Check FastICA method with logistic density
estX_logis = mlcaFP(xData = xmat, n.comp = 3, whiten = "none", restarts.pbyd = 30, distribution='logistic')
par(mfrow=c(1,3))
image(matrix(estX_logis$S[,1],33,33))
image(matrix(estX_logis$S[,2],33,33))
image(matrix(estX_logis$S[,3],33,33))
 
# Check FastICA method with tanh density
estX_tanh = mlcaFP(xData = xmat, n.comp = 3, whiten = "none", restarts.pbyd = 1, distribution='tanh')
par(mfrow=c(1,3))
image(matrix(estX_tanh$S[,1],33,33))
image(matrix(estX_tanh$S[,2],33,33))
image(matrix(estX_tanh$S[,3],33,33))

# Different measurements to check the recovered components
estX_tanh$M = est.M.ols(sData = estX_tanh$S,xData = simData$X)
frobICA(S1 = estX_tanh$S, S2 = smat, standardize=TRUE)


#####################
# Simulation (Examples to use our codes)
# Four different algorithms are allowed in one function: refinedICA(...)
source("./Code/RCode/refinedICA.R")
source('./Code/RCode/jngcaFunctions.R')

simData = SimFMRI123(var.inactive = 0, snr=0.2) 
xmat = simData$whitened.X
smat = simData$S
Q = 3

# examine the first three time points:
par(mfrow=c(1,3))
image(matrix(xmat[,1],33))
image(matrix(xmat[,2],33))
image(matrix(xmat[,3],33))

# e.g.: Do 10 Simulations 
itm = 10 

## Relax logistic 
## (Typically, this one is executed slowly, especially compared with others)
relax_logistic_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 2, method = "relax_logistic", lambda = 0.5, nu = 1, converge_plot = F)
  relax_logistic_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, relax_logistic_diff, ylim = c(0,1.2),
     main = "Relax-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Logistic Distance")

## Relax Softmax 
## (This runs fast, and it can allow sparsity results)
relax_softmax_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 5, method = "relax_laplace", lambda = 0.5, nu = 1, converge_plot = F)
  relax_softmax_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, relax_softmax_diff, ylim = c(0,1.2),
     main = "Relax-Softmax", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Softmax Distance")

## FastICA logistic 
## (This is just to run the FastICA codes with logistic density)
FastICA_logistic_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 1, method = "FastICA_logistic", lambda = 0.5, nu = 1)
  FastICA_logistic_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, FastICA_logistic_diff, ylim = c(0,1.2),
     main = "FastICA-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-Logistic Distance")

## FastICA tanh 
## (This is just to run the FastICA codes with tanh density)
FastICA_tanh_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 1, method = "FastICA_tanh", lambda = 0.5, nu = 1)
  FastICA_tanh_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, FastICA_tanh_diff, ylim = c(0,1.2),
     main = "FastICA-tanh", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-tanh Distance")

## For relax_laplace and relax_logistic, you can check the convergence plot by using converge_plot = T.
## There is no convergence plot for FastICA methods, because we just embedded their codes that didn't allow the convergence plot. 
estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 5, method = "relax_laplace", lambda = 0.5, nu = 1, converge_plot = T)
estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 2, method = "relax_logistic", lambda = 0.5, nu = 1, converge_plot = T)

## Combine the above plots into one. Need to change the margins for a good plot.
par(mfrow = c(2,2),mar=c(2,2,4,1))
plot(1:itm, relax_logistic_diff, ylim = c(0,1.2),
     main = "Relax-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Logistic Distance")
plot(1:itm, relax_softmax_diff, ylim = c(0,1.2),
     main = "Relax-Softmax", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Softmax Distance")
plot(1:itm, FastICA_logistic_diff, ylim = c(0,1.2),
     main = "FastICA-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-Logistic Distance")
plot(1:itm, FastICA_tanh_diff, ylim = c(0,1.2),
     main = "FastICA-tanh", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-tanh Distance")


## Show the recovered components images
## Some didn't recover the components correctly. 
par(mfrow=c(1,3))
image(matrix(estX$estS[,1],33,33))
image(matrix(estX$estS[,2],33,33))
image(matrix(estX$estS[,3],33,33))

## Make boxplot
library(ggplot2)
relax_logistic_diff_df = data.frame(Distance = relax_logistic_diff)
relax_logistic_diff_df$Method = "Relax-Logistic"
relax_softmax_diff_df = data.frame(Distance = relax_softmax_diff)
relax_softmax_diff_df$Method = "Relax-Softmax"
FastICA_logistic_diff_df = data.frame(Distance = FastICA_logistic_diff)
FastICA_logistic_diff_df$Method = "mlcaFP-Logistic"
FastICA_tanh_diff_df = data.frame(Distance = FastICA_tanh_diff)
FastICA_tanh_diff_df$Method = "mlcaFP-tanh"
diff = rbind(relax_logistic_diff_df,relax_softmax_diff_df,FastICA_logistic_diff_df,FastICA_tanh_diff_df)

ggplot(diff, aes(x = Method, y = Distance))+
  geom_boxplot(outlier.alpha = 0.5)


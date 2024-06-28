
# /////////////////////// Step 1: Loading Dataset and Library  /////////////////

# Load necessary libraries
# install.packages('depthTools')
library(depthTools)
library(tidyverse)
library(fda)
library(readr)

# path to the PERG data and participants info
RE_data_path <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/signal_RE_1.csv"
LE_data_path <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/signal_LE_1.csv"
participants_time <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/microsecond_df.csv"

time_microsecond_numeric <- read.csv(participants_time)
data_frame_RE <- read.csv(RE_data_path)
data_frame_LE <- read.csv(LE_data_path)

# Remove the first column which is time
data_RE <- data_frame_RE[,-1] 
data_LE <- data_frame_LE[,-1] 

# Convert the data frame to a matrix and transpose it
data_on_RE <- t(as.matrix(data_RE)) # transpose to match dimensions
data_on_LE <- t(as.matrix(data_LE))

# time sequence
time_micro <- data_frame_RE[,1]
time_micro

# Plot the change in signals data using matplot
matplot(x = time_micro, y = data_RE, type='l', lty=1, xlab='MicroSecond', ylab='Change in PERG signals', main='Right eye: RE_1')
matplot(x = time_micro, y = data_LE, type='l', lty=1, xlab='MicroSecond', ylab='Change in PERG signals', main='Left eye: LE_1')


# /////////////////////// Step 2: Smoothing functional data  /////////////////

# .....................B-spline ..................

# Create a B-spline basis for modeling changes in signals
def_basis <- create.bspline.basis(c(0,1499),norder=4)
plot(def_basis)

# Plotting and Evaluate the fitted data for diagnosis_levels = [2] and [3]for 
# 'Bilateral optic nerve atrophy','Autoimmune retinopathy', both LE and RE.

# WE do analysis separately for both RE_1(Signals of Right eye) and LE_1(Signals of Left eye)

# RE_1    (Signals of Right eye)

t_data_on_RE = t(data_on_RE)
t_data_on_RE[,2]

matrix_t_data_on_RE = as.matrix(t_data_on_RE[,2])
matrix_t_data_on_RE
funct_RE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_RE, basisobj = def_basis)
funct_RE_fitted <- eval.fd(time_micro, funct_RE)
# Plotting fitted data and actual data points in one pane
par(mfrow = c(1, 1))
plot(time_micro, funct_RE_fitted, type='l', ylim=c(-3.5, 3.5), xlab='MicroSecond', ylab='UN_fited and Fitted Signals for RE_1', main='RE_1 for Bilateral optic nerve atrophy')
points(time_micro, t_data_on_RE[,2], col='blue', pch=20)  

# LE_1  (Signals of Left eye)

t_data_on_LE = t(data_on_LE)
t_data_on_LE[,2]
matrix_t_data_on_LE = as.matrix(t_data_on_LE[,2])
matrix_t_data_on_LE
funct_LE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_LE, basisobj = def_basis)
funct_LE_fitted <- eval.fd(time_micro, funct_LE)
# Plotting fitted data and actual data points in one pane
par(mfrow = c(1, 1))
plot(time_micro, funct_LE_fitted, type='l', ylim=c(-4.0, 4.0), xlab='MicroSecond', ylab='UN_fited and Fitted Signals for LE_1', main='LE_1 for Bilateral optic nerve atrophyy')
points(time_micro, t_data_on_LE[,2], col='blue', pch=20)  


# In This section we will present with impact of the number of orders on fitting curve

# RE_1

# Create a for loop to iterate over each norder value for second diagnosis_levels = [2] and RE_1
norder_values <- 3:13  # range of norder values
for (norder in norder_values){
  def_basis_RE <- create.bspline.basis(c(0,1499), norder = norder)
  matrix_t_data_on_RE <- as.matrix(t_data_on_RE[,2])
  funct_RE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_RE, basisobj = def_basis_RE)
  funct_RE_fitted <- eval.fd(time_micro, funct_RE)
  
  # Plot fitted data and actual data points in each pane
  par(mfrow = c(1, 1))
  plot(time_micro, funct_RE_fitted, type = 'l', ylim = c(-3.5, 3.5), xlab = 'MicroSecond', ylab = 'UN_fited and Fitted Signals for RE_1')
  points(time_micro, t_data_on_RE[,2], col = 'blue', pch = 20)
  # Add title indicating the current value of norder
  title(main = paste("RE_1:norder =", norder))
}

# result: norder with 7 to 12 are good

# LE_1

norder_values <- 3:14
# Create a for loop to iterate over each norder value for second diagnosis_levels = [1] and LE_1
for (norder in norder_values){
   def_basis <- create.bspline.basis(c(0,1499), norder = norder)
   matrix_t_data_on_LE <- as.matrix(t_data_on_LE[,2])
   funct_LE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_LE, basisobj = def_basis)
   funct_LE_fitted <- eval.fd(time_micro, funct_LE)

   # Plot fitted data and actual data points in each pane
   par(mfrow = c(1, 1))
   plot(time_micro, funct_LE_fitted, type = 'l', ylim = c(-4.5, 4.5), xlab = 'MicroSecond', ylab = 'UN_fited and Fitted Signals for LE_1')
   points(time_micro, t_data_on_LE[,2], col = 'blue', pch = 20)
  # Add title indicating the current value of norder
   title(main = paste("LE_1: norder =", norder))
}

# result: norder with 14 is best



# ..................... Fourier basis............................

# ... RE_1 

# nbasis=4 .........
fourier_basis<-create.fourier.basis(rangeval = c(0,1499), nbasis=4, period = 255)
funct_RE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_RE, basisobj = fourier_basis)
funct_RE_fitted <- eval.fd(time_micro, funct_RE)
# Plotting fitted data and actual data points
plot(time_micro, funct_RE_fitted, type='l', ylim=c(-3.5, 3.5), xlab='MicroSecond', ylab='Fitted data RE_1')
points(time_micro, matrix_t_data_on_RE, col='red')
title(main = "RE_1: Fourier basis with nbasis = 4")

# nbasis=15 ........
fourier_basis<-create.fourier.basis(rangeval = c(0,1499), nbasis=15, period = 255)
funct_RE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_RE, basisobj = fourier_basis)
funct_RE_fitted <- eval.fd(time_micro, funct_RE)
# Plotting fitted data and actual data points
plot(time_micro, funct_RE_fitted, type='l', ylim=c(-3.5, 3.5), xlab='MicroSecond', ylab='Fitted data RE_1')
points(time_micro, matrix_t_data_on_RE, col='red')
title(main = "RE_1: Fourier basis with nbasis = 15")

# result: increasing the number of basis from 4 to 15 does not effect

# ... LE_1 

# nbasis=4 .........
fourier_basis<-create.fourier.basis(rangeval = c(0,1499), nbasis=4, period = 255)
funct_RE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_LE, basisobj = fourier_basis)
funct_RE_fitted <- eval.fd(time_micro, funct_RE)
# Plotting fitted data and actual data points
plot(time_micro, funct_RE_fitted, type='l', ylim=c(-3.5, 3.5), xlab='MicroSecond', ylab='Fitted data LE_1')
points(time_micro, matrix_t_data_on_LE, col='red')
title(main = "LE_1: Fourier basis with nbasis = 4")

# nbasis=15 ........
fourier_basis<-create.fourier.basis(rangeval = c(0,1499), nbasis=15, period = 255)
funct_RE <- Data2fd(argvals = time_micro, y = matrix_t_data_on_LE, basisobj = fourier_basis)
funct_RE_fitted <- eval.fd(time_micro, funct_RE)
# Plotting fitted data and actual data points
plot(time_micro, funct_RE_fitted, type='l', ylim=c(-3.5, 3.5), xlab='MicroSecond', ylab='Fitted data LE_1')
points(time_micro, matrix_t_data_on_LE, col='red')
title(main = "LE_1: Fourier basis with nbasis = 15")


# .....................B-spline with lambda penalty ..................

# Create a B-spline basis with lambda penalty term
norder = 12
nbasis = 255 + norder - 2   # use max number of basis
RE_basis = create.bspline.basis(c(0,1499), nbasis, norder, time_micro)
# Let’s fix a penalty using a derivative of order 4 and compute the smoothed fitting to the data, using lambda = 0.001
lam = 0.001
RE_fdPar = fdPar(RE_basis, Lfdobj = 4 ,lambda = lam)
RE_smooth = smooth.basis(time_micro, matrix_t_data_on_RE, RE_fdPar)$fd
# Evaluate the fitted Signals data at given time points
RE_fitted<-eval.fd(time_micro, RE_smooth)

matplot(time_micro, RE_fitted, type='l', lty=1, pch='o',ylim=c(-3.5, 3.5), xlab='MicroSecond',ylab='Fitted data RE_1')
points(time_micro, matrix_t_data_on_RE, col='purple',pch=20)
title(main = "B-spline basis with nbasis = 265, and lambda = 0.001")

# ... RE_1

# ...... GCV (Generalized Cross Validation)..........
# Let’s compute the GCV to see which values of lambda are optimal (or suboptimal).
loglam=seq(-6,0,by=0.25)
gcvsave=matrix(0,nrow=length(loglam),ncol = 1)

for (i in c(1:length(loglam))){
  lambdai=10^loglam[i]
  RE_fdPar=fdPar(RE_basis, Lfdobj = 4 ,lambda = lambdai)
  s = smooth.basis(time_micro, matrix_t_data_on_RE, RE_fdPar)
  gcvsave[i,] = s$gcv
}

plot(loglam, gcvsave, type='l', xlab = 'log(lambda)', ylab='GCV')
title(main = "GCV for RE_1")
# .......... best value for lambda ......
# The best value for lambda seems to be in between 0.0001 and 1. 
# Let’s try to fit again the data with lambda = 0.5
# Create a B-spline basis with lambda penalty term

norder = 12
nbasis = 255 + norder - 2   # use max number of basis
RE_basis = create.bspline.basis(c(0,1499), nbasis, norder, time_micro)
#Let’s fix a penalty using a derivative of order 4 and compute the smoothed fitting to the data, using lambda = 0.001
lam = 0.5
RE_fdPar = fdPar(RE_basis, Lfdobj = 4 ,lambda = lam)
RE_smooth = smooth.basis(time_micro, matrix_t_data_on_RE, RE_fdPar)$fd #extract the estimated functional data
# Evaluate the fitted Signals data at given time points
RE_fitted<-eval.fd(time_micro, RE_smooth) 

matplot(time_micro, RE_fitted, type='l', lty=1, pch='o',ylim=c(-3.5, 3.5), xlab='MicroSecond',ylab='Fitted data RE_1: Bilateral optic nerve atrophy')
points(time_micro, matrix_t_data_on_RE, col='purple',pch=20)
title(main = "RE_1: first col, B-spline, nbasis = 265, lambda = 0.5")


# ... LE_1

# ...... GCV (Generalized Cross Validation)..........
# Let’s compute the GCV to see which values of lambda are optimal (or suboptimal).
loglam=seq(-6,0,by=0.25)
gcvsave=matrix(0,nrow=length(loglam),ncol = 1)

for (i in c(1:length(loglam))){
  lambdai=10^loglam[i]
  RE_fdPar=fdPar(RE_basis, Lfdobj = 4 ,lambda = lambdai)
  s = smooth.basis(time_micro, matrix_t_data_on_LE, RE_fdPar)
  gcvsave[i,] = s$gcv
}

plot(loglam, gcvsave, type='l', xlab = 'log(lambda)', ylab='GCV')
title(main = "GCV for LE_1")
# .......... best value for lambda ......
# The best value for lambda seems to be in between 0.01 and 1. 
# Let’s try to fit again the data with lambda = 0.5
# Create a B-spline basis with lambda penalty term

norder = 12
nbasis = 255 + norder - 2   # use max number of basis
RE_basis = create.bspline.basis(c(0,1499), nbasis, norder, time_micro)
#Let’s fix a penalty using a derivative of order 4 and compute the smoothed fitting to the data, using lambda = 0.001
lam = 0.5
RE_fdPar = fdPar(RE_basis, Lfdobj = 4 ,lambda = lam)
RE_smooth = smooth.basis(time_micro, matrix_t_data_on_LE, RE_fdPar)$fd #extract the estimated functional data
# Evaluate the fitted Signals data at given time points
LE_fitted<-eval.fd(time_micro, RE_smooth)

matplot(time_micro, LE_fitted, type='l', lty=1, pch='o',ylim=c(-4.5, 4.5), xlab='MicroSecond',ylab='Fitted data LE_1: Bilateral optic nerve atrophy')
points(time_micro, matrix_t_data_on_LE, col='purple',pch=20)
title(main = "LE_1: first col, B-spline, nbasis = 265, lambda = 0.5")



# ..... Applying B-spline Smoothing on All observation in RE_1 and LE_1:


# .... All RE_1 Data:

# Create a B-spline basis with lambda penalty term " For All RE_1 Data
all_data_on_RE_mat = as.matrix(t_data_on_RE)

norder = 12
nbasis = 255 + norder - 2   # use max number of basis
RE_basis = create.bspline.basis(c(0,1499), nbasis, norder, time_micro)
#Let’s fix a penalty using a derivative of order 4 and compute the smoothed fitting to the data, using lambda = 0.001
lam = 0.5
RE_fdPar = fdPar(RE_basis, Lfdobj = 4 ,lambda = lam)
all_RE_smooth = smooth.basis(time_micro, all_data_on_RE_mat, RE_fdPar) 
all_RE_fd = all_RE_smooth$fd  #extract the estimated functional data
# Evaluate the fitted Signals data at given time points
all_RE_fitted<-eval.fd(time_micro, all_RE_fd)
matplot(time_micro, all_RE_fitted, type='l', lty=1, pch='o', xlab='MicroSecond',ylab='Fitted data RE_1')
title(main = "All RE_1:-- B-spline , nbasis = 265, lambda = 0.5")


# .... All LE_1 Data:

# Create a B-spline basis with lambda penalty term " For All LE_1 Data
all_data_on_LE_mat = as.matrix(t_data_on_LE)

norder = 12
nbasis = 255 + norder - 2   # use max number of basis
LE_basis = create.bspline.basis(c(0,1499), nbasis, norder, time_micro)
#Let’s fix a penalty using a derivative of order 4 and compute the smoothed fitting to the data, using lambda = 0.001
lam = 0.5
LE_fdPar = fdPar(LE_basis, Lfdobj = 4 ,lambda = lam)
all_LE_smooth = smooth.basis(time_micro, all_data_on_LE_mat, LE_fdPar) 
all_LE_fd = all_LE_smooth$fd  #extract the estimated functional data
# Evaluate the fitted Signals data at given time points
all_LE_fitted<-eval.fd(time_micro, all_LE_fd)
matplot(time_micro, all_LE_fitted, type='l', lty=1, pch='o', xlab='MicroSecond',ylab='Fitted data LE_1')
title(main = "all LE_1:-- B-spline , nbasis = 265, lambda = 0.5")



# /////////////////////// Step 3: Data registration  /////////////////


# ................  3.1 Landmark registration  ................

# Now let’s compute the Signals acceleration, that is the second derivative of 
# the fitted curves (smoothed data), let’s compute the mean and make a plot 
# (single functional data in  green, mean in black)

# Calculate the second derivative (acceleration) of the smoothed LE_1 and RE_1 Signals data:

# ..... All RE_1 Data:
accelfdUN = deriv.fd(all_RE_fd, 2)
accel=eval.fd(time_micro,accelfdUN)
accelmeanfdUN = rowMeans(accel)
#join the mean to the data matrix for the plot
a_all_RE = cbind(accel, accelmeanfdUN) 
# Plot the acceleration data
matplot(time_micro, a_all_RE, type='l', lty=1, pch='o', xlab='MicroSecond',ylab='Acceleration/Second Derivative', col=c(rep(3,25),1), ylim=c(-1,1))
title(main = "RE_1:  Acceleration of the Estimated Curve")

# ..... All LE_1 Data:
accelfdUN = deriv.fd(all_LE_fd, 2)
accel=eval.fd(time_micro,accelfdUN)
accelmeanfdUN = rowMeans(accel)
#join the mean to the data matrix for the plot
a_all_LE = cbind(accel, accelmeanfdUN) 
# Plot the acceleration data
matplot(time_micro, a_all_LE, type='l', lty=1, pch='o', xlab='MicroSecond',ylab='Acceleration/Second Derivative', col=c(rep(3,25),1), ylim=c(-1,1))
title(main = "LE_1:  Acceleration of the Estimated Curve") 

#  ... REsult:... since no Peak and vallay in plots so we do not need Landmark registration


# /////////////////////// Step 4: PCA)  //////////////////////////


## ............................. Principal Component Analysis ................

# ..... All RE_1 Data:

# Create a B-spline basis with lambda penalty term " For All RE_1 Data
norder = 9
nbasis = 255 + norder - 2   # use max number of basis
RE_basis = create.bspline.basis(c(0,1499), nbasis, norder, time_micro)
#Let’s fix a penalty using a derivative of order 4 and compute the smoothed fitting to the data, using lambda = 0.001
lam = 0.5
RE_fdPar = fdPar(RE_basis, Lfdobj = 4 ,lambda = lam)
all_RE_smooth = smooth.basis(time_micro, all_data_on_RE_mat, RE_fdPar) 
all_RE_fd = all_RE_smooth$fd  
all_RE_fitted<-eval.fd(time_micro, all_RE_fd)

mean_fitted<-rowMeans(all_RE_fitted)
binded_matrix<-cbind(all_RE_fitted, mean_fitted)
matplot(time_micro, binded_matrix, type='l', lty=1, pch='o',
        xlab='time_microsecond',ylab='RE_1: Signals', col=c(rep(3,26),1))

# ................... continuous registration .........

# Let’s register the data using continuous registration:
# Here, a 15x24 matrix of zeros is used to initialize the functional data object Wfd0CR. 
#The matrix dimensions match the number of basis functions and the number of curves.

wbasisCR = create.bspline.basis(c(0,1499), 15, 5)
Wfd0CR = fd(matrix(0,15,24),wbasisCR) #initialize with a set of null functions
WfdParCR = fdPar(Wfd0CR, 1, 1)     # use a penalty of order 1 but with lambda=1
target=mean.fd(all_RE_smooth$fd)
# Perform continuous registration
regList = register.fd(target,all_RE_smooth$fd, WfdParCR)

regfd=eval.fd(time_micro,regList$regfd)
mean_reg = rowMeans(regfd)

# plot the registered RE_1 data and the new mean
a_RE=cbind(regfd, mean_reg)
matplot(time_micro, a_RE, type='l', lty=1, pch='o',
        xlab='time_microsecond',ylab='registered RE_1 Signal',
        main='continuous registration', col=c(rep(3,24),1))

# Get the warping function from continuous registration for first diagnosis_levels = [1]
warpfd1=regList$warpfd
warpval1=eval.fd(time_micro,warpfd1)
plot(time_micro,warpval1[,1],type = 'l', col='blue', xlab='time_MicroSecond',
     ylab='RE_Signals(MicroSecond)', main='Warping function: First diagnosis')
lines(time_micro,time_micro,type='l',col='red')

# ...........  PCA ..........

#Let's compute the principal components of the registered data, using the function
# pca.fd(fdobj, nharm)
# where
# fdobj = functional data from which we want to compute the PC's
# nharm = number of PC's to be retained
# We keep 2 PC's and see how much variance do they account for

# Perform PCA on the registered temperature data
PCout=pca.fd(regList$regfd, nharm = 2)
print(PCout$varprop)
plot.pca.fd(PCout)

# ... result : We can see that the first PC is already explaining more than 94% of the variance.

# Plot the harmonics of principal components
plot.fd(PCout$harmonics, col=c(4,2), 
        main='blue=PC1, red=PC2', xlab = 'time_MicroSecond', 
        ylab = 'Principal components')


# Get the variance-maximizing principal components
varmax_pcout=varmx.pca.fd(PCout)
plot.pca.fd(varmax_pcout)
plot.fd(varmax_pcout$harmonics, col=c(4,2), 
        main='blue=PC1, red=PC2', xlab = 'time_MicroSecond', 
        ylab = 'PCs after VARMAX')




## //////////////////// ........... Step5: Functional  regression  /////////////////////////////


# Read visual acuity data from CSV file
visual_acuity_path <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/participants_info_24_Group.csv"

visual_acuity = read.csv(visual_acuity_path)

# va_re: "Visual acuity" for the right eye, measured on logMar (logarithm of 
# the minimum angle of resolution) scale. A logMAR value of 0 denotes "normal" vision,
# while values above 0 indicate a decrease in visual acuity.Conversely,
# negative logMAR values indicate better-than-normal visual acuity.
# va_le: "Visual acuity" for the Left eye

# ..................5.1- RE_1 and va_re .....

# We will use 'va_re_logMar' as 'response variable', and the
# "RE_1" profile as "covariate'.

visual_acuity$va_re_logMar <- as.numeric(gsub("[[:space:]]", "", visual_acuity$va_re_logMar))


# Adjust the margins
par(mar = c(10, 4, 4, 2) + 0.1)

# Create a bar plot with improved text
barplot(
  visual_acuity$va_re_logMar, 
  names.arg = visual_acuity$diagnosis1, 
  col = "skyblue", 
  las = 2,
  cex.names = 0.7, # Reduce the size of the names
  main = "RE_1: Visual Acuity by Diagnosis Type", 
  ylab = "Visual acuity of RE_1"
)
# ........
mat <- matrix(visual_acuity$va_re_logMar, ncol = 1, dimnames = list(visual_acuity$diagnosis1, NULL))
visual_mat = (mat)
visual_mat_1 <- as.numeric(visual_mat)
names(visual_mat_1) = c(visual_acuity$diagnosis1) 
# ............

smallbasis  <- create.fourier.basis(c(0,1499), 5)
# Build the functional data
RE_fd <- smooth.basis(time_micro, all_data_on_RE_mat, smallbasis)
RE_fd = RE_fd$fd

## Manual construction of xfdlist and betalist
RE_list = vector("list",2)
RE_list[[1]] = rep(1,24)
RE_list[[2]] = RE_fd

conbasis = create.constant.basis(c(0,1499))
betabasis = create.fourier.basis(c(0,1499),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

# Perform functional regression
model_func_regress_visual <- fRegress(visual_mat_1 ~ RE_fd)
visual_mat_1.fit3 <- model_func_regress_visual$yhatfdobj
# Plot the fit
plot(model_func_regress_visual$betaestlist[[2]], xlab="MicroSecond",
     ylab="Beta for RE Signals", main='Regression with input setup')
# Plot the data and the fit
plot(visual_mat_1.fit3, visual_mat_1, type="p", pch=20, col = "red", xlab='Predicted y', ylab='Observed y')
lines(visual_mat_1.fit3, visual_mat_1.fit3)

# Let's assess the "quality of fit"
visual_mat_1_hat1 = model_func_regress_visual$yhatfdobj
visual_mat_1_res1 = visual_mat_1 - visual_mat_1_hat1
SSE1.1 = sum((visual_mat_1_res1)^2)
SSE0 = sum((visual_mat_1 - mean(visual_mat_1))^2)
# We can now compute the squared multiple "correlation"
# and the usual F-ratio 
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/18)  # 24 - (5+1) = 18, df1 = nbase
RSQ1
Fratio1
#Compute the (approximate) pvalue for the test
# H0: all coefficients are identically 0
pvalue=1-pf(Fratio1,5,18)
pvalue


# ..................5.2- LE_1 and va_le .....

# We will use 'va_le_logMar' as 'response variable', and the
# "LE_1" profile as "covariate'.

visual_acuity$va_le_logMar <- as.numeric(gsub("[[:space:]]", "", visual_acuity$va_le_logMar))

# Create a bar plot
barplot(visual_acuity$va_le_logMar, names.arg = visual_acuity$diagnosis1, col = "skyblue", las = 2,
        main = "Visual Acuity by Diagnosis Type", ylab = "Visual acuity of LE_1")

# Adjust the margins
par(mar = c(10, 4, 4, 2) + 0.1)

# Create a bar plot with improved text
barplot(
  visual_acuity$va_le_logMar, 
  names.arg = visual_acuity$diagnosis1, 
  col = "skyblue", 
  las = 2,
  cex.names = 0.7, # Reduce the size of the names
  main = "LE_1: Visual Acuity by Diagnosis Type", 
  ylab = "Visual acuity of LE_1"
)

# ........
mat <- matrix(visual_acuity$va_le_logMar, ncol = 1, dimnames = list(visual_acuity$diagnosis1, NULL))
visual_mat = (mat)
visual_mat_1 <- as.numeric(visual_mat)
names(visual_mat_1) = c(visual_acuity$diagnosis1) 
# ............

smallbasis  <- create.fourier.basis(c(0,1499), 5)
# Build the functional data
RE_fd <- smooth.basis(time_micro, all_data_on_RE_mat, smallbasis)
RE_fd = RE_fd$fd

## Manual construction of xfdlist and betalist
RE_list = vector("list",2)
RE_list[[1]] = rep(1,24)
RE_list[[2]] = RE_fd

conbasis = create.constant.basis(c(0,1499))
betabasis = create.fourier.basis(c(0,1499),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

# Perform functional regression
model_func_regress_visual <- fRegress(visual_mat_1 ~ RE_fd)
visual_mat_1.fit3 <- model_func_regress_visual$yhatfdobj
# Plot the fit
plot(model_func_regress_visual$betaestlist[[2]], xlab="MicroSecond",
     ylab="Beta for LE Signals", main='LE_1: Regression with input setup') 
# Plot the data and the fit
plot(visual_mat_1.fit3, visual_mat_1, type="p", pch=20, col = "red", xlab='Predicted y', ylab='Observed y', main='LE_1 Regression')
lines(visual_mat_1.fit3, visual_mat_1.fit3)

# Let's assess the "quality of fit"
visual_mat_1_hat1 = model_func_regress_visual$yhatfdobj
visual_mat_1_res1 = visual_mat_1 - visual_mat_1_hat1
SSE1.1 = sum((visual_mat_1_res1)^2)
SSE0 = sum((visual_mat_1 - mean(visual_mat_1))^2)
# We can now compute the squared multiple "correlation"
# and the usual F-ratio 
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/18)  # 24 - (5+1) = 18, df1 = nbase
RSQ1
Fratio1
#Compute the (approximate) pvalue for the test
# H0: all coefficients are identically 0
pvalue=1-pf(Fratio1,5,18)
pvalue




# ................... Multiple Regression:  RE_fd + LE_fd ...................

# .... LE_1  ...

# ..................  5.3- RE_fd + LE_fd and'va_le_logMar' .....

# We will use 'va_le_logMar' as 'response variable', and the
# "LE_1" profile as "covariate'.

visual_acuity$va_le_logMar <- as.numeric(gsub("[[:space:]]", "", visual_acuity$va_le_logMar))

# Create a bar plot
barplot(visual_acuity$va_le_logMar, names.arg = visual_acuity$diagnosis1, col = "skyblue", las = 2,
        main = "Visual Acuity by Diagnosis Type", ylab = "Visual acuity of LE_1")

# ........
mat <- matrix(visual_acuity$va_le_logMar, ncol = 1, dimnames = list(visual_acuity$diagnosis1, NULL))
visual_mat = (mat)
visual_mat_1 <- as.numeric(visual_mat)
names(visual_mat_1) = c(visual_acuity$diagnosis1) 
# ............

smallbasis  <- create.fourier.basis(c(0,1499), 5)
# Build the functional data
RE_fd <- smooth.basis(time_micro, all_data_on_RE_mat, smallbasis)
RE_fd = RE_fd$fd

LE_fd <- smooth.basis(time_micro, all_data_on_LE_mat, smallbasis)
LE_fd = LE_fd$fd

## Manual construction of xfdlist and betalist
RE_list = vector("list",2)
RE_list[[1]] = rep(1,24)
RE_list[[2]] = RE_fd

LE_list = vector("list",2)
LE_list[[1]] = rep(1,24)
LE_list[[2]] = LE_fd

conbasis = create.constant.basis(c(0,1499))
betabasis = create.fourier.basis(c(0,1499),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

# Perform Multiple functional regression
model_func_regress_visual <- fRegress(visual_mat_1 ~ RE_fd + LE_fd)
visual_mat_1.fit3 <- model_func_regress_visual$yhatfdobj
# Plot the fit
plot(model_func_regress_visual$betaestlist[[2]], xlab="MicroSecond",
     ylab="Beta for LE Signals", main='Regression with input setup') 
# Plot the data and the fit
plot(visual_mat_1.fit3, visual_mat_1, type="p", pch=20, col = "red", xlab='Predicted y', ylab='Observed y', main='RE_fd + LE_fd and "va_le": Regression')
lines(visual_mat_1.fit3, visual_mat_1.fit3)

# Let's assess the "quality of fit"
visual_mat_1_hat1 = model_func_regress_visual$yhatfdobj
visual_mat_1_res1 = visual_mat_1 - visual_mat_1_hat1
SSE1.1 = sum((visual_mat_1_res1)^2)
SSE0 = sum((visual_mat_1 - mean(visual_mat_1))^2)
# We can now compute the squared multiple "correlation"
# and the usual F-ratio 
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/18)  # 24 - (5+1) = 18, df1 = nbase
RSQ1
Fratio1
#Compute the (approximate) pvalue for the test
# H0: all coefficients are identically 0
pvalue=1-pf(Fratio1,5,18)
pvalue



# ---------- RE_1

# ..................5.4- RE_fd + LE_fd and  'va_re_logMar' .....

# We will use 'va_le_logMar' as 'response variable', and the
# "LE_1" profile as "covariate'.

visual_acuity$va_re_logMar <- as.numeric(gsub("[[:space:]]", "", visual_acuity$va_re_logMar))

# Create a bar plot
barplot(visual_acuity$va_re_logMar, names.arg = visual_acuity$diagnosis1, col = "skyblue", las = 2,
        main = "Visual Acuity by Diagnosis Type", ylab = "Visual acuity of LE_1")

# ........
mat <- matrix(visual_acuity$va_re_logMar, ncol = 1, dimnames = list(visual_acuity$diagnosis1, NULL))
visual_mat = (mat)
visual_mat_1 <- as.numeric(visual_mat)
names(visual_mat_1) = c(visual_acuity$diagnosis1) 
# ............

smallbasis  <- create.fourier.basis(c(0,1499), 5)
# Build the functional data
RE_fd <- smooth.basis(time_micro, all_data_on_RE_mat, smallbasis)
RE_fd = RE_fd$fd

LE_fd <- smooth.basis(time_micro, all_data_on_LE_mat, smallbasis)
LE_fd = LE_fd$fd

## Manual construction of xfdlist and betalist
RE_list = vector("list",2)
RE_list[[1]] = rep(1,24)
RE_list[[2]] = RE_fd

LE_list = vector("list",2)
LE_list[[1]] = rep(1,24)
LE_list[[2]] = LE_fd

conbasis = create.constant.basis(c(0,1499))
betabasis = create.fourier.basis(c(0,1499),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

# Perform functional regression
model_func_regress_visual <- fRegress(visual_mat_1 ~ RE_fd + LE_fd)
visual_mat_1.fit3 <- model_func_regress_visual$yhatfdobj
# Plot the fit
plot(model_func_regress_visual$betaestlist[[2]], xlab="MicroSecond",
     ylab="Beta for LE Signals", main='Regression with input setup') 
# Plot the data and the fit
plot(visual_mat_1.fit3, visual_mat_1, type="p", pch=20, col = "red", xlab='Predicted y', ylab='Observed y', main='RE_fd + LE_fd and "va_re": Regression')
lines(visual_mat_1.fit3, visual_mat_1.fit3)

# Let's assess the "quality of fit"
visual_mat_1_hat1 = model_func_regress_visual$yhatfdobj
visual_mat_1_res1 = visual_mat_1 - visual_mat_1_hat1
SSE1.1 = sum((visual_mat_1_res1)^2)
SSE0 = sum((visual_mat_1 - mean(visual_mat_1))^2)
# We can now compute the squared multiple "correlation"
# and the usual F-ratio 
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/18)  # 24 - (5+1) = 18, df1 = nbase
RSQ1
Fratio1
#Compute the (approximate) pvalue for the test
# H0: all coefficients are identically 0
pvalue=1-pf(Fratio1,5,18)
pvalue



# ........................5-5: categorical response variable: ...........

# install.packages("fda.usc", repos = "http://cran.us.r-project.org")

library(fda.usc)
library(nnet)


# .......................... with PCA

# Define the path to the PERG data and participants info

visual_acuity_path_R <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/va_acuity_RE_1.csv"
visual_acuity_path_L <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/va_acuity_LE_1.csv"
visual_acuity_R1 = read.csv(visual_acuity_path_R)
visual_acuity_L1 = read.csv(visual_acuity_path_L)

# remove first col
visual_acuity_R1 <- visual_acuity_R1[,-1] 
visual_acuity_L1 <- visual_acuity_L1[,-1] 

RE_data_path <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/combined_RE_1.csv"
LE_data_path <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/combined_LE_1.csv"
data_frame_RE <- read.csv(RE_data_path)
data_frame_LE <- read.csv(LE_data_path)

# Extract functional covariates for the new observation
va_re_logMar_1 <- as.matrix(visual_acuity_R1)
va_le_logMar_1 <- as.matrix(visual_acuity_L1)
RE_1 <- as.matrix(data_frame_RE)
LE_1 <- as.matrix(data_frame_LE)

time_points <- time_micro

# Define a B-spline basis
basis <- create.bspline.basis(rangeval = c(0,1499), nbasis = 10)  


# Smooth the functional covariates for all subjects
va_re_fd <- smooth.basis(time_points, va_re_logMar_1, basis)$fd
va_le_fd <- smooth.basis(time_points, va_le_logMar_1, basis)$fd
RE_fd <- smooth.basis(time_points, RE_1, basis)$fd
LE_fd <- smooth.basis(time_points, LE_1, basis)$fd

# Perform functional PCA on each functional covariate
va_re_pca <- pca.fd(va_re_fd, nharm = 3)  # Adjust nharm as needed
va_le_pca <- pca.fd(va_le_fd, nharm = 3)  # Adjust nharm as needed
RE_pca <- pca.fd(RE_fd, nharm = 3)  # Adjust nharm as needed
LE_pca <- pca.fd(LE_fd, nharm = 3)  # Adjust nharm as needed


# Extract the scores (principal components) for each functional covariate
va_re_scores <- as.data.frame(va_re_pca$scores)
va_le_scores <- as.data.frame(va_le_pca$scores)
RE_scores <- as.data.frame(RE_pca$scores)
LE_scores <- as.data.frame(LE_pca$scores)


# Rename columns to ensure unique names
colnames(va_re_scores) <- paste0("va_re_PC", 1:ncol(va_re_scores))
colnames(va_le_scores) <- paste0("va_le_PC", 1:ncol(va_le_scores))
colnames(RE_scores) <- paste0("RE_PC", 1:ncol(RE_scores))
colnames(LE_scores) <- paste0("LE_PC", 1:ncol(LE_scores))

# Combine the principal components into a single data frame
combined_data <- cbind(va_re_scores, va_le_scores, RE_scores, LE_scores)

# Extract response variable (categorical)
diagnosis_levels <- c(    'Paracentral acute middle maculopathy type 2',
                          'Bilateral optic nerve atrophy',
                          'Autoimmune retinopathy',
                          'Foveal hypoplasia',
                          'Inherited optic atrophy',
                          'Blue cone monochromatism',
                          'Chorioretinopathy Birdshot type',
                          'AlbinoidismA',
                          'Infectious neuritis',
                          'Orbital ischemia',
                          'Optic neuropathy',
                          'Traumatic optic neuropathy',
                          'Sarcoidosis neuropathy',
                          'Cone-Rod dystrophy',
                          'Arterio-venous malformation in right thalamus',
                          'Albinism',
                          'Central nervous system disorder',
                          'Bilateral Optic Atrophy',
                          'External ophthalmoplegia',
                          'Congenital Achromatopsia',
                          'Macular dystrophy',
                          'Anisometropic amblyopia',
                          'Acute macular neuroretinopathy',
                          'Normal')
diagnosis_levels <- as.factor(diagnosis_levels)
diagnosis_levels

# Fit the multinomial logistic regression model
multinom_model <- multinom(diagnosis_levels ~ ., data = combined_data)

# Summary of the model
summary(multinom_model)

# Predicted values
predicted_values <- predict(multinom_model, newdata = combined_data)
predicted_values
# Create a confusion matrix
confusion_mtx <- confusionMatrix(predicted_values, diagnosis_levels)

# Print the confusion matrix
print(confusion_mtx)


# Extract evaluation metrics, handling cases where precision, recall are NaN
accuracy <- confusion_mtx$overall['Accuracy']

precision <- confusion_mtx$byClass['Pos Pred Value']
recall <- confusion_mtx$byClass['Sensitivity']
f1_score <- 2 * (precision * recall) / (precision + recall)

# Handle NA values in metrics
precision[is.na(precision)] <- 0
recall[is.na(recall)] <- 0
f1_score[is.na(f1_score)] <- 0

# Print evaluation metrics
cat("Accuracy:", accuracy, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-Score:", f1_score, "\n")

# Visualize the confusion matrix using a heatmap
confusion_mtx_df <- as.data.frame(confusion_mtx$table)
colnames(confusion_mtx_df) <- c("Reference", "Prediction", "Freq")
confusion_mtx_melt <- melt(confusion_mtx_df, id.vars = c("Reference", "Prediction"))

ggplot(data = confusion_mtx_melt, aes(x = Reference, y = Prediction, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Confusion Matrix", x = "Actual Diagnosis", y = "Predicted Diagnosis") +
  theme_minimal()




# ............. TEST   ..............

#  the path to the PERG data and participants info

visual_acuity_path_R_test <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/va_acuity_RE_1_test24.csv"
visual_acuity_path_L_test <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/va_acuity_LE_1_test24.csv"
visual_acuity_R1 = read.csv(visual_acuity_path_R_test)
visual_acuity_L1 = read.csv(visual_acuity_path_L_test)

# remove first col
visual_acuity_R1 <- visual_acuity_R1[,-1] 
visual_acuity_L1 <- visual_acuity_L1[,-1] 

RE_data_path_test <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/combined_RE_1_test24.csv"
LE_data_path_test <- "C:/00 Transferable Folders/00 Data Science/04. My Projects/GitHub Repository/08. Functional Data Analysis Project/Dataset/Merged_data/combined_LE_1_test24.csv"
data_frame_RE <- read.csv(RE_data_path_test)
data_frame_LE <- read.csv(LE_data_path_test)

# Extract functional covariates for the new observation
new_va_re_logMar_1 <- as.matrix(visual_acuity_R1)
new_va_le_logMar_1 <- as.matrix(visual_acuity_L1)
new_RE_1 <- as.matrix(data_frame_RE)
new_LE_1 <- as.matrix(data_frame_LE)

# Extract microseconds (time points) for the new observation
new_time_points <- time_micro

# Ensure the dimensions match
if(nrow(new_va_re_logMar_1) != length(new_time_points)) {
  new_va_re_logMar_1 <- t(new_va_re_logMar_1)
}
if(nrow(new_va_le_logMar_1) != length(new_time_points)) {
  new_va_le_logMar_1 <- t(new_va_le_logMar_1)
}
if(nrow(new_RE_1) != length(new_time_points)) {
  new_RE_1 <- t(new_RE_1)
}
if(nrow(new_LE_1) != length(new_time_points)) {
  new_LE_1 <- t(new_LE_1)
}

# Define the same B-spline basis used in the training data
basis <- create.bspline.basis(rangeval = c(min(new_time_points), max(new_time_points)), nbasis = 10)  # Adjust nbasis as needed

# Smooth the functional covariates for the new observation
new_va_re_fd <- smooth.basis(new_time_points, new_va_re_logMar_1, basis)$fd
new_va_le_fd <- smooth.basis(new_time_points, new_va_le_logMar_1, basis)$fd
new_RE_fd <- smooth.basis(new_time_points, new_RE_1, basis)$fd
new_LE_fd <- smooth.basis(new_time_points, new_LE_1, basis)$fd

# Project New Functional Data onto Principal Components
# Manually project the new functional data onto the principal components obtained from the PCA:

# Function to project new functional data onto existing PCA components
project_to_pca <- function(new_fd, pca_fd) {
  scores <- inprod(new_fd, pca_fd$harmonics)
  return(scores)
}

# Project new functional data onto the principal components
new_va_re_scores <- project_to_pca(new_va_re_fd, va_re_pca)
new_va_le_scores <- project_to_pca(new_va_le_fd, va_le_pca)
new_RE_scores <- project_to_pca(new_RE_fd, RE_pca)
new_LE_scores <- project_to_pca(new_LE_fd, LE_pca)

# Convert the scores to data frames
new_va_re_scores_df <- as.data.frame(new_va_re_scores)
new_va_le_scores_df <- as.data.frame(new_va_le_scores)
new_RE_scores_df <- as.data.frame(new_RE_scores)
new_LE_scores_df <- as.data.frame(new_LE_scores)

# Rename columns to ensure unique names
colnames(new_va_re_scores_df) <- paste0("va_re_PC", 1:ncol(new_va_re_scores_df))
colnames(new_va_le_scores_df) <- paste0("va_le_PC", 1:ncol(new_va_le_scores_df))
colnames(new_RE_scores_df) <- paste0("RE_PC", 1:ncol(new_RE_scores_df))
colnames(new_LE_scores_df) <- paste0("LE_PC", 1:ncol(new_LE_scores_df))

# Combine the principal components into a single data frame
new_combined_data <- cbind(new_va_re_scores_df, new_va_le_scores_df, new_RE_scores_df, new_LE_scores_df)

# Predict the diagnosis for the new observation
new_predicted_diagnosis <- predict(multinom_model, newdata = new_combined_data)

# Print the predicted diagnosis
print(new_predicted_diagnosis)

test_diagnosis_levels <- c(
  'Normal', 'Macular dystrophy', 'Macular dystrophy', 
  'Congenital achromatopsia', 'Congenital achromatopsia', 
  'Bilateral optic nerve atrophy', 'Autoimmune retinopathy', 
  'Albinism', 'Orbital ischemia', 'Autoimmune retinopathy', 
  'Normal', 'Normal', 'Normal', 'Normal', 'Macular dystrophy', 'Macular dystrophy', 
  'Congenital achromatopsia', 'Congenital achromatopsia', 
  'Bilateral optic nerve atrophy', 'Autoimmune retinopathy', 
  'Albinism', 'Orbital ischemia', 'Autoimmune retinopathy', 
  'Normal'
)



# /////////// ...................Step 6: Depth measures   ...................

# RE_1

# extract the RE signals and transform them in functional data
daytime <- time_micro + 0.5
time_basis15 <- create.fourier.basis(c(0,1499), 15)
smoothList <- with(data_RE, smooth.basis(daytime,
                                         all_data_on_RE_mat,
                                         time_basis15, fdnames=list("Microsecond", "Gender", "RE Signal")))

RE_1_fd<-smoothList$fd
plot(RE_1_fd)
b1<-boxplot(RE_1_fd, method = "MBD")
b2<-boxplot(RE_1_fd, method = "BD2")
b3<-boxplot(RE_1_fd, method = "Both")
b1$medcurve
b2$medcurve
b3$medcurve
DM<-b1$depth
# The Wilcoxon test in this case is testing
# H0: the means in the two groups are equal
mat_gender <- matrix(visual_acuity$sex, ncol = 1, dimnames = list(visual_acuity$id_record, NULL))
gender = mat_gender
index_gender<-which(gender=='Male')
D_Male<-DM[index_gender]
D_Female<-DM[-index_gender]
wilcox.test(D_Male,D_Female)
# ..........Result: The pvalue is quite small, meaning that there is a NO significant difference between the 
# Male and Female RE Signals.


# LE_1

# extract the RE signals and transform them in functional data
daytime <- time_micro + 0.5
time_basis15 <- create.fourier.basis(c(0,1499), 15)
smoothList <- with(data_RE, smooth.basis(daytime,
                                         all_data_on_LE_mat,
                                         time_basis15, fdnames=list("Microsecond", "Gender", "RE Signal")))

RE_1_fd<-smoothList$fd
plot(RE_1_fd)
b1<-boxplot(RE_1_fd, method = "MBD")
b2<-boxplot(RE_1_fd, method = "BD2")
b3<-boxplot(RE_1_fd, method = "Both")
b1$medcurve
b2$medcurve
b3$medcurve
DM<-b1$depth
# The Wilcoxon test in this case is testing
# H0: the means in the two groups are equal
mat_gender <- matrix(visual_acuity$sex, ncol = 1, dimnames = list(visual_acuity$id_record, NULL))
gender = mat_gender
index_gender<-which(gender=='Male')
D_Male<-DM[index_gender]
D_Female<-DM[-index_gender]
wilcox.test(D_Male,D_Female)
# ..........Result: The pvalue is quite small, meaning that there is a NO significant difference between the 
# Male and Female RE Signals.


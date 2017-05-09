## can we find that a composite likelihood with more blocks gives a higher likelihood than one with fewer blocks?
#
library(MASS)
#1. Exponential covariance matrix:

exp.param <- c(10, 1.2)
lat.data <- seq(10,80, by = 5) 
location.data <- expand.grid(lat.data, lat.data)
dist.mat <- geoDistance(location.data, location.data)
cov.mat <- expCov(dist.mat, exp.param, 1) #no nugget:
#
#2. simulate data:
n <- dim(location.data)[1]
mu <- array(0, dim=c(n,1))
out.data <- mvrnorm(1,mu,cov.mat)
#
# 3. calculate the likelihood under true:
m.error <- 0*diag(n) #identity matrix - measurement error
test.data <- cbind(location.data, out.data)

start1 <- Sys.time()
comp.NegLogLik <- composite.nll(exp.param, test.data,c(2,1),"Exponential", m.error, "NegLogLik") 
#nudge <- 0.00001
#adj.Cov.mat <- Cov.mat + nudge*diag(n_ind)
real.NegLogLik <- 1/2*log(det(cov.mat)) + 1/2*t(out.data)%*%solve(cov.mat)%*%out.data
real2.NegLogLik <- covariance.neg.lik.calc(exp.param, out.data, n, m.error, dist.mat, "Exponential")
#log-likelihood values agree 11/04/17
end1 <- Sys.time()
time1 <- end1 - start1

start2 <- Sys.time()
comp1.NegLogLik <- composite.nll(test.data,c(4,4),"Exponential",wrong.param, m.error, "NegLogLik") 
#Agrees with real.NegLogLik 
end2 <- Sys.time()
time2 <- end2 - start2


################################################################
################################################################
#
# Parameter estimation with composite likelihood function:
#
#1. Parameters:
exp.param <- c(0.5, 0.02)
lat.data <- seq(10,20, by = 1) 
location.data <- expand.grid(lat.data, lat.data)
dist.mat <- geoDistance(location.data, location.data)
cov.mat <- expCov(dist.mat, exp.param, 1) #no nugget:
#
#2. simulate data:
n <- dim(location.data)[1]
mu <- array(0, dim=c(n,1))
out.data <- mvrnorm(1,mu,cov.mat)
param.lower <- c(0.1,0.001)
param.upper <- c(1,0.1)
m.error <- 0*diag(n) #identity matrix - measurement error
test.data <- cbind(location.data, out.data)
#2.a)
#check that the output is the same:
full.loglik <- covariance.neg.lik.calc(cov.param = exp.param, data = out.data, n=n, m.error = m.error, dist.mat = dist.mat, model = "Exponential")
comp.loglik <- composite.nll(data = test.data, m=c(2,1),cov.model = "Exponential", Cov.par = exp.param, m.error = m.error, scoreFunc = "NegLogLik")
#
#3. MLE:
mle.output.full <- optim(exp.param, covariance.neg.lik.calc, data = out.data, n = n, m.error =  m.error, dist.mat = dist.mat, model = "Exponential", lower = param.lower, 
                               upper = param.upper, method="L-BFGS-B")
mle.output.full

mle.output.full <- optim(c(10,3), covariance.neg.lik.calc, data = out.data, n = n, m.error =  m.error, dist.mat = dist.mat, model = "Exponential", lower = param.lower, 
                         upper = param.upper, method="L-BFGS-B")
mle.output.full
#
#4. MLE from composite likelihood:
mle.output.partial <- optim(exp.param, composite.nll, data = test.data, m = c(2,1), m.error =  m.error, cov.model = "Exponential", scoreFunc = "NegLogLik", lower = param.lower, 
                         upper = param.upper, method="L-BFGS-B")
mle.output.partial

mle.output.full.NM <- optim(exp.param, covariance.neg.lik.calc, data = out.data, n = n, m.error =  m.error, dist.mat = dist.mat, model = "Exponential", method="Nelder-Mead")
mle.output.full.NM
mle.output.NM <- optim(exp.param, composite.nll, data = test.data, m = c(2,1), m.error =  m.error, cov.model = "Exponential", scoreFunc = "NegLogLik", method="Nelder-Mead")
mle.output.NM


####################################################################################
####################################################################################
#
# 1. Testing the likelihood functions:
#
exp.param <- c(10, 1.2)
lat.data <- seq(10,80, by = 10) 
location.data <- expand.grid(lat.data, lat.data)
dist.mat <- geoDistance(location.data, location.data)
cov.mat <- expCov(dist.mat, exp.param, 0) #no nugget:
#
# simulate data:
n <- dim(location.data)[1]
mu <- array(0, dim=c(n,1))
out.data <- mvrnorm(1,mu,cov.mat)
length(out.data)
m.error <- 0*diag(n) #identity matrix - measurement error
test.data <- cbind(location.data, out.data)
dim(test.data)
#
# calculate the likelihood under true:
Cov.mat <- expCov(dist.mat,exp.param, 0)  
#neg. log-likelihood:
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(test.data[,3])%*%solve(Cov.mat)%*%test.data[,3]
# Calculate Full Likelihood using composite code with 1x2 blocks:
full.NegLogLik <- covariance.neg.lik.calc(exp.param, out.data, 0, dist.mat, "Exponential")
# Calculate Composite Likelihood:
composite.NegLogLik <- composite.nll(exp.param, test.data, c(2,1), "Exponential", 0, "NegLogLik")  

#
#
# 2. Matern covariance:
mat.param <- c(0.4, 10, 1.2)
lat.data <- seq(10,80, by = 10) 
location.data <- expand.grid(lat.data, lat.data)
dist.mat <- geoDistance(location.data, location.data)
cov.mat <- maternCov(dist.mat, mat.param, 0) #no nugget:
#
# simulate data:
n <- dim(location.data)[1]
mu <- array(0, dim=c(n,1))
out.data <- mvrnorm(1,mu,cov.mat)
length(out.data)
m.error <- 0*diag(n) #identity matrix - measurement error
test.data <- cbind(location.data, out.data)
dim(test.data)
#
# calculate the likelihood under true:
Cov.mat <- maternCov(dist.mat,mat.param, 0)  
#neg. log-likelihood:
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(test.data[,3])%*%solve(Cov.mat)%*%test.data[,3]
# Calculate Full Likelihood using composite code with 1x2 blocks:
full.NegLogLik <- covariance.neg.lik.calc(mat.param, out.data, 0, dist.mat, "Matern")
# Calculate Composite Likelihood:
composite.NegLogLik <- composite.nll(mat.param, test.data, c(2,1), "Matern", 0, "NegLogLik")  

#
#
# 3. Matern covariance w/ nugget effect:
nug.param <- c(0.4, 10, 1.2, 0.1)
lat.data <- seq(10,80, by = 10) 
location.data <- expand.grid(lat.data, lat.data)
n <- dim(location.data)[1]
m.vec.error <- array(1, dim=c(n,1))
dist.mat <- geoDistance(location.data, location.data)
cov.mat <- maternCov(dist.mat, nug.param, m.vec.error) #no nugget:
#
# simulate data:
mu <- array(0, dim=c(n,1))
out.data <- mvrnorm(1,mu,cov.mat)
length(out.data)
test.data <- cbind(location.data, out.data)
dim(test.data)
#
# calculate the likelihood under true:
Cov.mat <- maternCov(dist.mat,nug.param, m.vec.error)  
#neg. log-likelihood:
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(test.data[,3])%*%solve(Cov.mat)%*%test.data[,3]
# Calculate Full Likelihood using composite code with 1x2 blocks:
full.NegLogLik <- covariance.neg.lik.calc(nug.param, out.data, m.vec.error, dist.mat, "Matern")
# Calculate Composite Likelihood:
 #have to give the vector of measurement errors to the composite likelihood.
composite.NegLogLik <- composite.nll(nug.param, test.data, c(2,1), "Matern", m.vec.error, "NegLogLik")  

##########################################################################################
##########################################################################################
#
# Testing likelihoods for geo. observations:
maternParam <- mle.output.matern.nug$par
Cov.mat <- maternCov(dist.mat,maternParam, meas.error)  
#neg. log-likelihood:
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(residual_Obs)%*%solve(Cov.mat)%*%residual_Obs
# Calculate Full Likelihood using composite code with 1x2 blocks:
full.NegLogLik <- covariance.neg.lik.calc(maternParam, residual_Obs, meas.error, dist.mat, "Matern")
# Calculate Composite Likelihood:
#have to give the vector of measurement errors to the composite likelihood.
res.data <- cbind(locations, residual_Obs)
composite.NegLogLik <- composite.nll(maternParam, res.data, c(2,1), "Matern", meas.error, "NegLogLik") 
#
# These agree as well.


#############################################################################################
#############################################################################################
#
# Testing composition blocks:

n_ind <- 3
index <- c(500,4200,5400,1000,15000,8500,900)
j <- 2
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[j],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=T) 

## calculating the composite likelihood:
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]

location.GCM <- Model_df.full[index,2:3]
data.GCM <- Model_df.full[index,5]
## residual:
rbf.data.GCM <- cbind(location.GCM, data.GCM)
colnames(rbf.data.GCM) <- c("lon","lat","val")
rbf.fit.GCM <- gausspr(val ~ lon + lat, data = rbf.data.GCM, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
#create residuals (Obs - Predicted):
fitted.val.GCM <- predict(rbf.fit.GCM, rbf.data.GCM[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
residual_GCM <- data.GCM - fitted.val.GCM #residuals from fitted Gaussian RBF model.
residual.GCM.data <- cbind(location.GCM, residual_GCM)

m.error.GCM.vec <- array(1, dim=c(n_ind,1)) #identity matrix - measurement error

comp.NegLogLik <- composite.nll(maternParam,residual.GCM.data,c(2,1),"Matern", m.error.GCM.vec, "NegLogLik") 
comp1.NegLogLik <- composite.nll(maternParam,residual.GCM.data,c(3,1),"Matern", m.error.GCM.vec, "NegLogLik") 








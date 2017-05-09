## Code accompanying the tex file: FittingNonZeroTrendGPs_sphere (Mar. 2017)

source("compNLL_Func_v2.R")
source("Code/Init.R")
library(Imap, kernlab, gstat, sp, geosphere, sphereplot)
set.seed(1)

#load packages: gstat, sp, flexclust, Imap, kernlab, geosphere and sphereplot

# User-defined functions:
#We need to define a Kernel that uses geodesic distances not Euclidean distances:
geoGaussianKernel=function(x,y){
  #calculate the geodesic distance between x and y
  x <- as.matrix(t(x))
  y <- as.matrix(t(y)) #these have to be passed in as matrices for the dimension call to work.
  r <- geoDistance(x,y)
  return(exp(-(r^2))) #ignore kernel parameter for now
}
class(geoGaussianKernel)="kernel"


# 1. Load the data:
# Observations
# -------------
# Observations need to have points in a list of 4 columns with headers "x,y,z,std" = lon, lat, observation, uncertainty
# Observations are stored in directory Observation_Data/<Obs_set> with filenames <obs>.txt
Obs_set='P3+_SST_anom'
obs = c('lambda_100', 'lambda_90', 'lambda_80', 'lambda_70', 'lambda_60', 'lambda_50', 'lambda_40', 'lambda_30', 'lambda_20', 'lambda_10')
j <- 1 #one time instance. 
# READ OBSERVATIONS
Obs <- read.table(paste0("Observation_data/",Obs_set,"/",obs[j],".txt"),header=T)

#Plot the observation locations on the sphere:
#sphereplot package:
n <- dim(Obs)[1]
rgl.sphgrid(radius = 1)
rgl.sphpoints(cbind(Obs[,1:2],array(1,dim=c(n,1))),deg=TRUE,col='black',cex=2)
rgl.postscript("locationSphere2.eps","eps")

#create a spatial points data frame object:
locations <- cbind(Obs$x,Obs$y)
z.dataframe <- data.frame(Obs$z)
y_obs.sp = SpatialPointsDataFrame(locations, z.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(y_obs.sp) ~ lon + lat  #coordinate data
summary(y_obs.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 


# 2. Estimate the Trend:
#
# i) Gaussian RBF (using geodesic distance):
rbf.data <- cbind(locations, Obs$z)
colnames(rbf.data) <- c("lon","lat","val")
rbf.fit.geo <- gausspr(val ~ lon + lat, data = rbf.data, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
alpha(rbf.fit.geo) #model parameters
#
#create residuals (Obs - Predicted):
fitted.val <- predict(rbf.fit.geo, rbf.data[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
residual_Obs <- Obs$z - fitted.val #residuals from fitted Gaussian RBF model.

#Create spatial data frame object:
residual.dataframe <- data.frame(residual_Obs)
residual.sp = SpatialPointsDataFrame(locations, residual.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp) ~ lon + lat  #coordinate data
summary(residual.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 

# ii) Linear Model:
lm.data <- cbind(locations, Obs$z)
colnames(lm.data) <- c("lon","lat","val")
lm.data <- data.frame(lm.data)
#
#Coefficient estimates:
lm.GPtrend <- lm(val ~ lon + lat, lm.data)
residual_Obs <- lm.GPtrend$residuals
lm.coef <- lm.GPtrend$coefficients

residual.dataframe <- data.frame(residual_Obs)
residual.sp = SpatialPointsDataFrame(locations, residual.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp) ~ lon + lat  #coordinate data
summary(residual.sp)

# 3. Estimating the Parameters in the Covariance Function (Using MLE):
#First, plot the empirical variogram:
vgm.sample <- variogram(residual_Obs ~ 1, data = residual.sp, cutoff = 5000, width = 500 )  
plot(vgm.sample)
dev.copy(pdf,'sampleVGM_residual_lm.pdf')
dev.off()

#MLE for covariance marix:
#
#Matern nugget:
meas.error <- Obs[,4]^2 
n <- length(residual_Obs)
dist.mat <- geoDistance(locations,locations)
param.lower.nug <- c(0.05,0.01,0.01,0.00001) #(nu,rho,sigma,nugget)
#needs upper limit for nu <= 1/2 too.
#needs upper limit for nu <= 1/2 too.
param.upper.nug <- c(1/2,5000,10,10) #for validity on the sphere.
ini.nug <- c(0.5,2500,1,0.01)

## Matern no nugget:
param.lower <- c(0.2,0.01,0.01) #(nu,rho,sigma)
#needs upper limit for nu <= 1/2 too.
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,1000) #for validity on the sphere.
ini <- c(0.5,2500,1)
##

## Exponential : nugget:
param.lower.exp.nug <- c(0.01,0.01,0.00001) #(rho,sigma,nugget)
param.upper.exp.nug <- c(10000,10,100) 
ini.exp.nug <- c(2000,2,0.1)
##

## Exponential : no nugget:
param.lower.exp <- c(0.01,0.01) #(rho,sigma)
param.upper.exp <- c(10000,1000) #for validity on the sphere.
ini.exp <- c(15,1)
##

mle.output.matern.nug <- optim(ini.nug, covariance.neg.lik.calc, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern", lower = param.lower.nug, 
                           upper = param.upper.nug, method="L-BFGS-B")
mle.output.matern.nug
mle.output.matern <- optim(ini, covariance.neg.lik.calc, data = residual_Obs, n = n, dist.mat = dist.mat, model = "Matern", lower = param.lower, 
                               upper = param.upper, method="L-BFGS-B")
mle.output.matern

mle.output.exp.nug <- optim(ini.exp.nug, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential", lower = param.lower.exp.nug, 
                        upper = param.upper.exp.nug, method="L-BFGS-B")
mle.output.exp.nug
mle.output.exp <- optim(ini, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential", lower = param.lower, 
                           upper = param.upper, method="L-BFGS-B")
mle.output.exp


mle.output.spherical <- optim(ini, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Spherical", lower = param.lower,
                              upper = param.upper, method="L-BFGS-B")
mle.output.spherical

#Plot the fitted model onto the empirical data points:
#Exp assumes that v=1/2
model.choice.mat <- vgm(psill =  mle.output.matern$par[2], "Exp", range = mle.output.matern$par[3], nugget = mle.output.matern$par[4]) 
model.choice.sph <- vgm(psill =  mle.output.spherical$par[2], "Sph", range = mle.output.spherical$par[3], nugget = mle.output.spherical$par[4]) 


fit.vgm.exp <- fit.variogram(vgm.sample, model = model.choice.sph)
plot(vgm.sample, model = fit.vgm.exp)
dev.copy(pdf,'fittedVGM_mle_lm.pdf')
dev.off()

################################################################################################
# 4. Composite likelihood calculation: Residual Observations
test.length <- round(runif(10,1,95))
pred.data <- cbind(locations[test.length,], t(t(residual_Obs[test.length])))
#names(pred.data) <- c("lat","long","val")
#Note we will get a singular covariance matrix (Cov.mat) when we include the one location that yields a unique predicted value
#because of linear dependence. 

maternParam <- mle.output.matern$par
dist.mat <- geoDistance(pred.data[,1:2],pred.data[,1:2]) #calculate the geodesic distance
Cov.mat <- maternCov(dist.mat,maternParam)  #using estimated Matern covariance parameters
#full negative log-likelihood:
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(pred.data[,3])%*%solve(Cov.mat)%*%pred.data[,3] 

full.NegLogLik <- composite.nll(pred.data,c(2,1),"Matern",maternParam) #Agrees with real.NegLogLik at least. 

# Calculate Composite Likelihood of Nieghbouring blocks:
ptm <- proc.time() #time the operation to compare with the full likelihood function.
composite.NegLogLik <- composite.nll(pred.data, c(3,2), "Matern", maternParam)  
proc.time() - ptm


# 5. Composite Likelihood : GCM Data

# Model settings
# --------------
# Models - gridded datasets of data values as a vector
# Model outputs are stored in directory Model_Data/<Model_set> with filenames <models>.txt
# Directory should also contain a file mask.txt with the indices of the region to be analysed
# IDL file um2krig will prepare UM data appropriately
Model_set='CO2_anom'
models = c('tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti')

# Specify grid dimensions for model data
# (lon_min,lon_max,lon_int,lat_min,lat_max,lat_int)
model_grid <- c(-180,178.75,1.25,-89.375,89.375,1.25) # HadCM3 ocean temps, shifted 1 hemisphere
# model_grid <- c(-180,176.25,3.75,-90,90,2.5)          # HadCM3 atmosphere temps, shifted 1 hemisphere
# model_grid <- c(-179.5,179.5,1.0,89.5,-89.5,-1.0)     # HadISST data
#load Read_model_data (coordinates still in degrees for consistency).
i <- 3
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[i],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=T) 
dev.copy(pdf,'GCM1_plot.pdf')
dev.off()

#
#
## calculating the composite likelihood:
maternParam <- mle.output.matern.nug$par
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
#Model_df.edit <- Model_df.full[1:25000,]
####
#
# Choosing points:
#1. random sample:
#n_ind <- 1000 #we don't want to compute this over the entire grid, at this point.
#index <- round(runif(n_ind,1,dim(Model_df.full)[1]))
#
#2. Regular sample:
n_ind <- 100
jump <- 270 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index <- seq(start, 27000, by = jump)
#plot the locations on the sphere:
#sphereplot package:
rgl.sphgrid(radius = 1)
rgl.sphpoints(cbind(Model_df.full[index,2:3],array(1,dim=c(length(index),1))),deg=TRUE,col='black',cex=2)
rgl.postscript("locationSphere2.eps","eps")
#
#
location.GCM <- Model_df.full[index,2:3]
data.GCM <- Model_df.full[index,5]
## residual:
rbf.data.GCM <- cbind(location.GCM, data.GCM)
colnames(rbf.data.GCM) <- c("lon","lat","val")
rbf.fit.GCM <- gausspr(val ~ lon + lat, data = rbf.data.GCM, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
alpha(rbf.fit.GCM) #model parameters
#
#create residuals (Obs - Predicted):
fitted.val.GCM <- predict(rbf.fit.GCM, rbf.data.GCM[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
residual_GCM <- data.GCM - fitted.val.GCM #residuals from fitted Gaussian RBF model.

#Create spatial data frame object:
residual.GCM.dataframe <- data.frame(residual_GCM)
residual.GCM.sp = SpatialPointsDataFrame(location.GCM, residual.GCM.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.GCM.sp) ~ lon + lat  #coordinate data
summary(residual.GCM.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 
#
#
residual.GCM.data <- cbind(location.GCM, residual_GCM)
#
GCM.param <- maternParam
#################################################################
#calculate the likelihood:
#
# Q: Are the three likelihoods the same?
ptm1 <- proc.time()
dist.mat <- geoDistance(location.GCM ,location.GCM ) #no NA returns here.
m.error.GCM.vec <- array(1, dim=c(n_ind,1)) #vector of ones.
Cov.mat <- maternCov(dist.mat, GCM.param, m.error.GCM.vec) 
#nudge <- 0.00001
#adj.Cov.mat <- Cov.mat + nudge*diag(n_ind)
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(residual_GCM)%*%solve(Cov.mat)%*%residual_GCM
proc.time() - ptm1

# Calculate Full Likelihood using composite code with 1x2 blocks:
ptm2 <- proc.time()
comp.NegLogLik <- composite.nll(GCM.param, residual.GCM.data,c(2,1),"Matern", m.error.GCM.vec, "NegLogLik") #Agrees with real.NegLogLik 
proc.time() - ptm2
#
#neg. log-likelihood:
full.NegLogLik <- covariance.neg.lik.calc(GCM.param, residual_GCM, m.error.GCM.vec, dist.mat, "Matern")
#
# Answer: Due to the "nudge" term because of the nearly computationally singular matrix, we find comp=full=/=real. This is to be expected.
#################################################################
# Calculate Composite Likelihood of Neighbouring blocks:
ptm3 <- proc.time()
comp.blocks <- c(5,5)
composite.NegLogLik <- composite.nll(residual.GCM.data, comp.blocks, "Matern", GCM.param, m.error.GCM.vec, "NegLogLik")  
proc.time() - ptm3

# Calculate Composite Hyvarinen score of Neighbouring blocks:
ptm3 <- proc.time()
comp.blocks <- c(15,15)
composite.NegLogLik <- composite.nll(residual.GCM.data, comp.blocks, "Matern", GCM.param, m.error.GCM.vec, "Hyvarinen")  
proc.time() - ptm3

###############################################################################
# 6. Comparing different block configurations:

# Calculate Full Likelihood using composite code with 1x2 blocks:
maternParam <- mle.output.matern.nug$par
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),] #this doesn't remove all NA?
n_ind <- 1000 #we don't want to compute this over the entire grid, at this point.
index <- round(runif(n_ind,1,dim(Model_df.full)[1]))

#Problem: Covariance matrix of some blocks comes out to be (computationally) singular.
#Soln: Added a small amount to the diagonal.
ptm1 <- proc.time()
full.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(3,2),"Matern",maternParam, m.error.GCM.vec, "NegLogLik")   #Agrees with real.NegLogLik at least. 
proc.time() - ptm1

ptm2 <- proc.time()
comp1.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(4,4),"Matern",maternParam, m.error.GCM.vec, "NegLogLik")   #Agrees with real.NegLogLik at least. 
proc.time() - ptm2

ptm3 <- proc.time()
comp2.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(10,10),"Matern",maternParam, m.error.GCM.vec, "NegLogLik")   #Agrees with real.NegLogLik at least. 
proc.time() - ptm3

################################################################
#
# COde to produce table of results:
#
Model_set='CO2_anom'
models = c('tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti')

# Specify grid dimensions for model data
# (lon_min,lon_max,lon_int,lat_min,lat_max,lat_int)
model_grid <- c(-180,178.75,1.25,-89.375,89.375,1.25) # HadCM3 ocean temps, shifted 1 hemisphere

maternParam <- mle.output.matern.nug$par
output.lik <- array(0, dim=c(8,4))
output.time <- array(0, dim=c(8,4))

n_ind <- 100
jump <- 270 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index <- seq(start, 27000, by = jump)

for(j in 1:8){ #loop over the GCMs
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
  
  m.error.GCM.vec <- 0*diag(n_ind) #identity matrix - measurement error
  
  # Calculate Full Likelihood using composite code with 1x2 blocks:
  dist.mat.gcm <- geoDistance(location.GCM, location.GCM)
  start1 <- Sys.time()
  k <- 12
  full.NegLogLik <- composite.nll(maternParam1,residual.GCM.data[1:k,],c(2,1),"Exponential", m.error.GCM.vec, "NegLogLik") 
  cov.mat2 <- expCov(dist.mat.gcm, maternParam1, m.error.GCM.vec)
  real2.NegLogLik <- 1/2*log(det(cov.mat2[1:3,1:3])) + 1/2*t(residual_GCM[1:3])%*%solve(cov.mat2[1:3,1:3])%*%residual_GCM[1:3]
  real.NegLogLik <- covariance.neg.lik.calc(maternParam1,residual_GCM[1:k], n_ind, m.error.GCM.vec, dist.mat.gcm[1:k,1:k], "Exponential") 
  end1 <- Sys.time()
  time1 <- end1 - start1
  
  start1 <- Sys.time()
  comp.NegLogLik <- composite.nll(maternParam,residual.GCM.data,c(2,2),"Matern", m.error.GCM.vec, "NegLogLik") 
  end1 <- Sys.time()
  time1 <- end1 - start1
  
  start2 <- Sys.time()
  comp1.NegLogLik <- composite.nll(maternParam, residual.GCM.data,c(3,3),"Matern", m.error.GCM.vec, "NegLogLik") 
  #Agrees with real.NegLogLik 
  end2 <- Sys.time()
  time2 <- end2 - start2
  
  start3 <- Sys.time()
  comp2.NegLogLik <- composite.nll(maternParam,residual.GCM.data,c(4,4),"Matern", m.error.GCM.vec, "NegLogLik") 
  #Agrees with real.NegLogLik 
  end3 <- Sys.time()
  time3 <- end3 - start3
  
  start4 <- Sys.time()
  comp3.NegLogLik <- composite.nll(maternParam,residual.GCM.data,c(5,5),"Matern", m.error.GCM.vec, "NegLogLik") 
  #Agrees with real.NegLogLik 
  end4 <- Sys.time()
  time4 <- end4 - start4
  
  
  #save into the table to output:
  output.lik[j,] <- c(comp.NegLogLik, comp1.NegLogLik, comp2.NegLogLik, comp3.NegLogLik)
  output.time[j,] <- c(time1, time2, time3, time4)
  
}






###################################################################
#
# Code to test method by assessing likelihood of kriging output:
#
# Choose some locations from the GCM output:
i <- 3
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[i],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=T) 
## calculating the composite likelihood:
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
#n_ind <- 500
#jump <- 54 #how regularly we want to retain locations in the grid
#start <- runif(1,1,jump)
#index <- seq(start, 27000, by = jump)
location.GCM <- Model_df.full[index,2:3] #locations to perform kriging on the sphere.
#
colnames(location.GCM) <- c("lon","lat")
location.GCM.points <- data.frame(location.GCM)
attach(location.GCM.points)
coordinates(location.GCM.points) = ~ lon + lat
model.vgm <- vgm(psill =  maternParam[3], "Exp", range = maternParam[2], nugget = maternParam[4]) 
proj4string(location.GCM.points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)

start.krige <- Sys.time()
gp.krige <- krige(residual_Obs ~ 1, residual.sp, location.GCM.points, model = model.vgm) 
end.krige <- Sys.time()
time.krige <- end.krige - start.krige

#combine location and krige:
krige.data <- cbind(location.GCM, gp.krige@data$var1.pred)
m.error.GCM.vec <- diag(n_ind) #identity matrix - measurement error

# Calculate Full Likelihood using composite code with 1x2 blocks:
start1 <- Sys.time()
comp.NegLogLik <- composite.nll(maternParam, krige.data, c(2,2), "Matern", m.error.GCM.vec, "NegLogLik") 
#Agrees with real.NegLogLik 
end1 <- Sys.time()
time1 <- end1 - start1

start2 <- Sys.time()
comp1.NegLogLik <- composite.nll(maternParam, krige.data,c(3,3),"Matern", m.error.GCM.vec, "NegLogLik") 
#Agrees with real.NegLogLik 
end2 <- Sys.time()
time2 <- end2 - start2

start3 <- Sys.time()
comp2.NegLogLik <- composite.nll(maternParam,  krige.data,c(4,4),"Matern",m.error.GCM.vec, "NegLogLik") 
#Agrees with real.NegLogLik 
end3 <- Sys.time()
time3 <- end3 - start3

start4 <- Sys.time()
comp3.NegLogLik <- composite.nll(maternParam,  krige.data,c(5,5),"Matern",m.error.GCM.vec, "NegLogLik") 
#Agrees with real.NegLogLik 
end4 <- Sys.time()
time4 <- end4 - start4

lik.row <- c(comp.NegLogLik, comp1.NegLogLik, comp2.NegLogLik, comp3.NegLogLik)
time.row <- c(time1, time2, time3, time4)

##
# CI metric alternative to composite likelihood:
jump <- 27 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index <- seq(start, 27000, by = jump)

CI.output <- array(0, dim=c(1,8))

for(j in 1:8){ #loop over the GCMs
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
  
  colnames(location.GCM) <- c("lon","lat")
  location.GCM.points <- data.frame(location.GCM)
  attach(location.GCM.points)
  coordinates(location.GCM.points) = ~ lon + lat
  model.vgm <- vgm(psill =  maternParam[3], "Exp", range = maternParam[2], nugget = maternParam[4]) 
  proj4string(location.GCM.points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
  
  gp.krige <- krige(residual_Obs ~ 1, residual.sp, location.GCM.points, model = model.vgm) 
  
  CI.krige <- CImetric(Model_df.full[index,5], gp.krige@data$var1.pred, gp.krige@data$var1.var ) #is this the variance or standard deviation?
  
  CI.output[j] <- CI.krige
  
}


#########################################################################

#########################################################################
## testing:

n.test <- dim(Model_df.edit)[1]
n.test2 <- 500
for(i in 1:n.test){
  for(j in 1:n.test){
    testing <- gdist(Model_df.edit[i,2],Model_df.edit[i,3],Model_df.edit[j,2],Model_df.edit[j,3], units="km")
    if(is.na(testing)){
      print(c(i,j))}
    }
  }

    
Model_df.full[397,]    
Model_df.full[23911,] 
test2 <- gdist(Model_df.full[index[8],2],Model_df.full[index[8],3],Model_df.full[index[21],2],Model_df.full[index[21],3])  #gives the relevant error. Changing either number fixes it...even if just by 0.001 in one direction
is.na(test2)

#which entries do not work? Why don't they work?
#1. (1,25014) and then +1 to both, i.e. (2,25015) etc...
#Actually there are more around 23000 with the 300s etc. What is happening here?




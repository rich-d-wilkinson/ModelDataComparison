## Code accompanying the tex file: FittingNonZeroTrendGPs_sphere (Mar. 2017)
#
# ---------------------------------------------------
# ---------- Pre-Amble: -----------------------------
# Load necessary functions:
source("compNLL_Func_v2.R") #this is not correct -change the filename
source("Code/Init.R")
class(geoGaussianKernel)="kernel"
#
#load packages: gstat, sp, flexclust, Imap, kernlab, geosphere and sphereplot
library(Imap, kernlab, gstat, sp, geosphere, sphereplot)
set.seed(1)
#
#
# -----------------------------------------------------
# --------- Initial Choices: --------------------------
# Choices that we can make, such as block size etc.
initial.MLE <- c(0.5,2500,1,0.01) #initial parameter values for optimisation.
matrixFunction <- c("Matern")
block.1 <- c(2,1) #Equivalent to the full negative log-likelihood
block.2 <- c(4,4)
block.3 <- c(9,9)
# For now just consider a subset of the full GCM dataset for computation time:
n_ind <- 200
jump <- 90 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index <- seq(start, 27000, by = jump)
#
# ------------------------------------------------------
# ------------- Main Body of Code: ---------------------
#
# 1. Load the data:
# i ) Observations
# Observations need to have points in a list of 4 columns with headers "x,y,z,std" = lon, lat, observation, uncertainty
# Observations are stored in directory Observation_Data/<Obs_set> with filenames <obs>.txt
Obs_set='P3+_SST_anom'
obs = c('lambda_100', 'lambda_90', 'lambda_80', 'lambda_70', 'lambda_60', 'lambda_50', 'lambda_40', 'lambda_30', 'lambda_20', 'lambda_10')
j <- 1 #choose the scaling of the stadard deviation.
# READ OBSERVATIONS
Obs <- read.table(paste0("Observation_data/",Obs_set,"/",obs[j],".txt"),header=T)
#
#create a spatial points data frame object:
locations <- cbind(Obs$x,Obs$y)
z.dataframe <- data.frame(Obs$z)
y_obs.sp = SpatialPointsDataFrame(locations, z.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(y_obs.sp) ~ lon + lat  #coordinate data
summary(y_obs.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 
#
#
# ii) Climate Models:
#
Model_set='CO2_anom'
models = c('tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti')
# Specify grid dimensions for model data : # (lon_min,lon_max,lon_int,lat_min,lat_max,lat_int)
model_grid <- c(-180,178.75,1.25,-89.375,89.375,1.25) # HadCM3 ocean temps, shifted 1 hemisphere
#
#
#Plot the observation locations on the sphere using sphereplot package:
n <- dim(Obs)[1]
rgl.sphgrid(radius = 1)
rgl.sphpoints(cbind(Obs[,1:2],array(1,dim=c(n,1))),deg=TRUE,col='black',cex=2)
rgl.postscript("locationSphere2.eps","eps")
#
#
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
#
#Create spatial data frame object:
residual.dataframe <- data.frame(residual_Obs)
residual.sp = SpatialPointsDataFrame(locations, residual.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp) ~ lon + lat  #coordinate data
summary(residual.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 
#
# # ii) Linear Model: (I believe RBFs are a better choice; this section is commented out)
# lm.data <- cbind(locations, Obs$z)
# colnames(lm.data) <- c("lon","lat","val")
# lm.data <- data.frame(lm.data)
# #
# #Coefficient estimates:
# lm.GPtrend <- lm(val ~ lon + lat, lm.data)
# residual_Obs <- lm.GPtrend$residuals
# lm.coef <- lm.GPtrend$coefficients
# 
# residual.dataframe <- data.frame(residual_Obs)
# residual.sp = SpatialPointsDataFrame(locations, residual.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
# #CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
# coordinates(residual.sp) ~ lon + lat  #coordinate data
# summary(residual.sp)
# #
#
# 3. Estimating the Parameters in the Covariance Matrix Function (Using MLE):
#First, plot the empirical variogram to inspect visually:
vgm.sample <- variogram(residual_Obs ~ 1, data = residual.sp, cutoff = 5000, width = 500 )  
plot(vgm.sample)
dev.copy(pdf,'sampleVGM_residual_rbf.pdf')
dev.off()
#
# MLE for covariance marix:
# 
# Covariance matrix function is one of Mater, Exponential and Spherical (w/ and w/o a nugget)
#
meas.error <- Obs[,4]^2 #measurement errors of observations
n <- length(residual_Obs)
dist.mat <- geoDistance(locations,locations) #create the matrix of distances between locations
#
covarianceMLE <- covariance.MLE(initial, residual_Obs, dist.mat, matrixFunction, meas.error)
#Plot the fitted model onto the empirical data points:
#Exp assumes that v=1/2 in the Matern cov. mat. function
model.choice.mat <- vgm(psill =  mle.output.matern$par[2], "Exp", range = mle.output.matern$par[3], nugget = mle.output.matern$par[4]) 
#
fit.vgm.exp <- fit.variogram(vgm.sample, model = model.choice.sph)
plot(vgm.sample, model = fit.vgm.exp)
dev.copy(pdf,'fittedVGM_mle_rbf.pdf')
dev.off()
#
#
# 4. Composite Likelihood Function for GCM:
#
output.lik <- array(0, dim=c(9,3))
output.time <- array(0, dim=c(9,3))

rbf.fit.GCM <- array(0, dim=c(n_ind, 8))

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
  rbf.fit.GCM[,j] <- gausspr(val ~ lon + lat, data = rbf.data.GCM, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
  #create residuals (Obs - Predicted):
  fitted.val.GCM <- predict(rbf.fit.GCM[,j], rbf.data.GCM[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
  residual_GCM <- data.GCM - fitted.val.GCM #residuals from fitted Gaussian RBF model.
  residual.GCM.data <- cbind(location.GCM, residual_GCM)
  m.error.GCM.vec <- array(1, dim=c(n_ind,1)) #identity matrix - measurement error
  #
  # composte likelihood calculations:
  start1 <- Sys.time()
  comp.NegLogLik <- composite.nll(covarianceMLE,residual.GCM.data,block.1,"Matern", m.error.GCM.vec, "NegLogLik") 
  end1 <- Sys.time()
  time1 <- end1 - start1
  # NOte: Block composition c(2,1) is equivalent to the usual log-likelihood
  start2 <- Sys.time()
  comp1.NegLogLik <- composite.nll(covarianceMLE, residual.GCM.data,block.2,"Matern", m.error.GCM.vec, "NegLogLik") 
  end2 <- Sys.time()
  time2 <- end2 - start2
  #
  start3 <- Sys.time()
  comp2.NegLogLik <- composite.nll(covarianceMLE, residual.GCM.data,block.3,"Matern", m.error.GCM.vec, "NegLogLik") 
  #Agrees with real.NegLogLik 
  end3 <- Sys.time()
  time3 <- end3 - start3
  #
  #save into the table to output:
  output.lik[j,] <- c(comp.NegLogLik, comp1.NegLogLik, comp2.NegLogLik)
  output.time[j,] <- c(time1, time2, time3)
}

# 5. Composite Likelihood Function using Kriged Output:
#
colnames(location.GCM) <- c("lon","lat")
location.GCM.points <- data.frame(location.GCM)
attach(location.GCM.points)
coordinates(location.GCM.points) = ~ lon + lat
#model.vgm <- vgm(psill =  maternParam[3], "Exp", range = maternParam[2], nugget = maternParam[4]) 
#I think this variogram model is the problem. We can pass a covariance function to the krige0 function in gstat:
proj4string(location.GCM.points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
start.krige <- Sys.time()
krige.output <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points, model = maternCov.krige, par.vec = covarianceMLE , m.error = meas.error )
#krige0 does not produce variance predictions at the locations. Will have to estimate this manually if we want to use the C.I. metric. 
krige.data <- cbind(location.GCM, krige.output)
end.krige <- Sys.time()
time.krige <- end.krige - start.krige
m.error.GCM.vec <- array(1, dim=c(n_ind,1)) #measurement error
# Composite likelihood calculation:
start1 <- Sys.time()
comp.krige.NegLogLik <- composite.nll(covarianceMLE, krige.data, block.1, "Matern", m.error.GCM.vec, "NegLogLik") 
end1 <- Sys.time()
time1 <- end1 - start1
#
start2 <- Sys.time()
comp1.krige.NegLogLik <- composite.nll(covarianceMLE, krige.data,block.2,"Matern", m.error.GCM.vec, "NegLogLik") 
end2 <- Sys.time()
time2 <- end2 - start2
#
start3 <- Sys.time()
comp2.krige.NegLogLik <- composite.nll(covarianceMLE,  krige.data,block.3,"Matern",m.error.GCM.vec, "NegLogLik") 
end3 <- Sys.time()
time3 <- end3 - start3
#
output.lik[9,] <- c(comp.krige.NegLogLik, comp1.krige.NegLogLik, comp2.krige.NegLogLik)
output.time[9,] <- c(time1, time2, time3)
#
#
# 6. CI metric as the score:
#
CI.output <- array(0, dim=c(1,8))
# Simulate from the fitted Gaussian Process:

krige.output <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points, model = maternCov.krige, par.vec = covarianceMLE , m.error = meas.error )
#krige0 does not produce variance predictions at the locations. Will have to estimate this manually if we want to use the C.I. metric. 
krige.data <- cbind(location.GCM, krige.output)

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
  
  fitted.val.GCM <- predict(rbf.fit.GCM[,j], rbf.data.GCM[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
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

# 7. multivariate Continuous Ranked Probability Score:



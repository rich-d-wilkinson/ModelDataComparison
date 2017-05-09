## Code accompanying the tex file: FittingNonZeroTrendGPs_sphere (Mar. 2017)

source("compNLL_Func.R")
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
rgl.postscript("locationSphere.pdf","pdf")

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
n <- length(residual_Obs)
dist.mat <- geoDistance(locations,locations)
param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000) #for validity on the sphere.
ini <- c(0.5,15,20) #how do we find initial values?
mle.output.matern <- optim(ini, matern.neg.lik.calc, data = residual_Obs, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern

#Plot the fitted model onto the empirical data points:
model.choice.mle <- vgm(psill =  2, "Exp", range = 2000) 
#This is singular if we include the nugget effect. 
fit.vgm.exp <- fit.variogram(vgm.sample, model = model.choice.mle)
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
i <- 1
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[i],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=T) 
dev.copy(pdf,'GCM1_plot.pdf')
dev.off()

#
#
## calculating the composite likelihood:
maternParam <- mle.output.matern$par
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
n_ind <- 10 #we don't want to compute this over the entire grid, at this point.
index <- round(runif(n_ind,1,dim(Model_df.full)[1]))
ptm1 <- proc.time()
dist.mat <- geoDistance(Model_df.full[index,2:3],Model_df.full[index,2:3])
Cov.mat <- maternCov(dist.mat,maternParam)  
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(Model_df.full[index,5])%*%solve(Cov.mat)%*%Model_df.full[index,5] 
proc.time() - ptm1

# Calculate Full Likelihood using composite code with 1x2 blocks:
ptm2 <- proc.time()
full.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(2,1),"Matern",maternParam) #Agrees with real.NegLogLik at least. 
proc.time() - ptm2

# Calculate Composite Likelihood of Nieghbouring blocks:
ptm3 <- proc.time()
comp.blocks <- c(3,2)
composite.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)], comp.blocks, "Matern", maternParam)  
proc.time() - ptm3

# 6. Comparing different block configurations:

# Calculate Full Likelihood using composite code with 1x2 blocks:
maternParam <- mle.output.matern$par
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
n_ind <- 1000 #we don't want to compute this over the entire grid, at this point.
index <- round(runif(n_ind,1,dim(Model_df.full)[1]))

ptm1 <- proc.time()
full.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(2,1),"Matern",maternParam) #Agrees with real.NegLogLik at least. 
proc.time() - ptm1

ptm2 <- proc.time()
comp1.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(4,4),"Matern",maternParam) #Agrees with real.NegLogLik at least. 
proc.time() - ptm2

ptm3 <- proc.time()
comp2.NegLogLik <- composite.nll(Model_df.full[index,c(2:3,5)],c(10,10),"Matern",maternParam) #Agrees with real.NegLogLik at least. 
proc.time() - ptm3

#########################################################################
## testing:

for(i in 1:25){
  for(j in 1:25){
    testing <- gdist(Model_df.full[index[i],2],Model_df.full[index[i],3],Model_df.full[index[j],2],Model_df.full[index[j],3], units="km")
    print(c(i,j))
    print(testing)
  }
}
    
Model_df.full[index[8],]    
Model_df.full[index[21],] 
test2 <- gdist(-177.5,-70.625,2.5,70.625)  #gives the relevant error. Changing either number fixes it...even if just by 0.001 in one direction
is.na(test2)


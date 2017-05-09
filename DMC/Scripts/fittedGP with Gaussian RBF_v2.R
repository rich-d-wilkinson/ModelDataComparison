## Estimating the mean function trend for the geospatial observation data:
# There might be a problem with the data being on the sphere, I'm not sure. For now, we just ignore it.
# We use the Gaussian Radial Basis Function commmand rbfdot within the kernlab package.
install.packages('kernlab')

#kernlab, flexclust, gstat, sp

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

#create a spatial points data frame object:
z.dataframe <- data.frame(Obs$z)
y_obs.sp = SpatialPointsDataFrame(cbind(Obs$x,Obs$y), z.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(y_obs.sp) ~ lon + lat  #coordinate data
summary(y_obs.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 


# 2. Fit the Gaussian RBF using kernlab:
#
# We fit the mean function with a Gausssian radial basis function: f(x) = sum_{i=1}^n alpha_i exp( - sigma L2_norm(x - x_i) )
covariate <- cbind(Obs$x, Obs$y)
response <- Obs$z
rbf.fit <- gausspr(covariate, response, kernel = "rbfdot", var = 1) #uses automatic sigma estimation? What does this mean?

##############
locations <- cbind(Obs$x,Obs$y)
rbf.data <- cbind(locations, Obs$z)
colnames(rbf.data) <- c("lon","lat","val")
rbf.fit2 <- gausspr(val ~ lon + lat, data = rbf.data, kernel = "rbfdot", var = 1) #uses automatic sigma estimation? What does this mean?
alpha(rbf.fit2)[1:5,]
#This gives a different answer. Not sure what is happening here. 
#############
alpha(rbf.fit)[1:5,] #model parameters


#3. Fitted variogram for the residuals of the regression: GP with mean zero and unknown covariance matrix function
#
fitted.val <- predict(rbf.fit, covariate, type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
residual_Obs <- Obs$z - fitted.val #residuals from fitted Gaussian RBF model.

#Create spatial data frame object:
residual.dataframe <- data.frame(residual_Obs)
residual.sp = SpatialPointsDataFrame(cbind(Obs$x,Obs$y), residual.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp) ~ lon + lat  #coordinate data
summary(residual.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 

####################################################
#sample variogram:
vgm.sample <- variogram(residual_Obs ~ 1, data = residual.sp, cutoff = 25000, width = 500 )  
plot(vgm.sample)
#Does this even make sense for the sample variogram? What are we looking at here?

######################################################################
#4. MLE to find the covariance matrix function parameters:
# maximum likelihood estimation using a Gaussian likelihood with the covariance matrix function of our choice, the data is the observation residuals. 
locations <- cbind(Obs$x,Obs$y)
n <- length(residual_Obs)
dist.mat <- geoDistance(locations,locations)
param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000) #for validity on the sphere.
ini <- c(0.5,15,20) #how do we find initial values?
mle.output.matern <- optim(ini, matern.neg.lik.calc, data = residual_Obs, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern

#we might also try exponential model:
param.lower <- c(0.01,0.01)
ini <- c(15,20) #how do we find initial values?
mle.output.matern <- optim(ini, sub.matern.neg.lik.calc, nu = 0.5, data = residual_Obs, n = n, dist.mat = dist.mat, lower = param.lower, method="L-BFGS-B")
mle.output.matern



################################################################################
################################################################################
#user-defined kernel with geodesic distance in Gaussian RBF:
#r is the distance between 
geoGaussianKernel=function(x,y){
    #calculate the geodesic distance between x and y
    x <- as.matrix(t(x))
    y <- as.matrix(t(y)) #these have to be passed in as matrices for the dimension call to work.
    r <- geoDistance(x,y)
    return(exp(-(r^2))) #ignore kernel parameter for now
}

class(geoGaussianKernel)="kernel"

# We fit the mean function with a Gausssian radial basis function: f(x) = sum_{i=1}^n alpha_i exp( - sigma L2_norm(x - x_i) )
locations <- cbind(Obs$x,Obs$y)
rbf.data <- cbind(locations, Obs$z)
colnames(rbf.data) <- c("lon","lat","val")
rbf.data.subset <- rbf.data[1:5,]
#rbf.fit <- gausspr(x = locations[1:5,], y = Obs$z[1:5], kernel = geoGaussianKernel, kpar = list(sigma = 1)) #uses automatic sigma estimation? What does this mean?
rbf.fit.geo <- gausspr(val ~ lon + lat, data = rbf.data.subset, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
#rbf.fit3 <- gausspr(covariate[1:5,], rbf.data, response[1:5], kernel = geoGaussianKernel, type="regression",kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
alpha(rbf.fit.geo) #model parameters

#Alternatively, give the kernel matrix to gausspr:
#don't know how this works.
#ker.mat <- geoGaussianKernel(covariate[1:5,],covariate[1:5,])
#rbf.fit <- gausspr(ker.mat, response[1:5], kernel = 'matrix', var = 0.5) #uses automatic sigma estimation? What does this mean?

fitted.val <- predict(rbf.fit.geo, rbf.data.subset[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
residual_Obs <- Obs$z[1:5] - fitted.val #residuals from fitted Gaussian RBF model.

#Create spatial data frame object:
residual.dataframe <- data.frame(residual_Obs)
residual.sp = SpatialPointsDataFrame(cbind(Obs$x,Obs$y)[1:5,], residual.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp) ~ lon + lat  #coordinate data
summary(residual.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 

#MLE for covariance marix:
locations <- cbind(Obs$x,Obs$y)[1:5,]
n <- length(residual_Obs)
dist.mat <- geoDistance(locations,locations)
param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000) #for validity on the sphere.
ini <- c(0.5,15,20) #how do we find initial values?
mle.output.matern <- optim(ini, matern.neg.lik.calc, data = residual_Obs, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern

#we might also try exponential model:

#####################
param.lower <- c(0.01,0.01)
ini <- c(15,20) #how do we find initial values?
mle.output.matern <- optim(ini, sub.matern.neg.lik.calc, nu = 0.5, data = residual_Obs, n = n, dist.mat = dist.mat, lower = param.lower, method="L-BFGS-B")
mle.output.matern
##########################################################
###############################################################################
#
# Fitting GP with a trend:
#
locations <- cbind(Obs$x,Obs$y)
lm.data <- cbind(locations, Obs$z)
colnames(lm.data) <- c("lon","lat","val")
lm.data <- data.frame(lm.data)
#
#1. find coefficient estimates
lm.GPtrend <- lm(val ~ lon + lat, lm.data)
lm.residuals <- lm.GPtrend$residuals
lm.coef <- lm.GPtrend$coefficients
#
#2. fit residuals to GP model
n <- length(lm.residuals)
dist.mat <- geoDistance(locations,locations)
param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000) #for validity on the sphere.
ini <- c(0.5,15,20) #how do we find initial values?
mle.output.matern <- optim(ini, matern.neg.lik.calc, data = lm.residuals, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern

#we might also try exponential model:
param.lower <- c(0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(10000,100000) 
ini <- c(15,20) #how do we find initial values?
mle.output.matern <- optim(ini.true, sub.matern.neg.lik.calc, data = lm.residuals, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern



ptm <- proc.time() #time the operation to compare with the full likelihood function.
test <- geoDistance(lm.data[,1:2],lm.data[,1:2])
proc.time() - ptm

ptm2 <- proc.time() #time the operation to compare with the full likelihood function.
test <- geoDistance(lm.data[,1:2],lm.data[,1:2])
proc.time() - ptm2


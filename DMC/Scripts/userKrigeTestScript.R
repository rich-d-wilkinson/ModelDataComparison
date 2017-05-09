##### 
# Replicating the krige0 example using meuse in sp package:
#
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

v = function(x, y = x) { 
  #print(coordinates(x[1:3,]))
  print(coordinates(y[1:3,]))
  b <- exp(-spDists(coordinates(x),coordinates(y))/500)
  print(dim(b))
  return(b)
  }
x <- krige0(zinc~1, meuse, meuse.grid, v)


#######################################################################
#
# Testing:
colnames(location.GCM) <- c("lon","lat")
location.GCM.points <- data.frame(location.GCM)
attach(location.GCM.points)
coordinates(location.GCM.points) = ~ lon + lat
proj4string(location.GCM.points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)


test.dist.krige <- function(data1,data2){
  #print(coordinates(data1[1:3,]))
  #print(coordinates(data2[1:3,]))
  d <- geoDistance(coordinates(data1), coordinates(data2))
  out <- exp(-1/500*d)
  return(out)
}

test.dist.krige.sig <- function(data1,data2,sigma.test, m.error){
  #print(coordinates(data1[1:3,]))
  #print(coordinates(data2[1:3,]))
  nugget <- 0
  d <- geoDistance(coordinates(data1), coordinates(data2))
  #if matrix is symmetric then we are in loop1 and add measurement error:
  out <- exp(-1/sigma.test*d)
  if( abs(d[1,2] - d[2,1]) < 0.0001){
    out <- out + nugget*m.error
    print("A")
  }
  else{ print("B") }
  print(dim(out))
  print(dim(d))
  return(out)
}

maternCov.krige <- function(data.obs, data.GCM, par.vec, m.error){
  n.1 <- dim(coordinates(data.obs))[1]
  n.2 <- dim(coordinates(data.GCM))[1]
  
  mat.error <- diag(m.error)
  if(length(par.vec) == 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    d <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    if( abs(d[1,2] - d[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2]) + nugget*mat.error
      index_0 <- which(d > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
      print("A")
    }
    else{
      out <- sigma^2*exp(-d.full/rho)#element-wise operations.
      print("B")
    }
    
  }
  else if (length(par.vec) > 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    nugget <- par.vec[3]
    d <- geoDistance(coordinates(data.obs), coordinates(data.GCM))

    if( abs(d[1,2] - d[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2]) + nugget*mat.error
      index_0 <- which(d > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
      print("A")
    }
    else{
      out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
      print("B")
    }
    
  }
  
  return(out)
}

expCov.krige <- function(data.obs, data.GCM, par.vec, m.error){
  n.1 <- dim(coordinates(data.obs))[1]
  n.2 <- dim(coordinates(data.GCM))[1]

  mat.error <- diag(m.error)
  if(length(par.vec) == 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    d.full <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    if( abs(d.full[1,2] - d.full[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d.full)[1], ncol = dim(d.full)[2]) + nugget*mat.error
      index_0 <- which(d.full > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*exp(-d.full[index_0]/rho)#element-wise operations.
      print("A")
    }
    else{
      out <- sigma^2*exp(-d.full/rho)#element-wise operations.
      print("B")
    }
    
  }
  else if (length(par.vec) > 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    nugget <- par.vec[3]
    d.full <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    if( abs(d.full[1,2] - d.full[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d.full)[1], ncol = dim(d.full)[2]) + nugget*mat.error
      index_0 <- which(d.full > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*exp(-d.full[index_0]/rho)#element-wise operations.
      print("A")
    }
    else{
      out <- sigma^2*exp(-d.full/rho)#element-wise operations.
      print("B")
    }
    
  }
  
  return(out)
}

test.krige <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points[1:5,], model = test.dist.krige)
test.krige.sig <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points[1:2,], model = test.dist.krige.sig, sigma.test = 500, m.error = meas.error)
#####################################################################################
## GCM version:
meas.error <- Obs[,4]^2 
krige.output <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points, model = expCov.krige, par.vec = maternParam[2:4] , m.error = meas.error )
krige.data <- cbind(location.GCM, krige.output)




#########################################################################
############################################################################
# testing parameters:

d.test <- geoDistance(location.GCM[1:10,], location.GCM[1:10,])
rho <- 1000
sigma <- 1.54
test.mat <- sigma^2*exp(-d.test/rho)


#########################################################################################################
#########################################################################################################
# What's wrong with the kriging function?
#
krige0(residual_Obs ~ 1, residual.sp, location.GCM.points[1:5,], model = maternCov.krige, par.vec = covarianceMLE$par , m.error = meas.error )
krige0(residual_Obs ~ 1, residual.sp, location.GCM.points[1:5,], model = maternCov.krige, par.vec = covarianceMLE$par , m.error = meas.error )

krige0(residual_Obs ~ 1, residual.sp, location.GCM.points[1:5,], computeVar = TRUE, model = maternCov.krige, par.vec = covarianceMLE$par , m.error = meas.error)
#dont know how computeVar works. 






###################################################################################################
####################################################################################################
#
# can we recover the covariance matrix parameters?
testing.par <- c(2500, 1)
jump <- 10 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index.test <- seq(start, 27000, by = jump)
n_ind = length(index.test)
j <- 1
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[j],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=F) 

## calculating the composite likelihood:
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
location.GCM.test <- Model_df.full[index.test,2:3]
#
colnames(location.GCM.test) <- c("lon","lat")
location.GCM.points.test <- data.frame(location.GCM.test)
attach(location.GCM.points.test)
coordinates(location.GCM.points.test) = ~ lon + lat
proj4string(location.GCM.points.test) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
#
testing.krige <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points.test, model = expCov.krige, par.vec = testing.par , m.error = meas.error )
dist.mat.test <- geoDistance(location.GCM.test, location.GCM.test)
#initial.MLE <- array(c(0.5,0.5,0.5,1000,2500,3500,1,1,1), dim=c(3,3)) #initial parameter values for optimisation.
initial.MLE <- array(c(1000,2500,3500,1,1,1), dim=c(3,2))
matrixFunction <- c("Exponential")
initial.MLE.test <- rbind(initial.MLE, testing.par)
#testing.krige.comp <- testing.krige[-which(is.na(testing.krige)),]
#dist.mat.test.comp <- geoDistance(location.GCM.test[-which(is.na(testing.krige)),], location.GCM.test[-which(is.na(testing.krige)),])

covarianceMLE.test <- covariance.MLE(initial.MLE.test, testing.krige, dist.mat.test, matrixFunction, meas.error)

#initial.MLE.NM <- array(c(log(0.4),log(0.4),log(0.4),log(2000),log(2500),log(3500),1,1,1), dim=c(3,3)) #initial parameter values for optimisation.
#testing.par.NM <- array(c(log(0.3), log(2500), 1), dim=c(1,3))
initial.MLE.NM <- array(c(log(2000),log(2500),log(3500),1,1.5,0.8), dim=c(3,2)) #initial parameter values for optimisation.
testing.par.NM <- array(c(log(2500), 1), dim=c(1,2))
matrixFunction <- c("Exponential")
initial.MLE.NM.test <- rbind(initial.MLE.NM, testing.par.NM)
covarianceMLE.NM.test <- covariance.MLE.NelderMead(initial.MLE.NM.test, testing.krige, dist.mat.test, matrixFunction, meas.error)
covarianceMLE.NM.log.test <- c(exp(covarianceMLE.NM.test$par[1]), covarianceMLE.NM.test$par[2]) 

###############################################
# using more data points:
jump <- 135 #how regularly we want to retain locations in the grid
start <- floor(runif(1,1, jump))
index.test <- seq(start, 27000, by = jump)
n_ind <- length(index)
##
#pseudo-observations:
#
Model_set='CO2_anom'
models = c('tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti')
# Specify grid dimensions for model data : # (lon_min,lon_max,lon_int,lat_min,lat_max,lat_int)
model_grid <- c(-180,178.75,1.25,-89.375,89.375,1.25) # HadCM3 ocean temps, shifted 1 hemisphere
j <- 1 #hopefully the NAs are the smae over the 8 GCM datasets.
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[j],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=F) 

#remove NA for now:
Model_df.full.pseudo <- Model_df[complete.cases(Model_df),]
n.pseudo <- 250
num.krigeLoc <- dim(Model_df.full.pseudo[-index,])[1] #don't want to choose any location twice.
vec.krigeLoc <- c(1:num.krigeLoc)
index.pseudo.test <- sample(vec.krigeLoc, n.pseudo, replace = FALSE)

#data:
location.data.test <- Model_df.full[index.pseudo.test,2:3]
location.krige.test <- Model_df.full[index.test,2:3]
obs.data.test <- Model_df.full[index.pseudo.test,5]
#

###
rbf.data.test <- cbind(location.data.test, obs.data.test)
colnames(rbf.data.test) <- c("lon","lat","val")
rbf.fit.test <- gausspr(val ~ lon + lat, data = rbf.data.test, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
#
#create residuals (Obs - Predicted):
fitted.val.test <- predict(rbf.fit.test, rbf.data.test[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
residual_Obs.test <- obs.data.test - fitted.val.test #residuals from fitted Gaussian RBF model.
#
residual.test.dataframe <- data.frame(residual_Obs.test)
residual.sp.test = SpatialPointsDataFrame(location.data.test, residual.test.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp.test) ~ lon + lat  #coordinate data
##
colnames(location.krige.test) <- c("lon","lat")
location.krige.points.test <- data.frame(location.krige.test)
attach(location.krige.points.test)
coordinates(location.krige.points.test) = ~ lon + lat
proj4string(location.krige.points.test) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
#
meas.error <- array(1, dim=c(n.pseudo,1))
testing.krige <- krige0(residual_Obs.test ~ 1, residual.sp.test, location.krige.points.test, model = expCov.krige, par.vec = testing.par , m.error = meas.error )

dist.mat.test <- geoDistance(location.krige.test, location.krige.test)
#initial.MLE <- array(c(0.5,0.5,0.5,1000,2500,3500,1,1,1), dim=c(3,3)) #initial parameter values for optimisation.
#initial.MLE <- array(c(1000,2500,3500,1,1,1), dim=c(3,2))
#matrixFunction <- c("Exponential")
#initial.MLE.test <- rbind(initial.MLE, testing.par)
#testing.krige.comp <- testing.krige[-which(is.na(testing.krige)),]
#dist.mat.test.comp <- geoDistance(location.GCM.test[-which(is.na(testing.krige)),], location.GCM.test[-which(is.na(testing.krige)),])

#covarianceMLE.test <- covariance.MLE(initial.MLE.test, testing.krige, dist.mat.test, matrixFunction, meas.error)

#Nelder-Mead:
#initial.MLE.NM <- array(c(log(0.4),log(0.4),log(0.4),log(2000),log(2500),log(3500),1,1,1), dim=c(3,3)) #initial parameter values for optimisation.
#testing.par.NM <- array(c(log(0.3), log(2500), 1), dim=c(1,3))
initial.MLE.NM <- array(c(log(2000),log(2500),log(3500),1,1.5,0.8), dim=c(3,2)) #initial parameter values for optimisation.
testing.par.NM <- array(c(log(2500), 1), dim=c(1,2))
matrixFunction <- c("Exponential")
initial.MLE.NM.test <- rbind(initial.MLE.NM, testing.par.NM)
covarianceMLE.NM.test <- covariance.MLE.NelderMead(initial.MLE.NM.test, testing.krige, dist.mat.test, matrixFunction, meas.error)
covarianceMLE.NM.log.test <- c(exp(covarianceMLE.NM.test$par[1]), covarianceMLE.NM.test$par[2]) 


###################################################################################################
####################################################################################################
#
# Checking the kriging function:
#Create spatial data frame object:
residual.krigeTest <- residual_Obs[49:50]
residual.krigeTest <- c(1,2)
residual.dataframe.krigeTest <- data.frame(residual.krigeTest)
locations.krigeTest <- array(c(-12,-12,80,50), dim=c(2,2))
residual.sp.krigeTest = SpatialPointsDataFrame(locations.krigeTest, residual.dataframe.krigeTest, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(residual.sp.krigeTest) ~ lon + lat  #coordinate data

location.GCM.krigeTest <- array(c(-12,-12,70,60), dim=c(2,2))
colnames(location.GCM.krigeTest) <- c("lon","lat")
location.GCM.sp.krigeTest <- data.frame(location.GCM.krigeTest)
attach(location.GCM.sp.krigeTest)
coordinates(location.GCM.sp.krigeTest) = ~ lon + lat
proj4string(location.GCM.sp.krigeTest) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)

param.krigeTest <- c(1000, 0.01)
set.seed(1)
testing.krige.1 <- krige0(residual.krigeTest ~ 1, residual.sp.krigeTest, location.GCM.sp.krigeTest, model = expCov.krige, par.vec = param.krigeTest, m.error = c(1,1)) #, computeVar = TRUE, fullCovariance = TRUE)

set.seed(100)
testing.krige.100 <- krige0(residual.krigeTest ~ 1, residual.sp.krigeTest, location.GCM.sp.krigeTest, model = expCov.krige, par.vec = param.krigeTest, m.error = c(1,1)) #, computeVar = TRUE, fullCovariance = TRUE)
#It looks like kriging only produces the mean of the GP, and will not produce the variance at the locations. 


############################################################
###########################################################
#
# How does outer work?
x <- array(c(10,-31,42,83,-60,30), dim=c(3,2))
y <- array(c(-33,11,2,-83,-14,55), dim=c(3,2))
func.op <- function(x,y,case){
  if(case == 1){
    x[,1]*y[,1] + x[,2]*y[,2]
  }
  else if(case == 2){
    x - y
  }
}
case.id <- 1
outer(x,y,"func.op", case = case.id)

#testing distance matrix calculation:
geo1 <- array(c(10,-31,42,83,-60,30), dim=c(3,2))
geo2 <- array(c(-33,11,2,-83,-14,55), dim=c(3,2))

outer(geo1, geo2, "geoDistance")


############################################################
############################################################
#
# Global Optimisation Search Procedures:
#
# Prelim:
testing.par <- c(2500, 1)
jump <- 270 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index.test <- seq(start, 27000, by = jump)
n_ind = length(index.test)
j <- 1
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[j],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=F) 

## calculating the composite likelihood:
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
location.GCM.test <- Model_df.full[index.test,2:3]
#
colnames(location.GCM.test) <- c("lon","lat")
location.GCM.points.test <- data.frame(location.GCM.test)
attach(location.GCM.points.test)
coordinates(location.GCM.points.test) = ~ lon + lat
proj4string(location.GCM.points.test) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
#
testing.krige <- krige0(residual_Obs ~ 1, residual.sp, location.GCM.points.test, model = expCov.krige, par.vec = testing.par , m.error = meas.error )
dist.mat.test <- geoDistance(location.GCM.test, location.GCM.test)
#initial.MLE <- array(c(0.5,0.5,0.5,1000,2500,3500,1,1,1), dim=c(3,3)) #initial parameter values for optimisation.
initial.MLE <- array(c(1000,2500,3500,1,1,1), dim=c(3,2))
matrixFunction <- c("Exponential")
initial.MLE.test <- rbind(initial.MLE, testing.par)

#
# 1. Simulated Annealing
lower.b <- c(500,0.01)
upper.b <- c(15000, 5)
covarianceMLE.test.sa <- covariance.MLE.SA(initial.MLE.test, testing.krige, lower.b, upper.b, dist.mat.test, matrixFunction, meas.error)



#
# 2. Particle Swarm Optimisation
covarianceMLE.test.pso <- covariance.MLE.PSO(array(1, dim=c(1,2)), testing.krige, lower.b, upper.b, dist.mat.test, matrixFunction, meas.error)





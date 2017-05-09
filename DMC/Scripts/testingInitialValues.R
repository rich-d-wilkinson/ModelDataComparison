#Testing initial values:

n <- length(residual_Obs)
dist.mat <- geoDistance(locations,locations)
param.lower <- c(0.01,0.01,0.01,0.00001) #(nu,rho,sigma,nugget)
#needs upper limit for nu <= 1/2 too.
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,1000,100) #for validity on the sphere.
#no nugget:
ini.3 <- c(0.5,20,1)
param.upper.3 <- c(1/2,10000,10000) #for validity on the sphere.
param.lower.3 <- c(0.01,0.01,0.01) #(nu,rho,sigma)
ini <- c(0.5, 15, 1, 0.1)
ini <- c(0.4,2500,1,0.1) #how do we find initial values?
ini.noise <- 0.01*Obs$std
ini.noise.const <- array(0.1,dim=c(1,n))
ini.noisy <- c(0.45,2500,1,ini.noise)
ini.noisy <- c(0.45,15,1,ini.noise.const)
mle.output.matern <- optim(ini.3, covariance.neg.lik.calc, data = residual_Obs, n = n, dist.mat = dist.mat, model = "Matern", lower = param.lower.3, 
                           upper = param.upper.3, method="L-BFGS-B")
mle.output.matern
mle.output.spherical <- optim(ini, covariance.neg.lik.calc, data = residual_Obs, n = n, dist.mat = dist.mat, model = "Spherical", lower = param.lower,
                              upper = param.upper, method="L-BFGS-B")
mle.output.spherical


#######################################
Model_set='CO2_anom'
models = c('tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti')

model_grid <- c(-180,178.75,1.25,-89.375,89.375,1.25) # HadCM3 ocean temps, shifted 1 hemisphere

i <- 3
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[i],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid, plot_it=T) 

#
## calculating the composite likelihood:
maternParam <- mle.output.matern$par

#remove NA for now:
dim(Model_df)
Model_df.full <- Model_df[complete.cases(Model_df),]
dim(Model_df.full)

#Can we run the distanve calculations?
n.test <- dim(Model_df.full)[1]
n.test2 <- 500
for(i in 1:ntest){
  for(j in 1:n.test){
    testing <- gdist(Model_df.full[i,2],Model_df.full[i,3],Model_df.full[j,2],Model_df.full[j,3], units="km")
    if(is.na(testing)){
      print(c(i,j))}
  }
}
#No, there are some terms that we can't calculate. Why? e.g. (1, 25014) and then add 1 to both sides.
Model_df.full[1,]
Model_df.full[25014,]

Model_df.full[25,]
Model_df.full[24783,]
gdist(Model_df.full[25,2],Model_df.full[25,3],Model_df.full[24783,2],Model_df.full[24783,3], units="km") #returns NA


######################################################
n_ind <- 250 #we don't want to compute this over the entire grid, at this point.
index <- round(runif(n_ind,1,dim(Model_df.full)[1]))
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
GCM.param <- mle.output.matern$par
#calculate the likelihood:
ptm1 <- proc.time()
dist.mat <- geoDistance(location.GCM ,location.GCM ) #no NA returns here.
m.error.GCM <- diag(n_ind) #identity matrix
Cov.mat <- maternCov(dist.mat,GCM.param, m.error.GCM)  
real.NegLogLik <- 1/2*log(det(Cov.mat)) + 1/2*t(residual_GCM)%*%solve(Cov.mat)%*%residual_GCM
proc.time() - ptm1

# Calculate Full Likelihood using composite code with 1x2 blocks:
ptm2 <- proc.time()
m.error.GCM.vec <- array(1, dim=c(1, n_ind)) #identity matrix
full.NegLogLik <- composite.nll(residual.GCM.data,c(2,1),"Matern",GCM.param, m.error.GCM.vec) #Agrees with real.NegLogLik 
proc.time() - ptm2



###########################################################
##########################################################
#
# Looking at MLE and the different GCM outputs:
# How flat is the likelihood? Do we find the parameter that maximises likelihood?
# What is the likelihood 
covarianceMLE.vec <- array(1, dim=c(8,2))


for(j in 1:8){ #loop over the GCMs 
  #1. Create the pseudo-obsevations for each GCM:
  Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[j],".txt"),
                              indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                              model_grid, plot_it=F) 
  
  #remove NA for now:
  Model_df.full.pseudo <- Model_df[complete.cases(Model_df),]
  location.pseudo <- Model_df.full.pseudo[index.pseudo,2:3]
  data.pseudo <- Model_df.full.pseudo[index.pseudo,5]
  
  pseudo.GCM <- data.pseudo #+ noise
  #
  #
  rbf.data.pseudo <- cbind(location.pseudo, pseudo.GCM)
  colnames(rbf.data.pseudo) <- c("lon","lat","val")
  rbf.fit.pseudo <- gausspr(val ~ lon + lat, data = rbf.data.pseudo, kernel = geoGaussianKernel, kpar=list(sigma=1)) #uses automatic sigma estimation? What does this mean?
  #create residuals (Obs - Predicted):
  fitted.pseudo <- predict(rbf.fit.pseudo, rbf.data.pseudo[,1:2], type = "response") #kernlab predict not gstat (pass gausspr object in first variable)
  residual_pseudo <- pseudo.GCM - fitted.pseudo #residuals from fitted Gaussian RBF model.
  #Create spatial data frame object:
  residual.dataframe.pseudo <- data.frame(residual_pseudo)
  residual.sp.pseudo = SpatialPointsDataFrame(location.pseudo, residual.dataframe.pseudo, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
  coordinates(residual.sp.pseudo) ~ lon + lat  #coordinate data
  #
  #meas.error.pseudo <- Obs[,4]^2 #measurement errors from the real observations (could change all of these to ones)
  meas.error.pseudo <- array(1, dim=c(n.pseudo,1))
  dist.mat.pseudo <- geoDistance(location.pseudo,location.pseudo) #create the matrix of distances between locations
  #
  # Newton-Raphson:
  cov.MLE <- covariance.MLE(initial.MLE[,2:3], residual_pseudo, dist.mat.pseudo, matrixFunction, meas.error.pseudo)
  covarianceMLE.vec[j,] <- cov.MLE$par
}

k <- 4
  Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[k],".txt"),
                              indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                              model_grid, plot_it=F) 
  
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
  

  dist.mat <- geoDistance(location.GCM,location.GCM) #create the matrix of distances between locations
  rho.testing <- seq(1000,5000,by=500)
  test.neglik.pseudo <- array(1, dim=c(8,1))
  for(k in 1:8){
     Cov.mat <- expCov(dist.mat, covarianceMLE.vec[k,], m.error.GCM.vec) 
    if (det(Cov.mat) < 0.000000000000001){
      print("Singular")
      dist.mat <- geoDistance(location.GCM + 0.001,location.GCM + 0.001) #create the matrix of distances between locations
      Cov.mat <- expCov(dist.mat, covarianceMLE.vec[k,], m.error.GCM.vec) 
    }
    test.neglik.pseudo[k] <- 1/2*log(det(Cov.mat)) + 1/2*t(residual_GCM)%*%solve(Cov.mat)%*%residual_GCM
  
} 
  
  
  
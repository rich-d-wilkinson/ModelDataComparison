############################################################################
# Pseudo-observations for checking the composite likelihood approach : ####
############################################################################
#
#This version has the RBF approximation to the mean removed to see if this affects any results.
#
# ---------------------------------------------------
# ---------- Pre-Amble: -----------------------------
# Load necessary functions:
source("compNLL_Func_v2_testing.R") #this is not correct -change the filename
source("Code/Init.R")
class(geoGaussianKernel)="kernel"
#
#load packages: gstat, sp, flexclust, Imap, kernlab, geosphere and sphereplot
library(Imap, kernlab, gstat, geosphere, sphereplot, plotrix)
set.seed(1)
#
#
# -----------------------------------------------------
# --------- Initial Choices: --------------------------
# Choices that we can make, such as block size etc.
K <- 25 #number of initial starting points for optimisation to try.
initial.MLE.nu <- array(0.5,dim = c(K,1)) #initial parameter values for optimisation.
initial.MLE.sigma <- array(1, dim = c(K,1))
initial.MLE.rho <- runif(K, 1000, 10000)
initial.MLE <- cbind(initial.MLE.nu, initial.MLE.rho, initial.MLE.sigma)
initial.MLE.NM <- cbind(log(initial.MLE.nu), log(initial.MLE.rho), initial.MLE.sigma )
matrixFunction <- c("Matern")
block.1 <- c(2,1) #Equivalent to the full negative log-likelihood
block.2 <- c(4,4)
block.3 <- c(9,9)
# For now just consider a subset of the full GCM dataset for computation time:

jump <- 135 #how regularly we want to retain locations in the grid
start <- floor(runif(1,1, jump))
index <- seq(start, 27000, by = jump)
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
n.pseudo <- 95
num.krigeLoc <- dim(Model_df.full.pseudo[-index,])[1] #don't want to choose any location twice.
vec.krigeLoc <- c(1:num.krigeLoc)
index.pseudo <- sample(vec.krigeLoc, n.pseudo, replace = FALSE)

#
#
Obs_set='P3+_SST_anom'
obs = c('lambda_100', 'lambda_90', 'lambda_80', 'lambda_70', 'lambda_60', 'lambda_50', 'lambda_40', 'lambda_30', 'lambda_20', 'lambda_10')
i <- 1 #choose the scaling of the stadard deviation.
# READ OBSERVATIONS
Obs <- read.table(paste0("Observation_data/",Obs_set,"/",obs[i],".txt"),header=T)
#
#output.lik.pseudo <- array(0, dim=c(8,3,8))
#output.time.pseudo <- array(0, dim=c(8,3,8))
full.neglik.pseudo <- array(0, dim=c(9,8)) #matrix should be largest along the diagonal. The final row is the kriged observations from the fitted GP
#
# ------------------------------------------------------
# ------------- Main Body of Code: ---------------------
#
# Generate the pseudo-observations for each GCM at the locations:
#
noise <- runif(n.pseudo, min = 0, max = 0.02)   #This is a small amount of noise we add to the observations to create the pseudo-observations
for(j in 1:2){ #loop over the GCMs 
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
  cbind.data.pseudo <- cbind(location.pseudo, pseudo.GCM)
  colnames(cbind.data.pseudo) <- c("lon","lat","val")
  
  dataframe.pseudo <- data.frame(pseudo.GCM)
  sp.pseudo = SpatialPointsDataFrame(location.pseudo, dataframe.pseudo, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
  coordinates(sp.pseudo) ~ lon + lat  #coordinate data
  #
  #meas.error.pseudo <- Obs[,4]^2 #measurement errors from the real observations (could change all of these to ones)
  meas.error.pseudo <- array(1, dim=c(n.pseudo,1))
  dist.mat.pseudo <- geoDistance(location.pseudo,location.pseudo) #create the matrix of distances between locations
  #
  # Newton-Raphson:
  covarianceMLE <- covariance.MLE(initial.MLE, pseudo.GCM, dist.mat.pseudo, matrixFunction, meas.error.pseudo)
  #print(covarianceMLE$par)
  #
  # Nelder-Mead:
  #covarianceMLE.NM <- covariance.MLE.NelderMead(initial.MLE.NM, residual_pseudo, dist.mat.pseudo, matrixFunction, meas.error.pseudo)
  #print(c(exp(covarianceMLE.NM$par[1]), exp(covarianceMLE.NM$par[2]), covarianceMLE.NM$par[3])) 
  #2. Calculate the composite likelihood for each GCM:

  for(k in 1:4){ #loop over the GCMs
    Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[k],".txt"),
                                indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                                model_grid, plot_it=F) 
    
    ## calculating the composite likelihood:
    #remove NA for now:
    Model_df.full <- Model_df[complete.cases(Model_df),]
    
    location.GCM <- Model_df.full[index,2:3]
    data.GCM <- Model_df.full[index,5]
    ## residual:
    nonZeroMean.data.GCM <- cbind(location.GCM, data.GCM)
    colnames(nonZeroMean.data.GCM) <- c("lon","lat","val")

    m.error.GCM.vec <- array(1, dim=c(n_ind,1)) #identity matrix - measurement error
    #
    # # composte likelihood calculations:
    # start1 <- Sys.time()
    # comp.NegLogLik <- composite.nll(covarianceMLE$par,residual.GCM.data,block.1,"Matern", m.error.GCM.vec, "NegLogLik") 
    # end1 <- Sys.time()
    # time1 <- end1 - start1
    # # NOte: Block composition c(2,1) is equivalent to the usual log-likelihood
    # start2 <- Sys.time()
    # comp1.NegLogLik <- composite.nll(covarianceMLE$par, residual.GCM.data,block.2,"Matern", m.error.GCM.vec, "NegLogLik") 
    # end2 <- Sys.time()
    # time2 <- end2 - start2
    # #
    # start3 <- Sys.time()
    # comp2.NegLogLik <- composite.nll(covarianceMLE$par, residual.GCM.data,block.3,"Matern", m.error.GCM.vec, "NegLogLik") 
    # #Agrees with real.NegLogLik 
    # end3 <- Sys.time()
    # time3 <- end3 - start3
    # #
    # #save into the table to output:
    # output.lik.pseudo[k,,j] <- c(comp.NegLogLik, comp1.NegLogLik, comp2.NegLogLik)
    # output.time.pseudo[k,,j] <- c(time1, time2, time3)
    #
    # Full log-likelihood calculation:
    dist.mat <- geoDistance(location.GCM,location.GCM) #create the matrix of distances between locations
    Cov.mat <- maternCov(dist.mat, covarianceMLE$par, m.error.GCM.vec) 
    full.neglik.pseudo[k,j] <- 1/2*log(det(Cov.mat)) + 1/2*t(data.GCM)%*%solve(Cov.mat)%*%data.GCM
  }
  
   
  colnames(location.GCM) <- c("lon","lat")
  location.GCM.points <- data.frame(location.GCM)
  attach(location.GCM.points)
  coordinates(location.GCM.points) = ~ lon + lat
  proj4string(location.GCM.points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
  
  krige.output <- krige0(pseudo.GCM ~ 1, sp.pseudo, location.GCM.points, model = maternCov.krige, par.vec = covarianceMLE$par , m.error = m.error.GCM.vec )
  
  full.neglik.pseudo[9,j] <- 1/2*log(det(Cov.mat)) + 1/2*t(krige.output)%*%solve(Cov.mat)%*%krige.output
}



## Code accompanying the tex file: FittingGPs_sphere (Mar. 2017)
#Version 2 includes MLE for Variogram model parameters.

source("compNLL_Func.R")
source("Code/Init.R")

set.seed(1)
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
z.dataframe <- data.frame(Obs$z)
y_obs.sp = SpatialPointsDataFrame(cbind(Obs$x,Obs$y), z.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(y_obs.sp) ~ lon + lat  #coordinate data
summary(y_obs.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 

# 2. Fitted Variogram:
vgm.sample <- variogram(Obs.z ~ 1, data = y_obs.sp, cutoff = 5000, width = 500 )  
plot(vgm.sample)
dev.copy(pdf,'sampleVGM.pdf')
dev.off()
#line of best fit:
# show.vgms()
model.choice1 <- vgm(psill =  9, "Exp", range = 5.5, nugget = 1.5) 
model.choice2 <- vgm(psill =  13, "Exp", range = 11, fit.nuget = FALSE) 
model.choice3 <- vgm(1, "Lin", 0)
#This is singular if we include the nugget effect. 
fit.vgm.exp <- fit.variogram(vgm.sample, model = model.choice3)
plot(vgm.sample, model = fit.vgm.exp)
dev.copy(pdf,'linearVGM.pdf')
dev.off()
#There is a problem here. The choice of variogram model is important because it defines the spatial correlation but I am unable to
#find a suitable variogram model. 

#2. ii) It is also possible (preferable?) to estimate the parmaeters in the variogram model using MLE 
#
locations <- cbind(Obs$x,Obs$y)
n <- length(lm.residuals)
dist.mat <- geoDistance(locations,locations)
param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000) #for validity on the sphere.
ini <- c(0.5,15,20) #how do we find initial values?
mle.output.matern <- optim(ini, matern.neg.lik.calc, data = lm.residuals, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern

# 3. Kriging: Observations
keep <- round(runif(10,1,95)) #index of observation for prediction
data_krige <- cbind(y_obs.sp@coords, y_obs.sp@data$Obs.z) 
colnames(data_krige) <- c("lon","lat","Obs.z")
data.pred <- data_krige[keep,] #locations to predict values
data.pred.points <- data.frame(data.pred)
attach(data.pred.points)
coordinates(data.pred.points) = ~ lon + lat
loc.data.pred <- data.frame(data.pred.points@coords)
attach(loc.data.pred)
coordinates(loc.data.pred) = ~ lon + lat
data.fit <- data_krige[-keep,] #locations and values kept for fitting the GP.
data.fit <- data.frame(data.fit)
attach(data.fit)
coordinates(data.fit) = ~ lon + lat
#define the coordinate systems:
proj4string(loc.data.pred) <- CRS("+init=epsg:4326") #this is a lon,lat grid projection of the Earth (I believe.)
proj4string(data.fit) <- CRS("+init=epsg:4326") #this is a lon,lat grid projection of the Earth (I believe.)
m2.krige <- krige(Obs.z ~ 1, data.fit, loc.data.pred, model = model.choice2) 
m2.krige #output the results.

# 4. Composite likelihood calculation: Observations
pred.data <- cbind(data.pred[1:5,1:2], t(t(m2.krige@data$var1.pred[1:5])))
names(pred.data) <- c("x.loc","y.loc","val")
#Note we will get a singular covairance matrix (Cov.mat) when we include the one location that yields a unique predicted value
#because of linear dependence. 

maternParam <- c(3,1,1)
dist.mat <- dist2(pred.data[,1:2],pred.data[,1:2],method="euclidean")
Cov.mat <- maternCov(dist.mat,maternParam)  
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

## calculating the composite likelihood:
maternParam <- c(3,1,1) #parameters - how should we select the variogram model?
#remove NA for now:
Model_df.full <- Model_df[complete.cases(Model_df),]
n_ind <- 10 #we don't want to compute this over the entire grid.
index <- round(runif(n_ind,1,dim(Model_df.full)[1]))
ptm1 <- proc.time()
dist.mat <- dist2(Model_df.full[index,2:3],Model_df.full[index,2:3],method="euclidean")
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





###################################
# testing variogram fitting and empirical:
sample.vgm.func <- function(cutoff_, width_){
  vgm.sample <- variogram(Obs.z ~ 1, data = y_obs.sp, cutoff = cutoff_, width = width_ )  
  plot(vgm.sample)
}
#1 . Changing the width
sample.vgm.func(5000,500)
sample.vgm.func(5000,1000)
sample.vgm.func(5000,100)

#2. Changing cutoff (keep width constant relative to cutoff)
sample.vgm.func(20000,2000)
sample.vgm.func(5000,500)
sample.vgm.func(500,50)
sample.vgm.func(100,10)

#3. Changing cutoff (keep width absolutely the same)
sample.vgm.func(20000,500)
sample.vgm.func(5000,500)

#4. Refining semivariogram plot:
sample.vgm.func(15000,1500)
dev.copy(pdf,'linearVGM.pdf')
dev.off()
sample.vgm.func(15000,500)
dev.copy(pdf,'linearVGM.pdf')
dev.off()
sample.vgm.func(15000,3000)
dev.copy(pdf,'linearVGM.pdf')
dev.off()

#it appears as though the space becomes less correlated until around 7500 where it flattens a bit, then starts decreasing again.
#One possibility is that the decrease at the end is because we are moving back closer to the points in the space (because spherical coordinates)
#Surely this is common? How do we deal with this? Surely the distance should be the shortest line between two points. Maybe there is some other
#reason for this parabola?
sample.vgm.func(8000,1000)
dev.copy(pdf,'sampleVGMIdeal.pdf')
dev.off()
#This gives a "standard" empirical variogram plot. 
#output for plots:

sample.vgm.func(15000,1500)
dev.copy(pdf,'sampleVGM1.pdf')
dev.off()
sample.vgm.func(15000,500)
dev.copy(pdf,'sampleVGM2.pdf')
dev.off()
sample.vgm.func(5000,500)
dev.copy(pdf,'sampleVGM3.pdf')
dev.off()






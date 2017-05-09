### Can we have a variable nugget term?
#
Obs_set='P3+_SST_anom'
obs = c('lambda_100', 'lambda_90', 'lambda_80', 'lambda_70', 'lambda_60', 'lambda_50', 'lambda_40', 'lambda_30', 'lambda_20', 'lambda_10')
j <- 1 #one time instance. 
# READ OBSERVATIONS
Obs <- read.table(paste0("Observation_data/",Obs_set,"/",obs[j],".txt"),header=T)

#create a spatial points data frame object:
locations <- cbind(Obs$x,Obs$y)
z.dataframe <- data.frame(Obs$z)
y_obs.sp = SpatialPointsDataFrame(locations, z.dataframe, coords.nrs = numeric(0), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#CRS is the coordinate reference system, this is the reference system of the location data. In our case this is longitude and latitude. 
coordinates(y_obs.sp) ~ lon + lat  #coordinate data
summary(y_obs.sp)  #is.projected = FALSE as expected - the locations are not projected onto a lower dimensional space. 


####
n <- length(Obs$z)
dist.mat <- geoDistance(locations[1:10,],locations[1:10,])
param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000) #for validity on the sphere.
ini <- c(0.5,15,20) #how do we find initial values?
mle.output.matern <- optim(ini, matern.neg.lik.calc, data = Obs$z, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B")
mle.output.matern
## Code to Load the observations and GCM data:
# 05.05.17 Phillip Paine (V1)

# ---------- Pre-Amble: -----------------------------
# Set working directory to DMC/Scripts.
# Load necessary functions:
source("Code/Init.R")

# Consider a subset of the full GCM dataset for computation time:
jump <- 270 #how regularly we want to retain locations in the grid
start <- runif(1,1,jump)
index <- seq(start, 27000, by = jump)
n_ind = length(index)
# ----------------------------------------------------
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

# ii) Climate Models:
#
Model_set='CO2_anom'
models = c('tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti')
# Specify grid dimensions for model data : # (lon_min,lon_max,lon_int,lat_min,lat_max,lat_int)
model_grid <- c(-180,178.75,1.25,-89.375,89.375,1.25) # HadCM3 ocean temps, shifted 1 hemisphere

GCM.data <- sapply(1:8, function(j){
    Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[j],".txt"),
                              indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                              model_grid, plot_it=F) 

   #remove NA:
   Model_df.full <- Model_df[complete.cases(Model_df),]
   location.GCM <- Model_df.full[index,2:3]
   z.GCM <- array(Model_df.full[index,5], dim=c(length(index),1))
   return(cbind(location.GCM,z.GCM))
} )





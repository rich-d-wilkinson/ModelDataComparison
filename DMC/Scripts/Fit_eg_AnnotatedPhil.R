start_time <- proc.time()
source("Code/Init.R")
set.seed(1)

# Example script

# ------------------------------------------
# INPUT DATA AND SETTINGS

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

# Observations
# -------------
# Observations need to have points in a list of 4 columns with headers "x,y,z,std" = lon, lat, observation, uncertainty
# Observations are stored in directory Observation_Data/<Obs_set> with filenames <obs>.txt
Obs_set='P3+_SST_anom'
obs = c('lambda_100', 'lambda_90', 'lambda_80', 'lambda_70', 'lambda_60', 'lambda_50', 'lambda_40', 'lambda_30', 'lambda_20', 'lambda_10')

# Mesh settings
# -------------
# Choose egg or sphere mesh
#manifold = "Egg" 
 manifold = "Sphere"

# Matlab is called for egg meshes to generate new mesh with specified anisotropy = "stretch"  (MATLAB file: Generate_Mesh_egg.m)
# stretch = sqrt of equatorial radius, default 0.05 (polar radius = 1)
# Although stretch=1 should give a sphere, this doesn't actually seem to work for a reason I haven't worked out, use sphere option.
#stretch <- 0.05 # Egg

# Sphere - writes different set of mesh files (MATLAB file: Generate_Mesh.m)
# For sphere option, don't need to vary radius, use smooth to vary region of influence, so this mesh is not updated.
stretch <- 1 # Sphere

# Calculation settings
# -------------------
# smoothness parameter, i.e. radius of influence. So higher is smoother. Default is smooth = 0.3.
smooth = 0.3

# These lines retained for reference only - context within greater system.
fit = "Obs" # "Sim" "Obs" or "Model"

# Output
# -------------------
mesh <- paste0(manifold,"-",stretch)
outdir <- paste0("Output/O-",Obs_set,"_M-",Model_set,"_ME-",mesh,"_SM-",smooth)

# ------------------------------------------
# END - INPUT DATA AND SETTINGS
# ------------------------------------------

# INITIALISATION

r = c(sqrt(stretch), sqrt(stretch),1) 
model_grid_rad <- model_grid*pi/180.
Lmodel_summary <- matrix(data=NA,nrow=length(obs),ncol=length(models))

# SET UP GLOBAL MESH

#Adjust radius in 3d space accordingly
if (manifold == "Egg") {
# Call Matlab to generate egg mesh with specified stretch value
if (.Platform$OS.type == "windows"){
    Matlab$startServer()
    matlab_server <<- Matlab(remote=F)
    Sys.sleep(10)
    open(matlab_server)
  } else if (.Platform$OS.type == "unix"){
    Matlab$startServer()
    matlab_server <<- Matlab(remote=F)
    Sys.sleep(10)
    open(matlab_server)
  }
  setVariable(matlab_server,stretch=stretch)
  cat("Generating egg mesh in Matlab",sep="\n")
  evaluate(matlab_server, paste("cd Code;Generate_Mesh_egg(stretch);"))
  close(matlab_server)
  cat("Finished generating egg mesh in MATLAB",sep="\n")
  
  Sphere_triang <- Load_egg_mesh()
} else {
  # mesh should not need updating for sphere as it does not change
  Sphere_triang <- Load_sphere_mesh()
}
p_lon_lat <- Convert_to_lon_lat(Sphere_triang@pars$p,1)
# ------------------------------------------

# READ OBSERVATIONS AND FIT SURFACE

for (j in 1:length(obs)) {  #fit for each time-point (each of lambda_100 etc...)
print(paste0("Fitting observations for ",obs[j]),quote=FALSE)
flush.console()

# READ OBSERVATIONS
Obs <- read.table(paste0("Observation_data/",Obs_set,"/",obs[j],".txt"),header=T)
#Convert to radians
Obs$x <- (Obs$x)*2*pi/360
Obs$y <- (Obs$y)*2*pi/360

# Organise observations
n_obs = nrow(Obs)      
locs = list(Obs$x,Obs$y)  
locs_xyz <- Convert_to_xyz(cbind(locs[[1]],locs[[2]]),r) #convert to cartesian coordinates
	
y_obs <- data.frame(x=locs_xyz[,1], y = locs_xyz[,2], z = locs_xyz[,3], val = as.vector(Obs$z),lon=locs[[1]],lat=locs[[2]])
Qobs <- sparseMatrix(i=1:n_obs,j=1:n_obs,x=1/Obs$std^2)
	
# Find observation locations in mesh
C_full <- FindC_sphere2(Sphere_triang@pars$p,Sphere_triang@pars$t,locs,r) 
# This is an updated, cool version of how we find points on the surface
# I got around the problem by projecting twice, first using the real z-axis and then switching the z-axis with the y-axis and taking
# the union of the two results. Works great.

# GENERATE FIT SURFACE
# Q is the precision matrix
Q <- Prec_from_SPDE_wrapper(M = Sphere_triang@pars$M,
                            K = Sphere_triang@pars$K,
                            intrinsic = 0, 
                            nu = 1,
                            desired_prec = 0.01*Imat(Sphere_triang@n),
                            l = smooth)
			    
# Structure
Inf_data <- list(C_full = C_full,
                 Qobs = Qobs,
                 y_obs = y_obs,
                 Q = Q)

# Now infer
Inf_data <- Infer(Inf_data)
#we output the mean, standard deviation of diagonal elements and precision matrix of weight given obs - which has dimension 2964x2964 because
# w | Y ~ N_(2964)(mean, covariance) - that is it describes the distribution on the grid. A draw from this is a draw from the discretised field. 

# Configure precision matrix
# Simulate from using this precision matrix
GMRF <- GMRF_init(as.matrix(rep(0,Sphere_triang@n)),Q,intrinsic = 0,perm=T,diag_offset=0.0001)
x_sim <- SampleGMRF(GMRF,reps=1,use_perm=T)

# EXTRACT RESULLTS
Res <- Inf_data$Res
Res <- cbind(Res,data.frame(x=p_lon_lat[,1],y=p_lon_lat[,2],val=x_sim))

# ------------------------------------------
# Read models and compare with observations

for (i in 1:length(models)) {
print(paste0("Reading model ",models[i]),quote=FALSE)
flush.console()

# Read_model_data generates the standard UM ocean grid and reads data into complete grid as vector
# with NAs for land
Model_df <- Read_model_data(modelpath = paste0("Model_data/",Model_set,"/",models[i],".txt"),
                            indexpath = paste0("Model_data/",Model_set,"/mask.txt"),
                            model_grid_rad, plot_it=T)
# apply mask to results - this needs to be done once for each set of observations
if (i == 1) {
   # Apply mask to FE results: Res
   # This could actually be done right at the top using p_lon_lat when the mesh is generated
   # Mask would have to be read separatately.
   # Slightly nervous that some processes seem to change the order of the data though.
   print("Mapping mask locations onto FE mesh",quote=FALSE)
   mask <- data.frame(Model_df$x,Model_df$y,Model_df$mask)
   names(mask) = c("grid_x","grid_y","mask")
   # Finds indices of Res (on FE grid) that are within mask (on model grid)
   # Find model grid cell for each element of Res
   Res$grid_x <- round((Res$x-model_grid_rad[1])/model_grid_rad[3])*model_grid_rad[3]+model_grid_rad[1]
   Res$grid_y <- round((Res$y-model_grid_rad[4])/model_grid_rad[6])*model_grid_rad[6]+model_grid_rad[4]
   # Round x and y values to pre-empt merge failures
   # This method might be flaky, time will tell! Merge needs exact values to work.
   Res$grid_x  <- round(Res$grid_x,digits=4)
   Res$grid_y  <- round(Res$grid_y,digits=4)
   mask$grid_x <- round(mask$grid_x,digits=4)
   mask$grid_y <- round(mask$grid_y,digits=4)
   # Lookup mask values for Res
   Res <- merge(Res,mask,by=c("grid_x","grid_y"),all.x=T,sort=F)
   mask_ind <- which(Res$mask == 1)
   # Locate any points that have failed and print out
   mask_fail_ind <- which(is.na(Res$mask))
   if (length(mask_fail_ind) != 0) {
      mask_fail <- cbind(Res$grid_x[mask_fail_ind]*180/pi,Res$grid_y[mask_fail_ind]*180/pi)
      print("Mask location failures",quote=FALSE)
      print(mask_fail)
      # don't get rid of the failures, just set the NAs to 0 so they will be ignored
      Res$mask[mask_fail_ind] <- 0
      }
} # end mask

# This function removes all the NAs from the vector
Model_df <- Model_df[complete.cases(Model_df),]

# Add the models I want to compare to the results by regressing onto the sphere. Here I just add Model_df
Res <- Add_model(Res,Model_df,Sphere_triang,name="Model",r) 

# ------------------------------------------

# CALCULATE LOG LIKELIHOOD BOTH INCLUDING AND EXCLUDING MASK AREAS

# This won't work because the zero at land really gives the field a low probability (not smooth at all)

Lmodel <- logLik_prop(Res$Model,Res$mean,Inf_data$Qtot)

Qpost_mask <- chol2inv(chol(chol2inv(chol(Inf_data$Qtot))[mask_ind,mask_ind]))
Lmodel_mask <- logLik_prop(Res$Model[mask_ind],
                              Res$mean[mask_ind],
                              Qpost_mask)

Lmodel_summary[j,i] <- Lmodel_mask

Null_Lmodel_mask <- logLik_prop(0*Res$Model[mask_ind],
                              Res$mean[mask_ind],
                                   Qpost_mask)
} # model loop/ i

# ------------------------------------------

# OUTPUT
print("Generating output",quote=FALSE)

# FE output for each set of observations

# This writes the FE output surface for the mask region only to file - this is probably not much use to most end users
if (!(file.exists(outdir))) {
   dir.create(outdir)
   }
# write.table(Res[mask_ind,],file = paste0("outdir,"/",obs[j],"_FE-output.txt"))

# Grid data for plotting/output for each set of observations

detail=200 # number of points in each grid direction
GG <- Triang_to_grid(Res,"mean",detail)
GGstd <- Triang_to_grid(Res,"std",detail)
GGmask <- Triang_to_grid(Res,"mask",detail)
GGtrue <- Triang_to_grid(Res,"val",detail)
GGcv <- GG; GGcv$z <- GG$z/GGstd$z
Gmodel <- Triang_to_grid(Res,"Model",detail)

# Create data frame with gridded results and remove NaNs
Res_grid <- data.frame(GG,GGstd,GGmask)
Res_grid <- Res_grid[complete.cases(Res_grid),]

# This writes the rectangular gridded output of mean, std and mask to file for plotting
write.table(Res_grid,file = paste0(outdir,"/",obs[j],"_grid_output.txt"))
} # fit loop/ j

# Results matrix
Results <- as.data.frame(Lmodel_summary)
colnames(Results) <- models
rownames(Results) <- obs

print("Results matrix",quote=FALSE)
print(Results)
flush.console()

write.table(Results,file = paste0(outdir,"/Results_matrix.txt"))
end_time <- proc.time()
print("Total time",quote=FALSE)
print(end_time - start_time)



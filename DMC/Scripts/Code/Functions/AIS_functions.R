Load_shapefiles <- function() {
  
  # load shapefiles
  grounding <- readShapeSpatial("./AIS_Data/Coastline/moa_groundingline.shp")
  grounding_fort <- fortify(grounding)
  grounding_fort$long <- grounding_fort$long/1000
  grounding_fort$lat <- grounding_fort$lat/1000
  grounding_sub <- grounding_fort[seq(1,nrow(grounding_fort),100),]
  
  coast <- readShapeSpatial("./AIS_Data/Coastline/moa_coastline.shp")
  coast_fort <- fortify(coast)
  coast_fort$long <- coast_fort$long/1000
  coast_fort$lat <- coast_fort$lat/1000
  coast_sub <- coast_fort[seq(1,nrow(coast_fort),100),]
  
  Islands <- read.table("./AIS_Data/Coastline/all_islands.txt",col.names=c("id","x","y"))
  rm(grounding_fort)
  
  
  shapefiles <- list(grounding_sub = grounding_sub,
                     coast_sub = coast_sub,
                     Islands = Islands)
  return(shapefiles)
}


Load_meshes <- function() {
  p_ICE  <- round(as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_pICE.csv',sep=",")))
  tri_ICE <- as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_tICE.csv',sep=","))
  M_ICE <-  as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_MICE.csv',sep=","))
  K_ICE <-  as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_KICE.csv',sep=","))
  
  p_SURF  <- round(as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_pSURF.csv',sep=",")))
  tri_SURF <- as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_tSURF.csv',sep=","))
  M_SURF <-  as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_MSURF.csv',sep=","))
  K_SURF <-  as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_KSURF.csv',sep=","))
  
  
  p_GIA <- round(as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_pGIA.csv',sep=",")))
  tri_GIA <- as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_tGIA.csv',sep=","))
  M_GIA <- as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_MGIA.csv',sep=","))
  K_GIA <- as.matrix(read.table('./AIS_ST_study/Triangulation_MATLAB/AIS_triangulation_KGIA.csv',sep=",")) 
  
  # put into sparse matrix format
  n_ICE <- nrow(p_ICE)
  n_SURF <- nrow(p_SURF)
  n_GIA <- nrow(p_GIA)
  
  M_ICE <- sparseMatrix(i=M_ICE[,1],j=M_ICE[,2],x=M_ICE[,3],dims=c(n_ICE,n_ICE))
  K_ICE <- sparseMatrix(i=K_ICE[,1],j=K_ICE[,2],x=K_ICE[,3],dims=c(n_ICE,n_ICE))
  M_SURF <- sparseMatrix(i=M_SURF[,1],j=M_SURF[,2],x=M_SURF[,3],dims=c(n_SURF,n_SURF))
  K_SURF <- sparseMatrix(i=K_SURF[,1],j=K_SURF[,2],x=K_SURF[,3],dims=c(n_SURF,n_SURF))
  M_GIA <- sparseMatrix(i=M_GIA[,1],j=M_GIA[,2],x=M_GIA[,3],dims=c(n_GIA,n_GIA))
  K_GIA <- sparseMatrix(i=K_GIA[,1],j=K_GIA[,2],x=K_GIA[,3],dims=c(n_GIA,n_GIA))
  
  Meshes <- list(ICE = initFEbasis(p=p_ICE, t = tri_ICE, M = M_ICE, K = K_ICE),
                 SURF = initFEbasis(p=p_SURF, t = tri_SURF, M = M_SURF, K = K_SURF),
                 GIA = initFEbasis(p=p_GIA, t = tri_GIA, M = M_GIA, K = K_GIA))
  return(Meshes)
  
}


Extract_RACMO_density <- function(shapefiles,load=F) {
  if (load==F) {
    RACMO <- read.table('./AIS_Data/RACMO_by_year/Racmo_prior_ann_anom_2003-2009_wrt20yr.dat',skip=2,col.names=c('x','y','2003','2004','2005','2006','2007','2008','2009'))
    load("./AIS_Data/Densities/rho_rates.RData")
    names(rho) <- c("x","y","rho")
    RACMO <- merge(RACMO,rho,all.x=T)
    missing_rho <- which(is.na(RACMO$rho)) 
    for(i in missing_rho) {
      sub_region <- subset(RACMO, (x < RACMO$x[i]+100) & (y < RACMO$y[i]+100) &
                             (x > RACMO$x[i]-200) & (y > RACMO$y[i]-200))
      RACMO$rho[i] <- mean(sub_region$rho,na.rm=T)
    }
    
    # Check whether within grounding line or not
    RACMO$in_land <- pnt.in.poly(cbind(RACMO$x,RACMO$y),shapefiles$grounding_sub[,1:2])$pip
    RACMO$IslandNum <- 0
    for( i in unique(shapefiles$Islands$id)) {
      my_sub <- subset(shapefiles$Islands,id==i)
      RACMO$in_land <- RACMO$in_land | pnt.in.poly(cbind(RACMO$x,RACMO$y),my_sub[,2:3])$pip
      myind <- which(pnt.in.poly(cbind(RACMO$x,RACMO$y),my_sub[,2:3])$pip == 1)
      RACMO$IslandNum[myind] <- i
    }
    save(RACMO,file='./AIS_Data/RACMO_by_year/RACMO_density_merged.Rdata')
  } else {
    load('./AIS_Data/RACMO_by_year/RACMO_density_merged.Rdata')
  }
  return(RACMO)
}



RACMO_basis <- function(RACMO,nEOF,nPKC,T,load=F) {
  
  if (load==F) {
    RACMO$rho <- RACMO$rho/1000 # convert to Mt per km2m
    RACMO[,3:(3+T-1)] <- RACMO[,3:(3+T-1)]/1000 # Convert to mweq
    RACMO_sub <- t(subset(RACMO,select = -c(x,y,rho,in_land,IslandNum)))
    RACMO_XY <- subset(RACMO,select = c(x,y))
    
    # Correct for density
    cat("Using interpolated density to correct RACMO",sep="\n")
    RACMO_sub <- (RACMO_sub) /(RACMO$rho) # Amplify to convert from mweq to actual height (density in Mt per km2m is in correct units)
    S <- Get_RACMO_basis(nEOF=6,nPKC=6,T=7)
    save(RACMO_XY,RACMO_sub,S,file='./AIS_Data/RACMO_by_year/RACMO_basis.Rdata')
  } else {
    load('./AIS_Data/RACMO_by_year/RACMO_basis.Rdata')
  }
  return(S)
}


LandAttribute <- function(df,shapefiles) {
  df$in_land <- pnt.in.poly(cbind(df$x,df$y),shapefiles$grounding_sub[,1:2])$pip
  df$IslandNum <- 0
  for( i in unique(shapefiles$Islands$id)) {
    my_sub <- subset(shapefiles$Islands,id==i)
    df$in_land <- df$in_land | pnt.in.poly(cbind(df$x,df$y),my_sub[,2:3])$pip
    myind <- which(pnt.in.poly(cbind(df$x,df$y),my_sub[,2:3])$pip == 1)
    df$IslandNum[myind] <- i
  }
  df$in_coast <- pnt.in.poly(cbind(df$x,df$y),shapefiles$coast_sub[,1:2])$pip
  df$floating <- df$in_coast  & !(df$in_land)
  df$in_coast[which(df$in_land == 1)] = 1  # For islands
  return(df)
}

Tessellate <- function(Meshes,load=T) {
  if(!load) {
    Tessellation <- list()
    for (meshname in names(Meshes)) {
      Mesh <- Meshes[[meshname]]
      Voronoi <- deldir(Mesh@pars$p[,1],
                        Mesh@pars$p[,2],
                        plotit='T',
                        sort=F,
                        rw=c(min(Mesh@pars$p[,1])-0.00001,
                             max(Mesh@pars$p[,1])+0.00001,
                             min(Mesh@pars$p[,2])-0.00001,
                             max(Mesh@pars$p[,2])+.00001))
      pol <- PolygonfromVoronoi(Voronoi,Mesh@pars$p)
      Tessellation[[meshname]] <- list(Voronoi = Voronoi, pol = pol)
    }
    save(Tessellation,file="./AIS_ST_study/Precomputed_ST/Voronoi_polygons.RData") 
    return(Tessellation)
  } else {
    load("./AIS_ST_study/Precomputed_ST/Voronoi_polygons.RData")  
  }
}


Attribute_variables <- function(Meshes,shapefiles,load=F) {
  if(!load) {
    Basins_ingo <- read.table("./AIS_Data/Basins/ingo_basins_1km_xyz.txt",col.names=c("y","x","Basin_ingo"))
    Basins_rignot11 <- read.table("./AIS_Data/Basins/Basins_Rignot_11.dat",col.names=c("x","y","z"))
    IV <- read.table("./AIS_Data/Velocities/velocities_full.dat",col.names=c("x","y","z"))
    load("./AIS_Data/Velocities/velocities_balance.RData")
    
    for (meshname in names(Meshes)) {
      df <- data.frame(x = Meshes[[meshname]]@pars$p[,1],
                       y = Meshes[[meshname]]@pars$p[,2],
                       n = 1:Meshes[[meshname]]@n)
      
      # Basins - Ingo
      df <- merge(df,Basins_ingo,by=c("x","y"),all.x=T)
      df <- arrange(df,desc(-n))
      df$Basin_ingo[is.na(df$Basin_ingo)] = 86
      df$Basin_ingo[df$Basin_ingo == 0] = 86
      # Manual fixing:
      if (any (df$x == -982 & df$y==-146)) df[df$x==-982 & df$y==-146,]$Basin_ingo <- 301
      if (any (df$x == -1506 & df$y == 105)) df[df$x==-1506 & df$y==105,]$Basin_ingo <- 324
      
      # Basins - Rignot11
      Basins_rignot11_cont <- function(s) { return(nn_grid_interp(s,df=Basins_rignot11,delta=5,miss_value=NA)) }
      df$Basins_rignot11  <- Basins_rignot11_cont(cbind(df$x,df$y))
      df$Basins_rignot11[is.na(df$Basins_rignot11)] = 0
      
      # In land? In sea? In island?
      df <- LandAttribute(df,shapefiles)
      
      # ICE velocities
      cat('Regressing Ice velocities ...',sep="\n");   flush.console()
      IV_cont <- function(s) { return(nn_grid_interp(s,df=IV,delta=10,miss_value=NA)) }
      df$Velocity <- IV_cont(cbind(df$x,df$y))
      
      # Balance velocities
      cat('Regressing Balance Ice velocities ...',sep="\n");   flush.console()
      BalVel_cont <- function(s) {return(nn_grid_interp(s,df=Bal_vel,delta=5,miss_value=NA)) }
      df$Velocityb <- BalVel_cont(cbind(df$x,df$y))
      
      # Join together
      cat('Replacing missing INSAR with Balance velocities ...',sep="\n");   flush.console()
      df$Velocity[which(is.na(df$Velocity))] <- df$Velocityb[which(is.na(df$Velocity))]    # Now find pmax
      cat('WARNING: Where we do not have velocity readings we are using 5m/a for prior on land ...',sep="\n");   flush.console()
      df$Velocity[is.na(df$Velocity) & (df$in_coast == 1) ] = 5
      df$Velocity[is.na(df$Velocity)] = 0
      
      Meshes[[meshname]]@pars$vars <- df
    }
    rm("Basins_ingo")
    rm("Basins_rignot11")
    save(Meshes,file='./AIS_ST_study/Precomputed_ST/Meshes.RData')  
  } else {
    load(file='./AIS_ST_study/Precomputed_ST/Meshes.RData')   
  }
  return(Meshes)
}


Preprocess_altimetry <- function(shapefiles,averaged = T, load=F) {
  
  if (!load) {
    # Load ICESAT
    if(doICESAT) {
      cat('Assembling ICESAT data ...',sep="\n"); flush.console()
      ICESAT_by_year <- NULL
      for ( year in c(2003:2009) ) {
        ICESAT_by_year <- rbind(ICESAT_by_year,cbind(read.table(paste('./AIS_Data/ICESAT_by_year/Annual_Icesat_dhdts_',
                                                                      year,".dat",sep=""),skip=1),
                                                     year=year)) 
      }
      names(ICESAT_by_year) <- c("x","y","z","std","year")
      ICESAT_by_year <- subset(ICESAT_by_year,std < 0.4)
      ICESAT_by_year$n = 1:nrow(ICESAT_by_year)
      ICESAT_by_year$t <-  ICESAT_by_year$year - ICESAT_by_year$year[1]
      save(ICESAT_by_year,file='./AIS_ST_study/Precomputed_ST/ICESAT_by_year.RData')  
      
      ICESAT_by_year <- Average_observations(ICESAT_by_year,box_size=20)
      save(ICESAT_by_year,file='./AIS_ST_study/Precomputed_ST/ICESAT_by_year_averaged.RData')  
    } 
    
    
    
    # Do ENVISAT
    if(doENVISAT) {
      ENVISAT_grouped <- read.table('./AIS_Data/Envisat_by_year/Envisat_trends_2003_2009.dat',skip=1)
      ENVISAT_by_year <- NULL
      for (t in t_axis) {
        ENVISAT_by_year <- rbind(ENVISAT_by_year,data.frame(x=ENVISAT_grouped[,1],y=ENVISAT_grouped[,2],
                                                            z=ENVISAT_grouped[,3+t], 
                                                            std=ENVISAT_grouped[,3+t+T],
                                                            t=t,year=2003+t))    
      }
      save(ENVISAT_by_year,file='./AIS_ST_study/Precomputed_ST/ENVISAT_by_year.RData')  
      ENVISAT_by_year <- Average_observations(ENVISAT_by_year,box_size=20)
      # Some of these are outside land, we need to filter
      ENVISAT_by_year <- LandAttribute(shapefiles,ENVISAT_by_year)
      ENVISAT_by_year <- subset(ENVISAT_by_year,in_coast==1)
      save(ENVISAT_by_year,file='./AIS_ST_study/Precomputed_ST/ENVISAT_by_year_averaged.RData')  
    }
  }
  
  
  if(!averaged) load(file='./AIS_ST_study/Precomputed_ST/ICESAT_by_year.RData')  
  if(averaged) load(file='./AIS_ST_study/Precomputed_ST/ICESAT_by_year_averaged.RData')  
  if(!averaged) load(file='./AIS_ST_study/Precomputed_ST/ENVISAT_by_year.RData')  
  if(averaged) load(file='./AIS_ST_study/Precomputed_ST/ENVISAT_by_year_averaged.RData')  
  return(list(ICESAT = ICESAT_by_year, 
              ENVISAT = ENVISAT_by_year))
}



Preprocess_GRACE <- function(load) {
  if(!load) {
    source('./AIS_ST_study/GRACE_preprocessing.R')
    print(ggplot() + geom_polygon(data=subset(GRACE_big,t==6),aes(x=x,y=y,fill=cmweq, group=id)) + scale_fill_gradient2(low = "red",mid="white",high="blue") + 
            geom_path(data =  grounding_sub,aes(long,lat),colour="black"))
  } else {
    load(file='./AIS_ST_study/Precomputed_ST/GRACE_aggregated.RData')  
    GRACE$z <-  as.numeric(GRACE$cmweq)*0.01*GRACE$area2    #Convert to Mt
    GRACE_big$z <-  as.numeric(GRACE_big$cmweq)*0.01*GRACE_big$area2    #Convert to Mt
    GRACE$std <-  as.numeric(GRACE$std)*0.01*GRACE$area2    #Convert to Mt
    GRACE_big$std <-  as.numeric(GRACE_big$std)*0.01*GRACE_big$area2    #Convert to Mt
  }
  
  Smooth_mat_a <- function(alpha,mGRACE) {
    return (sparseMatrix(1:mGRACE,1:mGRACE,x=1) + alpha * (D2 - D1))
  }
  alpha0 = 0.107
  gamma0 = 11.1
  
  return(list(GRACE=GRACE,GRACE_big=GRACE_big, GRACE_poly = GRACE_poly, Smooth_mat_a = Smooth_mat_a, alpha0 = alpha0, gamma0 = gamma0))  
}

Preprocess_GPS <- function() {
  # Load GPS
  GPS_per_year <-  read.table("./WAIS_Data/GPS rates/GPS_rates_WAIS_2003-2009_Thomas_King.dat",header=T)
  names(GPS_per_year) <- c("x","y","z","std")
  GPS_per_year$z <- GPS_per_year$z/1000     # convert to m
  GPS_per_year$std <- GPS_per_year$std/1000 # convert to m
  # REMOVE ANOMALOUS NEGATIVE OBS
  GPS_per_year <- GPS_per_year[-which.min(GPS_per_year$z),]
  
  GPS <- NULL
  for (t in t_axis) {
    GPS <- rbind(GPS,cbind(GPS_per_year,data.frame(t=rep(t,nrow(GPS_per_year)))))
  }
  
  # SINCE THIS MEAN WAS ESTIMATED OVER 7 YEARS I SHOULD IMITATE A TEMPORAL PROGRESSION BY AMPLIFYING IT PER YEAR
  # i.e. I MULTIPLY IT BY sqrt(N)
  GPS$std <- GPS$std*sqrt(T)
  
  if (!useGPS) {
    GPS$z <- 0
    GPS$std <- 10
  }
  return(list(GPS = GPS))
}



Find_Smooth_mat <- function(GRACE,alpha0){
  mGRACE <- nrow(GRACE)
  X <- pdist(cbind(GRACE$x,GRACE$y),cbind(GRACE$x,GRACE$y))
  X2 <- X*(X<450)
  GRACE_neighb <- neighb_from_prec(X2)
  GRACE$est_smooth <- 0
  Smooth_mat <- matrix(0,mGRACE,mGRACE)
  # Numerical diffusion from http://pauli.uni-muenster.de/tp/fileadmin/lehre/NumMethoden/WS0910/ScriptPDE/Heat.pdf
  for(i in 1:mGRACE) {
    nn <- length(GRACE[GRACE_neighb[[i]],]$x)
    Smooth_mat[i,i] <- (1-nn*alpha0)
    Smooth_mat[i,GRACE_neighb[[i]]] <- alpha0
  }
  D2 <- Smooth_mat
  D2[D2 > 0] <- 1
  D2 <- D2 - diag(diag(D2))
  D1 <- diag(rowSums(D2))
  Smooth_mat_a <- function(alpha) {
    return (sparseMatrix(1:mGRACE,1:mGRACE,x=1) + alpha * (D2 - D1))
  }
  P <- Smooth_mat_a(alpha0)
  PNA <- P
  PNA[PNA == 0] <- NA
  dense_mascons <- which(apply(PNA,1,min,na.rm=T) < 0.1)
  # Where mascons are dense assume they read averages
  for (i in dense_mascons) {
    P[i,which(P[i,] != 0) ] = 1/ length(which(P[i,] != 0))
  }
  return(P)
}


Find_C_GRACE <- function(Gravimetry,Meshes,consts,RACMO,S,load=T) {
  if (!load) {
    C_GRACE <- list()
    mGRACE_1yr <- nrow(subset(Gravimetry$GRACE,t==0))
    P <- Find_Smooth_mat(subset(Gravimetry$GRACE,t==0),Gravimetry$alpha0)
    C_GRACE$GIA <- consts$rho_ROCK* P %*% FindC_polyaverage(Meshes[["GIA"]]@pars$p,
                                                            Meshes[["GIA"]]@pars$t,
                                                            Gravimetry$GRACE_poly[1:mGRACE_1yr],
                                                            plotit=F,
                                                            method="C",
                                                            ds=400)   
    C_GRACE$ICE <- consts$rho_SURF_ICE_everywhere* P %*% FindC_polyaverage(Meshes[["ICE"]]@pars$p,
                                                                           Meshes[["ICE"]]@pars$t,
                                                                           Gravimetry$GRACE_poly[1:mGRACE_1yr],
                                                                           plotit=F,
                                                                           method="C",
                                                                           ds=400)    
    C_GRACE$ICE[,which(!Meshes[["ICE"]]@pars$vars$in_land)] = 0
    
    RACMO_XY <- subset(RACMO,select = c(x,y))
    C_GRACE$SURF1 <- P %*% FindC_polyaverage_EOF(X = cbind(RACMO_XY,S),
                                                 polygons = Gravimetry$GRACE_poly[1:mGRACE_1yr],
                                                 mulfun=RACMO$rho*100*as.numeric(RACMO$in_land))    
    
    
    C_GRACE$SURF2 <- P %*% FindC_polyaverage(Meshes[["SURF"]]@pars$p,
                                             Meshes[["SURF"]]@pars$t,
                                             Gravimetry$GRACE_poly[1:mGRACE_1yr],
                                             plotit=F,
                                             method="C",
                                             ds=400,
                                             mulfun=RACMO$rho*100*as.numeric(RACMO$in_land))    
    C_GRACE$SURF <- cBind(C_GRACE$SURF1,C_GRACE$SURF2)
    C_GRACE$FIRN <- P %*% Zeromat(mGRACE_1yr,Meshes[["SURF"]]@n + ncol(S))
    save(C_GRACE,file='./AIS_ST_study/Precomputed_ST/C_GRACE.Rdata')
  } else {
    load(file='./AIS_ST_study/Precomputed_ST/C_GRACE.Rdata')
  }
  return(C_GRACE)
}




ENVISAT_spat_res_structure <- function(All_data,removeS = T,removeFE = T) {
  Qobs = All_data$Qobs
  y_tot = All_data$y_tot
  C_full = All_data$C_full
  Q_full = All_data$Q_full
  S = All_data$S
  Meshes = All_data$Meshes
  shapefiles = All_data$shapefiles
  
  fns_names <- paste("S",(1:dim(S)[2]),sep="")
  x_rep <- data.frame(n = 1:nrow(Q_full), 
                      proc = c(rep(
                        c(rep("ICE",n_ICE),
                          fns_names,rep("SURF",Meshes[["SURF"]]@n),
                          fns_names,rep("FIRN",Meshes[["SURF"]]@n)),T),
                               c(rep("ICEa",n_ICE)),c(rep("ICEb",n_ICE))))
  remove_fns = 1:dim(S)[2]
  S <- as.matrix(S[,-remove_fns])
  
  indices_to_remove=NULL
  if(removeS) {
    indices_to_remove <- c(indices_to_remove,which(x_rep$proc %in% c(fns_names[remove_fns])))
  } 
  
  if (removeFE) {
    indices_to_remove <- c(indices_to_remove,which(x_rep$proc %in% c("SURF","FIRN")))
  }
  
  if (!is.null(indices_to_remove)) {
    Q_full <- Q_full[-indices_to_remove,-indices_to_remove]
    C_full <- C_full[,-indices_to_remove]
  }
  
  n_FIRN <- n_SURF <- Meshes[["SURF"]]@n + ncol(S)
  n_tot <- n_ICE + n_SURF + n_FIRN
  intbeta0 <- ((T*n_tot+1):(T*n_tot + n_ICE))
  intbeta1 <- ((T*n_tot + n_ICE + 1):(T*n_tot + 2*n_ICE))
  
  # Remove ICESAT
  indices_to_remove <- which(y_tot$Obs == "ICESAT")
  y_tot <- y_tot[-indices_to_remove,] 
  Qobs <- Qobs[-indices_to_remove,-indices_to_remove] 
  C_full <- C_full[-indices_to_remove,]
  
  # No small scale variation
  ybar = t(C_full)%*%Qobs%*%y_tot$z
  Qtot <- t(C_full)%*%Qobs%*%C_full + Q_full 
  cat("Doing Cholesky and Takahashi",sep="\n")
  X <- cholPermute(Qtot)
  Partial_Cov <- Takahashi_Davis(Qtot,cholQp = X$Qpermchol,P = X$P)
  x_margvar <- diag(Partial_Cov)
  x_mean <- as.vector(cholsolve(Qtot,ybar,perm=T,cholQp = X$Qpermchol, P = X$P))
  
  y_tot$z2 <- as.vector(y_tot$z - C_full %*% matrix(x_mean))
  
  if (!removeFE) { str1=""
  } else {
    str1="_noSMB" 
  }
  if (removeS) { str2=""
  } else {
    str2="_withS" 
  }
  
  
  
  
  
  g <- PlotAntarctica(LinePlotTheme(),shapefiles) 
  g <- g + geom_point(data=subset(y_tot,Obs=="ENVISAT" & t==3),aes(x,y,colour=pmin(pmax(z2,-0.1),0.1)),size=3) + 
    scale_colour_gradient2(low=muted("red"),mid="light yellow",high=muted("blue"),guide=guide_legend(title="m/yr")) +
    coord_fixed() 
  #coord_fixed(xlim=c(-1500,-500),ylim=c(-1000,200))
  png(filename=paste("../RATES/Publications/2013_STPaper1/ENV2006",str1,str2,".png", sep=""),width=1000,height=1000); print(g); dev.off()
  
  
  g <- PlotAntarctica(LinePlotTheme(),shapefiles) 
  g <- g + geom_point(data=subset(y_tot,Obs=="ENVISAT" & t==4),aes(x,y,colour=pmin(pmax(z2,-0.1),0.1)),size=3) + 
    scale_colour_gradient2(low=muted("red"),mid="light yellow",high=muted("blue"),guide=guide_legend(title="m/yr")) +
    coord_fixed() 
  #coord_fixed(xlim=c(-1500,-500),ylim=c(-1000,200))
  png(filename=paste("../RATES/Publications/2013_STPaper1/ENV2007",str1,str2,".png", sep=""),width=1000,height=1000); print(g); dev.off()
  
  
}

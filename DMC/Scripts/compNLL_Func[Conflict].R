#composite likelihood function:
#Function to calculate the composite likelihood
# Code Assumptions:
# 1. Grid of coordinates is a square 
# 2.
# Input:
# data : this should contain the coordinates of the grid and the values at each node
#       The grid of coordinates should be a square.
# m : vector containing (m_x, m_y) - the number of blocks in the x and y direction respectively
# Cov.fnc : Covariance matrix function. One of: Matern, ... .
# Cov.par : parameters that define the covariance matrix function - the dimension will depend on cov.fnc
# Output:
#
#
install.packages("Imap")

composite.nll <- function(data, m, Cov.fnc, Cov.par){
  n <- length(data[,1])
  delta <- 0.1
  p_begin.tilde <- min(data[,1],data[,2])
  p_end.tilde <- max(data[,1],data[,2]) + delta #need to ensure that no block falls outside of the grid. 
  p_width <- p_end.tilde - p_begin.tilde
  p_interval.x <- p_width/m[1]
  p_interval.y <- p_width/m[2]
  x_span <- seq(from = p_begin.tilde, to = p_end.tilde, by = p_interval.x)
  y_span <- seq(from = p_begin.tilde, to = p_end.tilde, by = p_interval.y)
  M <- m[1]*m[2]
  #now we need to assign the locations to the relevant block. The blocks are across x and then up y etc.
  x_span.loc <- array(0,dim=c(0))
  y_span.loc <- array(0,dim=c(0))
  for(i in 1:n){
    x_span.loc <- cbind(x_span.loc,findInterval(data[i,1],x_span)) #find interval is left-closed and right open intervals
    y_span.loc <- cbind(y_span.loc,findInterval(data[i,2],y_span))
  }
  D <- x_span.loc + m[1]*(y_span.loc - 1)    #identifying block that contains each observation.
  #create blocks that contain the indices of the data in each block:
  D_index <- append(0,which(D %in% 1))
  for(i in 2:M){
    D_index <- append(D_index,0)
    D_index <- append(D_index,which(D %in% i))
  }
  D_index <- append(D_index,0) #needs zeros at the beginning and end
  
  #Note that the indices of a new block begins when the index is zero. 
  
  D_Change <- blockChangeInd(D_index,M)
  N_arrow <- neighbours.grid(m[1],m[2])$N_arrow #this tells us which D_Change to call.
  #This code fails if any of the blocks are empty: So we create N_adj_arrow which removes all empty blocks
  block_list <- seq(1,M,1) #all possible block numbers
  remove_D <- block_list[!(block_list %in% D)] #identify blocks that do not appear in D - i.e. no observations lie within the block.
  if( length(remove_D) > 0){
    N_adj_arrow <- N_arrow[,!apply(N_arrow, 2, function(x) any(x %in% remove_D))]
  }
  else{
    N_adj_arrow <- N_arrow
  }

  #Evaluate likelihood:
  #browser()
  nll.out <- 0
  for(i in 1:dim(N_adj_arrow)[2]){ 
          nll.out <- nll.out + nll.calc(data, D_index[(D_Change[N_adj_arrow[1,i]]+1):(D_Change[N_adj_arrow[1,i]+1]-1)],
                                    D_index[(D_Change[N_adj_arrow[2,i]]+1):(D_Change[N_adj_arrow[2,i]+1]-1)],Cov.par)
  }
 return(nll.out)
}

#This version calculates the likelihood conditional on the observations within a certain distance. This is slightly different
#to the implementation of the composite likelihood above. 
## I'm not sure the conditional density observation approach is correct here. Leave this for now. 
# See Block Composite Likelihood paper and teh composite likelihood review paper by Cox et al. 
#composite.nll2 <- function( , ){
#  
#}

#######################################################################
# Auxillary functions:
# 1. i) Matern covariance matrix function.

maternCov <- function(d,par.vec){
  nu <- par.vec[1]
  rho <- par.vec[2]
  sigma <- par.vec[3]
  out <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d/rho)^nu*besselK(sqrt(2*nu)*d/rho,nu)
  #if any element d=0 then we will get NA because of division by zero. We want to replace
  #NA with sigma^2 (nugget effect)
  out[out == "NaN"] <- sigma^2 #need more robust solution to this problem.
  return(out)
}




# 1. ii) 

# 2. Calculate covariance function for 2 given blocks: Uses geodesic distance.
blockCovariance <- function(D1,D2,data,par.vec){
  #distance calculation:
  #1. sites for each block?
  S1 <- data[D1,1:2] #coordinates of sites in block 1
  S2 <- data[D2,1:2] 
  S_all <- rbind(S1,S2)
  #2.distances of sites
  dist.mat <- geoDistance(S_all,S_all)
  #3.covariance matrix:
  cov.mat <- maternCov(dist.mat,par.vec) #covariance matrix
}

# 3. Index where a new block starts:
blockChangeInd <- function(D_index, M){
  D_change <- array(0,dim=c(1,(M+1)))
  k <- 1
  for(i in 1:(M+1)){
    while(D_index[k] > 0){
      k <- k + 1
    }
    D_change[i] <- k
    k <- k + 1
  }
  return(D_change)
}

# 4. calculate the negative log-likelihood
nll.calc <- function(data, D1, D2, cov.par){
  Cov <-  blockCovariance(D1,D2,data,cov.par)
  #need to combine the covariance matrices together so that it is square ofc. 
  #then calculate the likelihood:
  data_nll <- append(data[D1,3],data[D2,3])
  out <- 1/2*log(det(Cov)) + 1/2*t(data_nll)%*%solve(Cov)%*%data_nll
}

#This code returns the pairs in N_k(arrow) from the composite likelihood paper.
neighbours.grid <- function(Nx,Ny){ #Nx and Ny are the number of blocks in the x and y directions.  
  ret<-c()
  orig.mat <- matrix(1:(Nx*Ny),Nx,Ny)
  temp.mat<-cbind(NA,rbind(NA,orig.mat,NA),NA)
  addresses <- expand.grid(x = 1:Nx, y = 1:Ny)
  for(i in 1:-1)
    for(j in 1:-1)
      if(i!=0 || j !=0)
        ret<-rbind(ret,temp.mat[addresses$x+i+1+nrow(temp.mat)*(addresses$y+j)]) 
  
  #ret is a matrix indicating which blocks are neighbours. Now we want to create N_k(arrow) 
  #We can do this by adding NA to the block pairs that are duplicates (i.e. j < i)
  N_arrow <- array(0,c(2,0))
  for(k in 1:(Nx*Ny)){
    for(j in 1:8){
      if(!is.na(ret[j,k]) && ret[j,k] > k)
        N_arrow <- cbind(N_arrow, c(k,ret[j,k]))
    }
  }
  
 list(ret=ret, N_arrow=N_arrow)
}



Convert_to_lon_lat2 <- function(X) {
  r <- sqrt(X[,1]^2 + X[,2]^2 + X[,3]^2)
  lat <- acos(X[,3]/r)
  lon <- atan2(X[,2], X[,1])
  return(cbind(lon,lat))
}

rad2deg <- function(X){ #we want to constrain angles to between certain values: can we do this without a loop?
  new_X <- X*180/pi
  n <- dim(new_X)[1]
  for(i in 1:n){
    if(new_X[i,1] < 0){
      new_X[i,1] <- 360 + new_X[i,1] 
    }
    
    new_X[i,2] <- new_X[i,2] - 90 #we want (-90,90) for the paleo code.
    
  }
  return(new_X)
}

Convert_to_xyz <- function(lonlat,r) {
  x = r[1] * cos(lonlat[,2]) * cos(lonlat[,1])
  y = r[2] * cos(lonlat[,2]) * sin(lonlat[,1])
  z = r[3] *sin(lonlat[,2])
  return(cbind(x,y,z))
}

#distGeo in the geosphere package will only compare row-wise, and return a vector. We want a matrix of distances over all combinations of rows.
geoDistance <- function(data1, data2){
  #data1 and data2 should be nx2 matrices of lon, lat coordinates:
  n_1 <- dim(data1)[1]
  n_2 <- dim(data2)[1]
  out.mat <- array(0,c(n_1,n_2))
  for(i in 1:n_1){
    for(j in 1:n_2){
      #Problem: for a small selection of coordinates gdist returns NA because the angle estimate never converges. For now,
      #we just add this small fix, we can deal with this later. is this now too slow?
      
      out.mat[i,j] <- gdist(data1[i,1],data1[i,2],data2[j,1],data2[j,2],units="km") #units = kilometers. gdist is in the Imap package.
      #The gdist function is faster than the alternative distGeo function.
      #out.mat[i,j] <- 1/1000*distGeo(data1[i,1:2],data2[j,1:2]) #distGeo is in the geosphere package. - converted to kilometers
      #default: radius of the Earth.
      if(is.na(out.mat[i,j]) ){
        out.mat[i,j] <- gdist(data1[i,1],data1[i,2] + 0.0001,data2[j,1],data2[j,2],units="km")
        #A warning will still be produced ofc.
      }
      
      
    }
  }
  return(out.mat)
}

matern.neg.lik.calc <- function(cov.param, data, n, dist.mat){
  
  Cov <- maternCov(dist.mat, cov.param)
  det.Cov <- det(Cov)
  if( det.Cov <= 0 ){
    neg.log.likelihood <- 100000000000000000
  }
  else{
    for( i in 1:n){
      neg.log.likelihood <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data  #negative log-likelihood
    }
  }
  return(neg.log.likelihood)
}

sub.matern.neg.lik.calc.constrain <- function(cov.param.2, nu, data, n, dist.mat){
  cov.param <- c(nu, cov.param.2)
  if( cov.param.2[1] > 0 && cov.param.2[2] > 0){ #rho and sigma should be strictly positive
    Cov <- maternCov(dist.mat, cov.param)
    det.Cov <- det(Cov)
    if( det.Cov <= 0 ){
      neg.log.likelihood <- 100000000000000000
    }
    else{
      for( i in 1:n){
        neg.log.likelihood <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data  #negative log-likelihood
      }
    }
  }
  else neg.log.likelihood <- 100000000000000000
  #print(list(neg.log.likelihood, cov.param, det.Cov))
  return(neg.log.likelihood)
}

sub.matern.neg.lik.calc <- function(cov.param.2, nu, data, n, dist.mat){
  cov.param <- c(nu, cov.param.2)
    Cov <- maternCov(dist.mat, cov.param)
    det.Cov <- det(Cov)
    if( det.Cov <= 0 ){
      neg.log.likelihood <- 100000000000000000
    }
    else{
      for( i in 1:n){
        neg.log.likelihood <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data  #negative log-likelihood
      }
    }
  
  #print(list(neg.log.likelihood, cov.param, det.Cov))
  return(neg.log.likelihood)
}

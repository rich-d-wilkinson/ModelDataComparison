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

composite.nll <- function(Cov.par, data, m, cov.model, m.error, scoreFunc){
  #m.error should be a vector not a matrix of error measurements.
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
 # print(N_arrow)
 # print(D_index)
  if(scoreFunc == "NegLogLik"){
    nll.out <- 0
    for(i in 1:dim(N_adj_arrow)[2]){ 
            nll.out <- nll.out + nll.calc(data, D_index[(D_Change[N_adj_arrow[1,i]]+1):(D_Change[N_adj_arrow[1,i]+1]-1)],
                                      D_index[(D_Change[N_adj_arrow[2,i]]+1):(D_Change[N_adj_arrow[2,i]+1]-1)],cov.model,Cov.par, m.error)
     }
    score.out <- nll.out
  }
  
  else if(scoreFunc == "Hyvarinen"){
    nhyv.out <- 0
    for(i in 1:dim(N_adj_arrow)[2]){ 
      nhyv.out <- nhyv.out + nhyv.calc(data, D_index[(D_Change[N_adj_arrow[1,i]]+1):(D_Change[N_adj_arrow[1,i]+1]-1)],
                                    D_index[(D_Change[N_adj_arrow[2,i]]+1):(D_Change[N_adj_arrow[2,i]+1]-1)],cov.model,Cov.par, m.error)
    }
    score.out <- nhyv.out
  }
 return(score.out)
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

maternCov <- function(d,par.vec,m.error){
  #print(d[1:3,1:3])
  #print(par.vec)
  if(length(par.vec) == 3){ #no nugget term
    nu <- par.vec[1] #smoothness (defines model sub-clas)
    rho <- par.vec[2] #range of spatial dependence
    sigma <- par.vec[3] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2])
    index_0 <- which(d > 0, arr.in = TRUE) 
    out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
  }
  else if (length(par.vec) > 3){ #nugget
    nu <- par.vec[1] #smoothness (defines model sub-clas)
    rho <- par.vec[2] #range of spatial dependence
    sigma <- par.vec[3] #partial sill
    nugget <- par.vec[4] #nugget (sill = nugget + partial sill)
    mat.error <- diag(m.error) #convert vector to diagonal matrix
    out <- matrix(sigma^2, nrow = dim(d)[1], ncol = dim(d)[2]) + nugget*mat.error
    index_0 <- which(d > 0, arr.in = TRUE) 
    out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
  }

  return(out)
}

# 1. ii) Exponential Covarianec Matrix function:

expCov <- function(d,par.vec, m.error){
  if(length(par.vec) == 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2])
    index_0 <- which(d > 0, arr.in = TRUE) 
    out[index_0] <- sigma_sSq*exp(-d[index_0]/rho)#element-wise operations.
  }
  else if (length(par.vec) > 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    nugget <- par.vec[3] #nugget (sill = nugget + partial sill)
    mat.error <- diag(m.error) #convert vector to diagonal matrix
    out <- matrix(sigma^2, nrow = dim(d)[1], ncol = dim(d)[2]) + nugget*mat.error
    index_0 <- which(d > 0, arr.in = TRUE) 
    out[index_0] <- sigma^2*exp(-d[index_0]/rho) #element-wise operations.
  }

  return(out)
}

# 1. iii) Spherical Covarianec Matrix function:

sphericalCov <- function(d,par.vec){
  nu <- par.vec[1] #smoothness (defines model sub-clas)
  rho <- par.vec[2] #range of spatial dependence
  sigma <- par.vec[3] #partial sill
  nugget <- par.vec[4:length(par.vec)] #nugget (sill = nugget + partial sill)
  sigma_sSq <- sigma^2 + nugget #measurement error. 
  out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2])
  index_0 <- which(d > 0 & d <= rho, arr.in = TRUE) #&& is not vectorised, we require & here.
  index_1 <- which(d > rho, arr.in = TRUE) 
  out[index_0] <- sigma^2*(1 - 3/2*d[index_0]/rho + 1/2*(d[index_0]/rho)^3)
  out[index_1] <- 0
  return(out)
}

testingCov <- function(d, par.vec){  #testing the covariance matrix functions in the composite likelihood.
  rho <- par.vec[1] #range of spatial dependence
  sigma <- par.vec[2] #partial sill
  out <- matrix(sigma^2, nrow = dim(d)[1], ncol = dim(d)[2])
  index_0 <- which(d > 0, arr.in = TRUE) 
  #out[index_0] <- rho #element-wise operations.
  #out[index_0] <- sigma^2*exp(-d[index_0]/rho) #element-wise operations.
  out[index_0] <- 1/d[index_0]
  return(out)
}

# 2. Calculate covariance function for 2 given blocks: Uses geodesic distance.
blockCovariance <- function(D1,D2,data,model,cov.param, m.error){
  #distance calculation:
  #1. sites for each block?
  S1 <- data[D1,1:2] #coordinates of sites in block 1
  S2 <- data[D2,1:2] 
  S_all <- rbind(S1,S2)
  #2.distances of sites
  dist.mat <- geoDistance(S_all,S_all)
  m.error.blocks <- c(m.error[D1],m.error[D2])
  #3.covariance matrix:
  #switch(model, Matern <- { Cov = maternCov(dist.mat, cov.param)}, Spherical <- {Cov = sphericalCov(dist.mat, cov.param)})
  if( model == "Matern"){
    Cov <- maternCov(dist.mat, cov.param, m.error.blocks)
  }
  else if (model == "Spherical"){
    Cov <- sphericalCov(dist.mat, cov.param, m.error.blocks)
  }
  else if (model == "Exponential"){
    Cov <- expCov(dist.mat, cov.param, m.error.blocks)
  }
  else if (model == "Test"){
    Cov <- testingCov(dist.mat, cov.param)
  }
  return(Cov)
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

# 4. a) calculate the negative log-likelihood
nll.calc <- function(data, D1, D2, cov.model, cov.par, m.error){
  Cov <-  blockCovariance(D1,D2,data,cov.model,cov.par, m.error)
  #need to combine the covariance matrices together so that it is square ofc. 
  #then calculate the likelihood:
  data_nll <- append(data[D1,3],data[D2,3])
  det.Cov <- det(Cov)
  #print("A")
  #print(cov.par)
  #print(det.Cov)
  if( det.Cov < 0.0000000000000000000001 ){
    nudge <- 0.0000000000000000000001
    adj.Cov <- Cov + nudge*diag(dim(Cov)[1])
    det.adj.Cov <- det(adj.Cov)
    #print("B")
    #print(det.adj.Cov)
    out <- 1/2*log(det.adj.Cov) + 1/2*t(data_nll)%*%solve(adj.Cov)%*%data_nll
  }
  else{
    out <-  1/2*log(det.Cov) + 1/2*t(data_nll)%*%solve(Cov)%*%data_nll
  }
  #print(out)

  return(out)
  
}

# 4. b) calculate the negative Hyvarinen Score:
nhyv.calc <- function(data, D1, D2, cov.model, cov.par, m.error){
  Cov <-  blockCovariance(D1,D2,data,cov.model,cov.par, m.error)
  #need to combine the covariance matrices together so that it is square ofc. 
  #then calculate the likelihood:
  data_nll <- append(data[D1,3],data[D2,3])
  condNo.Cov <- kappa(Cov) #don't need to calculate determinant in Hyvarinen score method.
  #print("A")
  #print(det.Cov)
  if( condNo.Cov > 1e12){
    nudge <- 0.01
    adj.Cov <- Cov + nudge*diag(dim(Cov)[1])
    det.adj.Cov <- det(adj.Cov)
    #print("B")
    #print(det.adj.Cov)
    prec.adj.Cov <- solve(adj.Cov)
    out <- - norm_vec(prec.adj.Cov%*%data_nll)^2 + 2*mat.trace(prec.adj.Cov) #negative score
  }
  else{
    prec.Cov <- solve(Cov)
    out <- - norm_vec(prec.Cov%*%data_nll)^2 + 2*mat.trace(prec.Cov) #negative score
  }
  #print(out)
  
}

norm_vec <- function(x) sqrt(sum(x^2))

mat.trace <- function(A) sum(diag(A))

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
  #     #Problem: for a small selection of coordinates gdist returns NA because the angle estimate never converges. For now,
  #     #we just add this small fix, we can deal with this later. is this now too slow?
  
        #haversine formulation:
         out.mat[i,j] <- gcd.hf(data1[i,1], data1[i,2],data2[j,1],data2[j,2]) #in km
  #     #IMAP:
  #     #out.mat[i,j] <- gdist(data1[i,1],data1[i,2],data2[j,1],data2[j,2],units="km") #units = kilometers. gdist is in the Imap package.
  #     #The gdist function is faster than the alternative distGeo function.
  #     #out.mat[i,j] <- 1/1000*distGeo(data1[i,1:2],data2[j,1:2]) #distGeo is in the geosphere package. - converted to kilometers
  #     #default: radius of the Earth.
  #     
  #     if(is.na(out.mat[i,j]) ){
  #       out.mat[i,j] <- gdist(data1[i,1] + 0.00001,data1[i,2] + 0.00001,data2[j,1],data2[j,2],units="km")
  #       #A warning will still be produced ofc.
  #     }
  #     #out.mat[j,i] <- out.mat[i,j]
  #     
    }
  }
  return(out.mat)
}


covariance.neg.lik.calc <- function(cov.param, data, m.error, dist.mat, model){
 
  #switch(model, Matern <- { Cov <- maternCov(dist.mat, cov.param)}, Spherical <- {Cov <- sphericalCov(dist.mat, cov.param)})
    if( model == "Matern"){
      Cov <- maternCov(dist.mat, cov.param, m.error)
    }
    else if (model == "Spherical"){
      Cov <- sphericalCov(dist.mat, cov.param, m.error)
    }
    else if (model == "Exponential"){
      Cov <- expCov(dist.mat, cov.param, m.error)
    }
    det.Cov <- det(Cov)
    if( det.Cov < 0.0000001 ){
      nudge <- 0.1
      adj.Cov <- Cov + nudge*diag(dim(Cov)[1])
      det.adj.Cov <- det(adj.Cov)
      #print(det.adj.Cov)
      #print("B")
      prec.adj.Cov <- solve(adj.Cov)
      neg.log.likelihood <-  1/2*log(det.adj.Cov) + 1/2*t(data)%*%prec.adj.Cov%*%data  #negative log-likelihood
    }
    else{
        neg.log.likelihood <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data  #negative log-likelihood
     
    }
    #print(list(neg.log.likelihood, cov.param, det.Cov))
  
  
  
  return(neg.log.likelihood)
}

covariance.neg.lik.calc.NM <- function(cov.param, data, m.error, dist.mat, model){
  #For Nelder-Mead optimisaion:
  
  #Exponentiate parameters to get around constraint problem:
  cov.param.pos <- cbind(exp(cov.param[1]), exp(cov.param[2]), cov.param[3])
  
  if(cov.param.pos[1] > 1/2){
    neg.log.likelihood <- 100000000
  }
  else{
    #switch(model, Matern <- { Cov <- maternCov(dist.mat, cov.param)}, Spherical <- {Cov <- sphericalCov(dist.mat, cov.param)})
    if( model == "Matern"){
      Cov <- maternCov(dist.mat, cov.param.pos, m.error)
    }
    else if (model == "Spherical"){
      Cov <- sphericalCov(dist.mat, cov.param.pos, m.error)
    }
    else if (model == "Exponential"){
      Cov <- expCov(dist.mat, cov.param.pos, m.error)
    }
    det.Cov <- det(Cov)
    if( det.Cov < 0.0000001 ){
      nudge <- 0.1
      adj.Cov <- Cov + nudge*diag(dim(Cov)[1])
      det.adj.Cov <- det(adj.Cov)
      #print(det.adj.Cov)
      #print("B")
      prec.adj.Cov <- solve(adj.Cov)
      neg.log.likelihood <-  1/2*log(det.adj.Cov) + 1/2*t(data)%*%prec.adj.Cov%*%data  #negative log-likelihood
    }
    else{
      neg.log.likelihood <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data  #negative log-likelihood
      
    }
    #print(list(neg.log.likelihood, cov.param, det.Cov))
    
  }
  
  return(neg.log.likelihood)
}

covariance.neg.lik.calc.exp.NM <- function(cov.param, data, m.error, dist.mat, model){
  #For Nelder-Mead optimisaion:
  
  #Exponentiate parameters to get around constraint problem:
  cov.param.pos <- cbind(exp(cov.param[1]), cov.param[2])
  
 
    #switch(model, Matern <- { Cov <- maternCov(dist.mat, cov.param)}, Spherical <- {Cov <- sphericalCov(dist.mat, cov.param)})
    if( model == "Matern"){
      Cov <- maternCov(dist.mat, cov.param.pos, m.error)
    }
    else if (model == "Spherical"){
      Cov <- sphericalCov(dist.mat, cov.param.pos, m.error)
    }
    else if (model == "Exponential"){
      Cov <- expCov(dist.mat, cov.param.pos, m.error)
    }
    det.Cov <- det(Cov)
    if( det.Cov < 0.0000001 ){
      nudge <- 0.1
      adj.Cov <- Cov + nudge*diag(dim(Cov)[1])
      det.adj.Cov <- det(adj.Cov)
      #print(det.adj.Cov)
      #print("B")
      prec.adj.Cov <- solve(adj.Cov)
      neg.log.likelihood <-  1/2*log(det.adj.Cov) + 1/2*t(data)%*%prec.adj.Cov%*%data  #negative log-likelihood
    }
    else{
      neg.log.likelihood <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data  #negative log-likelihood
      
    }
    #print(list(neg.log.likelihood, cov.param, det.Cov))
    

  
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

CImetric <- function(GCM, krigingVal, krigingSD){ #count of GCM output within kriging confidence interbval
  lower <- krigingVal - 2*krigingSD
  upper <- krigingVal + 2*krigingSD
  output <- 1/length(GCM)*sum(Indicator(GCM, lower, upper))
  return(output)
}

Indicator<-function(x, min=0, max=Inf)
{ 
  as.numeric(x >= min) * as.numeric(x <= max)
}



maternCov.krige <- function(data.obs, data.GCM, par.vec, m.error){
  n.1 <- dim(coordinates(data.obs))[1]
  n.2 <- dim(coordinates(data.GCM))[1]
  
  mat.error <- diag(m.error)
  if(length(par.vec) == 3){
    nu <- par.vec[1]
    rho <- par.vec[2] #range of spatial dependence
    sigma <- par.vec[3] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    d <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    if( abs(d[1,2] - d[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2]) 
      index_0 <- which(d > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
     
    }
    else{
      out <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d/rho)^nu*besselK(sqrt(2*nu)*d/rho,nu) #element-wise operations.
   
    }
    
  }
  else if (length(par.vec) > 3){
    nu <- par.vec[1]
    rho <- par.vec[2] #range of spatial dependence
    sigma <- par.vec[3] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    nugget <- par.vec[4]
    d <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    
    if( abs(d[1,2] - d[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d)[1], ncol = dim(d)[2]) + nugget*mat.error
      index_0 <- which(d > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d[index_0]/rho)^nu*besselK(sqrt(2*nu)*d[index_0]/rho,nu) #element-wise operations.
  
    }
    else{
      out <- sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d/rho)^nu*besselK(sqrt(2*nu)*d/rho,nu) #element-wise operations.
    }
    
  }
  
  return(out)
}

expCov.krige <- function(data.obs, data.GCM, par.vec, m.error){
  n.1 <- dim(coordinates(data.obs))[1]
  n.2 <- dim(coordinates(data.GCM))[1]
  
  
  if(length(par.vec) == 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    d.full <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    if( abs(d.full[1,2] - d.full[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d.full)[1], ncol = dim(d.full)[2]) 
      index_0 <- which(d.full > 0, arr.in = TRUE) 
      out[index_0] <- sigma_sSq*exp(-d.full[index_0]/rho)#element-wise operations.

    }
    else{
      out <- sigma^2*exp(-d.full/rho)#element-wise operations.

    }
    
  }
  else if (length(par.vec) > 2){
    rho <- par.vec[1] #range of spatial dependence
    sigma <- par.vec[2] #partial sill
    sigma_sSq <- sigma^2  #measurement error. 
    nugget <- par.vec[3]
    mat.error <- diag(m.error)
    d.full <- geoDistance(coordinates(data.obs), coordinates(data.GCM))
    if( abs(d.full[1,2] - d.full[2,1]) < 0.0001){
      out <- matrix(sigma_sSq, nrow = dim(d.full)[1], ncol = dim(d.full)[2]) + nugget*mat.error
      index_0 <- which(d.full > 0, arr.in = TRUE) 
      out[index_0] <- sigma^2*exp(-d.full[index_0]/rho)#element-wise operations.

    }
    else{
      out <- sigma^2*exp(-d.full/rho)#element-wise operations.

    }
    
  }
  
  return(out)
}

mCRPS.eval <- function(gcm.data, gp.data){ #the data are assumed to be the residuals so there is no mean term.
  sum <- 0
  n <- length(gcm.data)
  for(i in 1:n){
    sum <- sum - (gp.data[i,1] - gcm.data[i])^2 + 1/2*(gp.data[i,1] - gp.data[i,2])^2
  }
  out <- 1/n*sum
  return(out)
}

#We need to define a Kernel that uses geodesic distances not Euclidean distances:
geoGaussianKernel=function(x,y){
  #calculate the geodesic distance between x and y
  x <- as.matrix(t(x))
  y <- as.matrix(t(y)) #these have to be passed in as matrices for the dimension call to work.
  r <- geoDistance(x,y)
  return(exp(-(r^2))) #ignore kernel parameter for now
}

covariance.MLE <- function(initial, residual_Obs, dist.mat, matrixFunction, meas.error){
  #M is the number of different starting points to use to find the highest likelihood (problem is that rho doesn't move) 
  M <- dim(initial)[1]
  negLogLik <- 1000000000 #initialise the current best choice of parameters.
  for(m in 1:M){
    if( matrixFunction == "Matern_nug"){
      # i) Matern with nugget:
      param.lower <- c(0.05,0.01,0.01,0.00001) #(nu,rho,sigma,nugget)
      #needs upper limit for nu <= 1/2 too.
      param.upper <- c(1/2,5000,10,10) #for validity on the sphere.
      mle.output <- optim(initial[m,], covariance.neg.lik.calc, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern", lower = param.lower, 
                                     upper = param.upper, method="L-BFGS-B")
    }
    
    else if (matrixFunction == "Matern"){
      # ii) Matern without nugget:
      param.lower <- c(0.2,500,0.01) #(nu,rho,sigma)
      #needs upper limit for nu <= 1/2 too.
      param.upper <- c(1/2,10000,50) #for validity on the sphere.
      mle.output <- optim(initial[m,], covariance.neg.lik.calc, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern", lower = param.lower, 
                                  upper = param.upper, method="L-BFGS-B")
    }
    
    else if (matrixFunction == "Exponential_nug"){
      # iii) Exponential with nugget:
      param.lower <- c(0.01,0.01,0.00001) #(rho,sigma,nugget)
      param.upper <- c(10000,10,100) 
      mle.output <- optim(initial[m,], covariance.neg.lik.calc, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential", lower = param.lower, 
                               upper = param.upper, method="L-BFGS-B")
    }
    
    else if (matrixFunction == "Exponential"){
      # iv) Exponential without nugget:
      param.lower <- c(0.01,0.01) #(rho,sigma)
      param.upper <- c(10000,1000) #for validity on the sphere.
      mle.output <- optim(initial[m,], covariance.neg.lik.calc, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential", lower = param.lower, 
                                 upper = param.upper, method="L-BFGS-B")
    }
    
    else if (matrixFunction == "Spherical"){
      # mle.output.spherical <- optim(initial, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Spherical", lower = param.lower,
      #                              upper = param.upper, method="L-BFGS-B")
      # mle.output.spherical
    }
    #print(mle.output$value)
    if(mle.output$value < negLogLik){
      mle.output.best <- mle.output
      negLogLik <- mle.output$value
    }
  }
  
  return(mle.output.best)
  
}

covariance.MLE.NelderMead <- function(initial, residual_Obs, dist.mat, matrixFunction, meas.error){
  #Nedler-Mead version:
  
  #M is the number of different starting points to use to find the highest likelihood (problem is that rho doesn't move) 
  M <- dim(initial)[1]
  #To use Nelder-Mead we need to make sure that we are only passsing positive parameters. 

  negLogLik <- 1000000000 #initialise the current best choice of parameters.
  for(m in 1:M){
    if( matrixFunction == "Matern_nug"){
      # i) Matern with nugget:

      mle.output <- optim(initial[m,], covariance.neg.lik.calc.NM, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern", method="Nelder-Mead")
    }
    
    else if (matrixFunction == "Matern"){
      # ii) Matern without nugget:

      mle.output <- optim(initial[m,], covariance.neg.lik.calc.NM, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern", method="Nelder-Mead") # ,control=list(trace=TRUE))
    }
    
    else if (matrixFunction == "Exponential_nug"){
      # iii) Exponential with nugget:

      mle.output <- optim(initial[m,], covariance.neg.lik.calc.exp.NM, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential", method="Nelder-Mead")
    }
    
    else if (matrixFunction == "Exponential"){
      # iv) Exponential without nugget:

      mle.output <- optim(initial[m,], covariance.neg.lik.calc.exp.NM, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential", method="Nelder-Mead") #,control=list(trace=TRUE))
    }
    
    else if (matrixFunction == "Spherical"){
      # mle.output.spherical <- optim(initial, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Spherical", lower = param.lower,
      #                              upper = param.upper, method="L-BFGS-B")
      # mle.output.spherical
    }

    if(mle.output$value < negLogLik){
      mle.output.best <- mle.output
      negLogLik <- mle.output$value
    }
  }
  
  return(mle.output.best)
  
}


covariance.MLE.SA <- function(initial, residual_Obs, lower.b, upper.b, dist.mat, matrixFunction, meas.error){
  #Simulated Annealing version:
  
  #M is the number of different starting points to use to find the highest likelihood (problem is that rho doesn't move) 
  M <- dim(initial)[1]
  #To use Nelder-Mead we need to make sure that we are only passsing positive parameters. 
  
  negLogLik <- 1000000000 #initialise the current best choice of parameters.
  for(m in 1:M){
    if( matrixFunction == "Matern_nug"){
      # i) Matern with nugget:
      
      mle.output <- GenSA(initial[m,], covariance.neg.lik.calc, lower = lower.b, upper = upper.b, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern")
    }
    
    else if (matrixFunction == "Matern"){
      # ii) Matern without nugget:
      
      mle.output <- GenSA(initial[m,], covariance.neg.lik.calc, lower = lower.b, upper = upper.b , data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern") # ,control=list(trace=TRUE))
    }
    
    else if (matrixFunction == "Exponential_nug"){
      # iii) Exponential with nugget:
      
      mle.output <- GenSA(initial[m,], covariance.neg.lik.calc, lower = lower.b, upper = upper.b , data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential")
    }
    
    else if (matrixFunction == "Exponential"){
      # iv) Exponential without nugget:
      
      mle.output <- GenSA(initial[m,], covariance.neg.lik.calc, lower = lower.b, upper = upper.b , data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential") #,control=list(trace=TRUE))
    }
    
    else if (matrixFunction == "Spherical"){
      # mle.output.spherical <- optim(initial, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Spherical", lower = param.lower,
      #                              upper = param.upper, method="L-BFGS-B")
      # mle.output.spherical
    }
    
    if(mle.output$value < negLogLik){
      mle.output.best <- mle.output
      negLogLik <- mle.output$value
    }
  }
  
  return(mle.output.best)
  
}

covariance.MLE.PSO <- function(initial, residual_Obs, lower.b, upper.b, dist.mat, matrixFunction, meas.error){
  #Particle Swarm Optimisation approach:
  

  negLogLik <- 1000000000 #initialise the current best choice of parameters.

    if( matrixFunction == "Matern_nug"){
      # i) Matern with nugget:
      
      mle.output <- psoptim(initial, covariance.neg.lik.calc, lower = lower.b, upper = upper.b, data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern")
    }
    
    else if (matrixFunction == "Matern"){
      # ii) Matern without nugget:
      
      mle.output <- psoptim(initial, covariance.neg.lik.calc, lower = lower.b, upper = upper.b , data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Matern") # ,control=list(trace=TRUE))
    }
    
    else if (matrixFunction == "Exponential_nug"){
      # iii) Exponential with nugget:
      
      mle.output <- psoptim(initial, covariance.neg.lik.calc, lower = lower.b, upper = upper.b , data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential")
    }
    
    else if (matrixFunction == "Exponential"){
      # iv) Exponential without nugget:
      
      mle.output <- psoptim(initial, covariance.neg.lik.calc, lower = lower.b, upper = upper.b , data = residual_Obs, m.error =  meas.error, dist.mat = dist.mat, model = "Exponential") #,control=list(trace=TRUE))
    }
    
    else if (matrixFunction == "Spherical"){
      # mle.output.spherical <- optim(initial, covariance.neg.lik.calc, data = residual_Obs, n = n, m.error =  meas.error, dist.mat = dist.mat, model = "Spherical", lower = param.lower,
      #                              upper = param.upper, method="L-BFGS-B")
      # mle.output.spherical
    }
    

  
  return(mle.output)
  
}

deg2rad <- function(deg) return(deg*pi/180) #convert degrees to radians

#haversine formula for great-circle distance: Takes degrees and then converts to radians.
gcd.hf <- function(long1, lat1, long2, lat2) {
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

##
#Covariance matrix calculation by Rich.

create_K <- function(location1, location2){
  
  klist <- lapply(1:2,FUN=create_k)
  out <- abs(outer(location1,location2, FUN="-"))
  tmp<-lapply(1:2, FUN=function(i) klist[[i]](out[,i,,i])     )
  K = tmp[[1]]*tmp[[2]]
  return(K)
}

krige.sim <- function(Obs, dataframe, new.locations, model, cov.model, par.vecs, meas.error, M){ #draw realisations from a fitted Gaussian Process

  colnames(new.locations) <- c("lon","lat")
  location.points <- data.frame(new.locations)
  detach(location.points)
  attach(location.points)
  coordinates(location.points) = ~ lon + lat
  proj4string(location.points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84") #this is a lon,lat grid projection of the Earth (I believe.)
  
  GP.mean <- krige0(Obs ~ 1, dataframe, location.points, model, par.vec = par.vecs, m.error = meas.error)
  dist.mat <- geoDistance(new.locations, new.locations)
  new.m.error <- array(1, dim=c(dim(new.locations)[1],1))
  if( cov.model == "Matern"){
    GP.Cov <- maternCov(dist.mat, par.vec, new.m.error)
  }
  else if (cov.model == "Spherical"){
    GP.Cov <- sphericalCov(dist.mat, par.vec, new.m.error)
  }
  else if (cov.model == "Exponential"){
    GP.Cov <- expCov(dist.mat, par.vec, new.m.error)
  }
  output <- mvrnorm(M, GP.mean, GP.Cov)
  t.output <- t(output)
  return(t.output)

}


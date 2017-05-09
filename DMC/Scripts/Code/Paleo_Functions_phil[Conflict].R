Prec_from_lattice <- function(Grid,ds)
{
  n <- nrow(Grid)
  #### set up distance and neighbourhood (W, based on sharing a common border) matrices
  distance <-array(0, c(n,n))
  W <- as.matrix(dist(Grid))==ds
  n.neighbours <- as.numeric(apply(W, 1, sum))
  Q <- diag(n.neighbours) - W
  return(as(Q,"dgCMatrix"))
}

# interpolates a continuous data set onto a grid
nn_grid_interp <- function(s,df,delta=10,miss_value=NA) {
  if (length(delta) == 1) {
    s_rnd <- data.frame(round(s/delta)*delta)
  } else {
    s_rnd <- data.frame(round(s[,1]/delta[1])*delta[1],round(s[,2]/delta[2])*delta[2])
  }
  s_rnd[,3] <- 1:nrow(s_rnd)
  names(s_rnd) = c("x","y","n")    
  df_sub <- subset(df,(x >= min(s_rnd$x) &
                         x <= max(s_rnd$x) &
                         y >= min(s_rnd$y) &
                         y <= max(s_rnd$y))) 
  df <- merge(df_sub,s_rnd,all.y=T)               
  df <- arrange(df,n)
  if (any(is.na(df$z)))  df[is.na(df$z),]$z <- miss_value
  return(df$z)
}

cholPermute <- function(Q)  {
  n <- dim(Q)[1]
  
  e <-tryCatch({ symchol <- Cholesky(Q)},error= function(temp) {print("Cholesky failed, coercing to symmetric")},finally="Cholesky successful")
  if (class(e) == "character")  {
    symchol <- Cholesky(forceSymmetric(Q))
  }
  
  
  j <- 1:n
  i <- symchol@perm + 1
  P <- sparseMatrix(i,j,x=rep(1,n))
  if (class(e) == "character")  {
    Qpermchol <- t(chol(forceSymmetric(t(P)%*%Q%*%P)))
  } else { Qpermchol <- t(chol(t(P)%*%Q%*%P)) }
  return(list(Qpermchol=Qpermchol,P=P))
}




## Solve Qx = y
cholsolve <- function(Q,y,perm=F,cholQ = matrix(1,0,0),cholQp = matrix(1,0,0),P=NA)  {
  if (perm == F) {
    if (dim(cholQ)[1] == 0) {
      e <-tryCatch({L <- t(chol(Q))},error= function(temp) {print("Cholesky failed, coercing to symmetric")},finally="Cholesky successful")
      if (class(e) == "character") {
        L <- t(chol(forceSymmetric(Q))) }
    }  else {
      L <- cholQ
    }
    
    v <- solve(L,y)
    x <- solve(t(L),v)
  }
  if (perm == T) {
    if (dim(cholQp)[1] == 0) {
      QP <- cholPermute(Q)
      Lp <- QP$Qpermchol
      P <- QP$P
    } else {
      Lp <- cholQp
    }
    
    v <- solve(Lp,t(P)%*%y)
    w <- solve(t(Lp),v)
    x <- P%*%w
  }
  return(x)
}





Takahashi_Davis <- function(Q,return_perm_chol = 0,cholQp = matrix(0,0,0),P=0) {
  
  n <- dim(Q)[1]
  
  if (dim(cholQp)[1] == 0) {
    symchol <- Cholesky(forceSymmetric(Q))
    j <- 1:n
    i <- symchol@perm + 1
    P <- sparseMatrix(i,j,x=rep(1,n))
    Lperm <- L <- t(chol(t(P)%*%Q%*%P))
  } else {
    Lperm <- L <- cholQp
    P <- P
  }
  
  d = diag (L)
  Ld <- L%*%sparseMatrix(i=1:n,j=1:n,x=1/d)
  L <- tril(Ld,-1)
  U = t(L)
  d = d^2
  D = sparseMatrix(i=1:n,j=1:n,x=d)
  
  R = as(L+D,"dgTMatrix")
  R = sparseMatrix(i=R@i+1,j=R@j+1,x=1) ;
  Zpattern = R + t(R) - sparseMatrix(i=1:n,j=1:n,x=1)
  
  Z <- sparseinv_wrapper(L,d,t(U),Zpattern)
  if (return_perm_chol == 0) {
    return(P%*%Z%*%t(P))
  } else {
    return(list(S=P%*%Z%*%t(P),Lp = Lperm,P=P))
  }
  
}

sparseinv_wrapper <- function(L,d,U,Zpattern) {
  
  n <- dim(L)[1]
  Lp <- L@p
  Li <- L@i
  Lx <- L@x
  
  Up <- U@p
  Uj <- U@i
  Ux <- U@x
  
  Zpatp <- Zpattern@p
  Zpati <- Zpattern@i
  znz = Zpatp [n+1]
  
  if (.Platform$OS.type == "windows") {
    dyn.load("./Code/sparseinv/R_WINDOWS/sparseinvR.dll")
  } else {
    dyn.load("./Code/sparseinv/R_UNIX/sparseinvR.so")
  }
  X <- .C("sparseinv",as.integer(n),as.integer(Lp),as.integer(Li),as.double(Lx),as.double(d),as.integer(Up),as.integer(Uj),as.double(Ux),as.integer(Zpatp),as.integer(Zpati),result = double(znz))
  
  Lower <- L + sparseMatrix(i=1:n,j=1:n,x = d)
  Zp <- Lower@p
  Zi <- Lower@i
  Zx <- X$result[-which(X$result == 0)]
  Z <- sparseMatrix(p = Zpatp, i =Zpati, x = X$result,index1=F)
  return(Z)
}

logLik_prop <- function(x,mu,Q) {
  return( 0.5*logdet(chol(Q)) - as.vector(0.5*t(x - mu) %*% Q %*% (x - mu) ))
}

CARfit <- function(x,y,yobs) {
  ds1 <- mean(diff(x))
  ds2 <- mean(diff(y))
  if(ds1 !=  ds2) {
    y  = y*ds1/ds2
    yobs$y  = yobs$y*ds1/ds2
  }
  ds <- ds1
  Grid <- as.data.frame(expand.grid(x, y))
  Q <- Prec_from_lattice(Grid,ds)
  n <- nrow(Grid)
  names(Grid) <- c("x","y")
  Grid$z <- 1:nrow(Grid)
  m <- nrow(yobs)
  Qobs <- sparseMatrix(i=1:m,j=1:m,x=1/(yobs$std^2))
  Ci <- 1:nrow(yobs)
  Cj <- nn_grid_interp(yobs,Grid,delta=ds)
  C <- sparseMatrix(i = Ci, j = Cj, x=1,,dims=c(m,n))  
  
  ybar = t(C)%*%Qobs%*%(yobs$z) 
  Qtot <- t(C)%*%Qobs%*%C + Q 
  X <- cholPermute(Qtot)
  Partial_Cov <- Takahashi_Davis(Qtot,cholQp = X$Qpermchol,P = X$P)
  Grid$std <- as.vector(sqrt(diag(Partial_Cov)))
  Grid$mean <- as.vector(cholsolve(Qtot,ybar,perm=T,cholQp = X$Qpermchol, P = X$P))
  Grid$y <- Grid$y * ds2/ds1
  return(list(data=Grid,Q = Qtot))
}

Infer <- function(Inf_mats) {
  C_full = Inf_mats$C_full
  y_obs = Inf_mats$y_obs
  Qobs = Inf_mats$Qobs
  Q = Inf_mats$Q
  
  ybar = t(C_full)%*%Qobs%*%matrix(y_obs$val) 
  Qtot <- t(C_full)%*%Qobs%*%C_full + Q 
  X <- cholPermute(Qtot)
  Partial_Cov <- Takahashi_Davis(Qtot,cholQp = X$Qpermchol,P = X$P)
  Resstd <- as.vector(sqrt(diag(Partial_Cov)))
  Resmean <- as.vector(cholsolve(Qtot,ybar,perm=T,cholQp = X$Qpermchol, P = X$P))
  Res = list(mean = Resmean, std=Resstd)
  Inf_mats$Qtot = Qtot
  Inf_mats$Res = Res
  return(c(Inf_mats))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

marg_prec_from_kappa <- function(kappa_l,nu) {
  return(1/(gamma(nu)/(gamma(nu+1)*4*pi*kappa_l^(2*nu))))
}


Prec_from_SPDE <- function(M,K,tau,kappa,alpha=1)  {
  if (!(alpha %in% c(1,2,3,4))) {
    stop("alpha > 4 not implemented yet")
  }
  n <- nrow(M)
  if(class(tau) == "numeric") {
    tau <- sparseMatrix(i=1:n,j=1:n,x=tau)
  }
  if(class(kappa) == "numeric") {
    kappa <- sparseMatrix(i=1:n,j=1:n,x=kappa)
  }
  
  M_approx <- sparseMatrix(i=1:n,j=1:n,x=rowSums(M))
  M_approx_inv <- sparseMatrix(i=1:n,j=1:n,x=1/rowSums(M))
  M_kappa2 <- kappa%*%M%*%kappa
  G <- (M_kappa2 + K)
  if (alpha == 1) {
    Q <- tau%*%G%*%tau
  } else if (alpha ==2) {
    Q <- tau%*%G%*%M_approx_inv%*%G%*%tau
  } else if (alpha == 3) {
    Q <-  tau%*%G%*%M_approx_inv%*%G%*%M_approx_inv%*%G%*%tau
  } else if (alpha == 4) {
    Q <-  tau%*%G%*%M_approx_inv%*%G%*%M_approx_inv%*%G%*%M_approx_inv%*%G%*%tau
  }
  return(Q)
}

Prec_from_SPDE_wrapper <- function(M,K,intrinsic,nu,desired_prec,l) {
  kappa_l = sqrt(8*nu)/l
  marg_prec <- marg_prec_from_kappa(kappa_l,nu)
  tau <- sqrt(desired_prec/marg_prec)
  Q <- Prec_from_SPDE(M,K,tau=tau,kappa=kappa_l,alpha=nu+1)
}


GMRF_init <- function(mu = 0, Q = 1, intrinsic = 0,perm=T,diag_offset=0.00000000001, Chol=0) {
  n = length(mu)
  if (intrinsic == 0)   {
    Quncon <- Q
    this_class <- 'GMRF'
  }
  # Create a new IGMRF of first order
  if (intrinsic == 1)    {           # See Rue & Held pg. 108
    #Quncon <- (Q + 1*matrix(1,n,n))
    Quncon <- Q + diag_offset*sparseMatrix((1:n),(1:n),x=1)
    this_class <- 'IGMRF1'
  }
  
  if (intrinsic == 2) {
    A = matrix(1,2,n)
    A[2,] <- 1:n
    #Quncon <- (Q + t(A)%*%A)
    Quncon <- Q + diag_offset*sparseMatrix((1:n),(1:n),x=1)
    this_class <- 'IGMRF2'
  }
  
  if(!is.list(Chol)) {
    cholQ <- t(chol(Quncon))
    
    if (perm == T) {
      # Find permutation matrices
      # The default used is AMD
      # See http://www.cise.ufl.edu/research/sparse/cholmod/CHOLMOD/Doc/UserGuide.pdf
      QP <- cholPermute(Quncon)
      obj <- list(mu = mu,Q = Q,cholQ = cholQ,cholQ.perm = QP$Qpermchol,P = QP$P,Quncon = Quncon)
      # Here fill-in is defined as the ratio of non-zeros in L to the non-zeros in the lower triangular part of Q
      #cat(paste("Fill-in before permutation=" ,length(which(t(cholQ)!=0))/length(which(tril(Q)!=0))),sep="\n")
      #cat(paste("Fill-in after permutation (AMD)=" ,length(which(t(QP$Qpermchol)!=0))/length(which(tril(Q)!=0))),sep="\n")
    }
    
    else {
      obj <- list(mu = mu,Q = Q,cholQ = cholQ,Quncon = Quncon)
    }
  } else {
    obj <- list(mu = mu,Q = Q,cholQ.perm = Chol$Qpermchol,P = Chol$P)
  }
  class(obj) <- this_class
  return(obj)
}


## Sample from the GMRF (both intrinsic and non-intrinsic)
SampleGMRF <- function(G,reps=1,use_perm=F) {
  n = length(G$mu)
  x <- matrix(0,n,reps)
  z <- matrix(rnorm(n*reps),n,reps)
  
  # Algorithm 2.4, Rue and Held
  if (use_perm ==F) {
    for (i in (1:reps)) {
      v <- solve(t(G$cholQ),z[,i])
      x[,i] <- matrix(G$mu + v)
    }
  }
  
  if (use_perm == T) {
    for (i in (1:reps)) {
      #zbar <- t(G$P)%*%G$cholQ%*%z[,i]
      #v1 <- solve(G$cholQ.perm,zbar)
      #v2 <- solve(t(G$cholQ.perm),v1)
      #x[,i] <- matrix(G$mu + G$P%*%v2)
      
      v <- G$P %*% solve(t(G$cholQ.perm), z[,i])
      x[,i] <- matrix(G$mu + v)
      
    }
  }
  
  if(class(G) == "IGMRF1")  {    # Algorithm 2.6, Rue and Held - correct
    A = matrix(1,1,n)
    P1 <- solve(G$cholQ,t(A))
    V <- solve(t(G$cholQ),P1)
    W <- A%*%V
    cholW = chol(W)
    P2 <- solve(cholW,t(V))
    U  <- solve(t(cholW),P2)
    for (i in (1:reps)) {
      cvec <- A %*% x[,i]                                                     
      x[,i] <- x[,i] - matrix(t(U) %*% cvec)
    }
  }
  
  if(class(G) == "IGMRF2")  {    # Algorithm 2.6, Rue and Held - correct
    A = matrix(1,2,n)
    A[2,] <- 1:n
    P1 <- solve(G$cholQ,t(A))
    V <- solve(t(G$cholQ),P1)
    W <- A%*%V
    U  <- solve(W,t(V))
    for (i in (1:reps)) {
      cvec <- A %*% x[,i]
      x[,i] <- x[,i] - matrix(t(U) %*% cvec)
    }
  }
  
  
  return(x)
  
}


Convert_to_lon_lat <- function(X,r) {
  lat <- asin(X[,3]/r)
  lon <- atan2(X[,2], X[,1])
  return(cbind(lon,lat))
}

Convert_to_lon_lat2 <- function(X) {
  #We update the code so that it can take a vector:
  n <- dim(X)[1]
  r <- sqrt(X[,1]^2 + X[,2]^2 + X[,3]^2)
  lat <- acos(X[,3]/r)
  lon <- atan2(X[,2], X[,1])
  return(cbind(lon,lat))
}

Convert_to_xyz <- function(lonlat,r) {
  x = r[1] * cos(lonlat[,2]) * cos(lonlat[,1])
  y = r[2] * cos(lonlat[,2]) * sin(lonlat[,1])
  z = r[3] *sin(lonlat[,2])
  return(cbind(x,y,z))
}


## Find C matrix when observations are isolated points
FindC_sphere <- function(p,tri,locs,r) {
  if (length(locs[[1]]) > 0)  {
    # See http://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates     
    
    protate <- function(p,angle) {
      p <- p + angle
      p[p < -pi] = p[p < -pi] + 2*pi
      p[p > pi] = p[p  > pi] - 2*pi
      return(p)
    }
    
    
    p_lon_lat <- Convert_to_lon_lat(p,r[3])
    locs_xyz <- Convert_to_xyz(cbind(locs[[1]],locs[[2]]),r)
    cat("Observations close to East-West boundary will cause error",sep="\n")
    
    
    
    
    # OK so the way to sort this is to do the below both like this and for z-axis interchanged
    # with x-axis or y-axis. This will work if we only do the latitudes < 80 for both axes.
    
    # We need to do some geomatric acrobatics to find which triangle the output is at lon boundaries
    remove_tri <- NULL
    count =0
    for (i in 1:nrow(tri)) {
      this_p <- p_lon_lat[tri[i,],]
      if (diff(range(this_p[,1])) > r[1]) remove_tri <- c(remove_tri,i)

    }
    tri_temp <- tri[-remove_tri,]
    t_num <- tsearch(p_lon_lat[,1], p_lon_lat[,2], tri_temp, locs[[1]], locs[[2]], bary = FALSE) # Find triangles in lat/lon
    for (i in 1:length(t_num)) {
       t_num[i] <- (1:nrow(tri))[-remove_tri][t_num[i]] # Correct indices for what we removed
    }
    
    # Now repeat the same exercise with everything rotated by 90 degrees
    remove_tri <- NULL
    for (i in 1:nrow(tri)) {
      this_p <- p_lon_lat[tri[i,],]
      this_p[,1] <- protate(this_p[,1],pi/2)
      if (diff(range(this_p[,1])) > r[1]) remove_tri <- c(remove_tri,i)
    }
    tri_temp2 <- tri[-remove_tri,]
    t_num2 <- tsearch(protate(p_lon_lat[,1],pi/2), p_lon_lat[,2], tri_temp2, protate(locs[[1]],pi/2), locs[[2]], bary = FALSE) # Find triangles in lat/lon
    for (i in 1:length(t_num2)) {
      t_num2[i] <- (1:nrow(tri))[-remove_tri][t_num2[i]] # Correct indices for what we removed
    }
    t_num[is.na(t_num)] = t_num2[is.na(t_num)]

#     # Version 2
#     for (i in 1:length(locs[[1]])) {
#       this_obs <- c(locs[[1]][i],locs[[2]][i])
#       this_p <- cbind(p_lon_lat[,1] - this_obs[1],p_lon_lat[,2] - this_obs[2])
#       this_p[this_p[,1] < -pi,1] = this_p[this_p[,1] < -pi,1] + 2*pi
#       this_p[this_p[,1] > pi,1] = this_p[this_p[,1]  > pi,1] - 2*pi
#       this_p[this_p[,2] < -pi/2,2] = this_p[this_p[,2] < -pi/2,2] + pi
#       this_p[this_p[,2] > pi/2,2] = this_p[this_p[,2]  > pi/2,2] - pi
#       this_p_xyz <- Convert_to_xyz(this_p,r)
#       this_p2 <- Convert_to_lon_lat(this_p_xyz,r[3])
#     }
#     
#     
#     # Version 2
#     t_num_new<- rep(0,length(locs[[1]]))
#     for (i in 1:length(locs[[1]])) {
#       this_obs <- c(locs[[1]][i],locs[[2]][i])
#       #Perform an albers projection centred at this coordinate
#       p_proj <- mapproject(p_lon_lat[,1],p_lon_lat[,2],projection="trapezoidal",par=c(locs[[1]][i],locs[[2]][i]))
#       obs_proj <- mapproject(locs[[1]][i],locs[[2]][i],projection="trapezoidal",par=c(locs[[1]][i],locs[[2]][i]))
#       p_proj <- mapproject(p_lon_lat[,1],p_lon_lat[,2],projection="rectangular",par=c(locs[[2]][i]))
#       obs_proj <- mapproject(locs[[1]][i],locs[[2]][i],projection="rectangular",par=c(locs[[2]][i]))
#       t_num_new[i] <- tsearch(p_proj$x,p_proj$y, tri, obs_proj$x,obs_proj$y, bary = FALSE) # Find triangles in lat/lon
#     }
#     tri_df <- as.data.frame(tri)
#     t_num <- rep(0,length(locs[[1]]))
#     # Use 3D information
#     for (i in 1:length(locs[[1]])) {
#       dist <- data.frame(dist = sqrt((p[,1] - locs_xyz[i,1])^2 + (p[,2] - locs_xyz[i,2])^2 + (p[,3] - locs_xyz[i,3])^2),
#                          n = 1:nrow(p))
#       c = sort(dist)[1,2]
#       tris <- subset(tri_df, (V1 == c[1] | V2 == c[1] | V3 == c[1])) # These are the possible triangles holding the points. Now search
#       in_tri <- rep(0,6)
#       for(j in 1:nrow(tris)) {
#         
#         # Project triangle onto a plane
#         this_p <- p[as.numeric(tris[j,]),]
#         d1 = sqrt((this_p[1,1] - this_p[1,2])^2 + (this_p[2,1] - this_p[2,2])^2  + (this_p[3,1] - this_p[3,2])^2)
#         d2 = sqrt((this_p[1,1] - this_p[1,3])^2 + (this_p[2,1] - this_p[2,3])^2  + (this_p[3,1] - this_p[3,3])^2)
#         d3 = sqrt((this_p[1,2] - this_p[1,3])^2 + (this_p[2,2] - this_p[2,3])^2  + (this_p[3,2] - this_p[3,3])^2)
#         p_proj <- matrix(0,3,2)
#         p_proj[2,] = c(d1,0);
#         xx = (d2^2 - d3^2 + d1^2)/(2*d1);
#         p_proj[3,] = c(xx,sqrt(d2^2 - xx^2));
#         
#         # Project obs onto the plane
#         this_obs <- Convert_to_xyz(cbind(locs[[1]][i],locs[[2]][i]),r)
#         T <- rbind(t(p_proj),c(0,0,0)) %*% solve(t(this_p)) # transofmration matrix
#         obs_proj <- (T %*% t(this_obs))[1:2,]
#         
#         in_tri[j] <- pnt.in.poly(t(matrix(obs_proj)),p_proj)$pip
#       }
#       t_num[i] <- as.numeric(row.names(tris))[which(in_tri == 1)[1]]
#     }
#     
    
  # GOOD BUT STILL HAVE MULTIPLE NAs at the top edges for some reason
    
   z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
    b <- matrix(0,3,1)
    A <- matrix(0,3,3)
    A[,1] = 1
    
    for (i in 1:length(t_num)) {
      t_num_i <- t_num[i]
      this_tri <- tri[t_num_i,]
      this_p <- p[this_tri,]
      
      # Correct for longitude boundary
      if (diff(range(this_p[,1])) > r[1]) {
        this_p[,1] <- protate(this_p[,1],pi/2)
      }
      
      # Project triangle onto a plane
      d1 = sqrt((this_p[1,1] - this_p[1,2])^2 + (this_p[2,1] - this_p[2,2])^2  + (this_p[3,1] - this_p[3,2])^2)
      d2 = sqrt((this_p[1,1] - this_p[1,3])^2 + (this_p[2,1] - this_p[2,3])^2  + (this_p[3,1] - this_p[3,3])^2)
      d3 = sqrt((this_p[1,2] - this_p[1,3])^2 + (this_p[2,2] - this_p[2,3])^2  + (this_p[3,2] - this_p[3,3])^2)
      p_proj <- matrix(0,3,2)
      p_proj[2,] = c(d1,0);
      xx = (d2^2 - d3^2 + d1^2)/(2*d1);
      p_proj[3,] = c(xx,sqrt(d2^2 - xx^2));
      
      A[,2:3] <- p_proj
      Ainv <- solve(A) 
      
      
      # Project obs onto the plane
      this_obs <- Convert_to_xyz(cbind(locs[[1]][i],locs[[2]][i]),r)
      T <- rbind(t(p_proj),c(0,0,0)) %*% solve(t(this_p)) # transofmration matrix
      obs_proj <- (T %*% t(this_obs))[1:2,]
      
      for (j in 1:3) {
        b[,] <- 0
        b[j] <- 1
        i_ind[i,j] <- i
        j_ind[i,j] <- this_tri[j]
        z[i,j] <- matrix(c(1,obs_proj[1],obs_proj[2]),1,3)%*%Ainv%*%b
      }
    }
    
    C <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                      dims = c(length(locs[[1]]),dim(p)[1]))
    return(C)
  }  else {
    return(matrix(1,0,0))
  }
  
}


## Find C matrix when observations are isolated points
FindC_sphere2 <- function(p,tri,locs,r) {
  if (length(locs[[1]]) > 0)  {
    # See http://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates     
    
    protate <- function(p,angle) {
      p <- p + angle
      p[p < -pi] = p[p < -pi] + 2*pi
      p[p > pi] = p[p  > pi] - 2*pi
      return(p)
    }
    locs_xyz <- Convert_to_xyz(cbind(locs[[1]],locs[[2]]),r)
    t_num <- matrix(length(locs[[1]]),2)
    
    t_try <- matrix(0,length(locs[[1]]),2)
    for(j in 1:2) {
      if (j == 1) {
        p_lon_lat <- Convert_to_lon_lat(p,r[3])
        locs_vec <-  Convert_to_lon_lat(locs_xyz,r[3])
      } else {
        p_lon_lat <- Convert_to_lon_lat(cbind(p[,1],p[,3],p[,2]),r[3])
        locs_vec <-  Convert_to_lon_lat(cbind(locs_xyz[,1],locs_xyz[,3],locs_xyz[,2]),r[3])
      }
      locs[[1]] <- locs_vec[,1]
      locs[[2]] <- locs_vec[,2]
      
        # OK so the way to sort this is to do the below both like this and for z-axis interchanged
        # with x-axis or y-axis. This will work if we only do the latitudes < 80 for both axes.
        
        # We need to do some geomatric acrobatics to find which triangle the output is at lon boundaries
        remove_tri <- NULL
        count =0
        for (i in 1:nrow(tri)) { # Find which triangles to remove
	      # print(paste0("FindC_sphere2: ",i),quote=FALSE)
	      # flush.console()
          this_p <- p_lon_lat[tri[i,],]
          if (diff(range(this_p[,1])) > r[1]) remove_tri <- c(remove_tri,i)
        }
        tri_temp <- tri[-remove_tri,]
        t_num <- tsearch(p_lon_lat[,1], p_lon_lat[,2], tri_temp, locs[[1]], locs[[2]], bary = FALSE) # Find triangles in lat/lon
        for (i in 1:length(t_num)) {
          t_num[i] <- (1:nrow(tri))[-remove_tri][t_num[i]] # Correct indices for what we removed
        }
        
        # Now repeat the same exercise with everything rotated by 90 degrees
        remove_tri <- NULL
        for (i in 1:nrow(tri)) {
            # print(paste0("FindC_sphere2 - 90: ",i),quote=FALSE)
            # flush.console()
         this_p <- p_lon_lat[tri[i,],]
          this_p[,1] <- protate(this_p[,1],pi/2)
          if (diff(range(this_p[,1])) > r[1]) remove_tri <- c(remove_tri,i)
        }
        tri_temp2 <- tri[-remove_tri,]
        t_num2 <- tsearch(protate(p_lon_lat[,1],pi/2), p_lon_lat[,2], tri_temp2, protate(locs[[1]],pi/2), locs[[2]], bary = FALSE) # Find triangles in lat/lon
        for (i in 1:length(t_num2)) {
          t_num2[i] <- (1:nrow(tri))[-remove_tri][t_num2[i]] # Correct indices for what we removed
        }
        t_num[is.na(t_num)] = t_num2[is.na(t_num)]
        t_try[,j] <- t_num
    
    }
    
    # There are very minor differences between the two due to the project but we don't care about points which are largely on the boundary so we take 
    # the first column and only those of the second which are NA in the first
    which_na <- is.na(t_try[,1])
    t_try[which_na,1] <- t_try[which_na,2]
    t_num <- t_try[,1]
    # GOOD BUT STILL HAVE MULTIPLE NAs at the top edges for some reason
    
    z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
    b <- matrix(0,3,1)
    A <- matrix(0,3,3)
    A[,1] = 1
    
    # Convert back to default
    p_lon_lat <- Convert_to_lon_lat(p,r[3])
    locs_vec <-  Convert_to_lon_lat(locs_xyz,r[3])
    locs[[1]] <- locs_vec[,1]
    locs[[2]] <- locs_vec[,2]
    
    for (i in 1:length(t_num)) {
      t_num_i <- t_num[i]
      this_tri <- tri[t_num_i,]
      this_p <- p[this_tri,]
      
      # Correct for longitude boundary
         print(paste0("FindC_sphere - lon bound: ",i),quote=FALSE)
         flush.console()
      if (diff(range(this_p[,1])) > r[1]) {
        this_p[,1] <- protate(this_p[,1],pi/2)
      }
      
      # Project triangle onto a plane
      d1 = sqrt((this_p[1,1] - this_p[1,2])^2 + (this_p[2,1] - this_p[2,2])^2  + (this_p[3,1] - this_p[3,2])^2)
      d2 = sqrt((this_p[1,1] - this_p[1,3])^2 + (this_p[2,1] - this_p[2,3])^2  + (this_p[3,1] - this_p[3,3])^2)
      d3 = sqrt((this_p[1,2] - this_p[1,3])^2 + (this_p[2,2] - this_p[2,3])^2  + (this_p[3,2] - this_p[3,3])^2)
      p_proj <- matrix(0,3,2)
      p_proj[2,] = c(d1,0);
      xx = (d2^2 - d3^2 + d1^2)/(2*d1);
      p_proj[3,] = c(xx,sqrt(d2^2 - xx^2));
      
      A[,2:3] <- p_proj
      Ainv <- solve(A) 
      
      
      # Project obs onto the plane
      this_obs <- Convert_to_xyz(cbind(locs[[1]][i],locs[[2]][i]),r)
      T <- rbind(t(p_proj),c(0,0,0)) %*% solve(t(this_p)) # transformation matrix
      obs_proj <- (T %*% t(this_obs))[1:2,]
      
      for (j in 1:3) {
        b[,] <- 0
        b[j] <- 1
        i_ind[i,j] <- i
        j_ind[i,j] <- this_tri[j]
        z[i,j] <- matrix(c(1,obs_proj[1],obs_proj[2]),1,3)%*%Ainv%*%b
      }
    }
    
    C <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                      dims = c(length(locs[[1]]),dim(p)[1]))
    return(C)
  }  else {
    return(matrix(1,0,0))
  }
  
}


OverlayPlot <- function(GG,GGstd=NULL,leg_title="",do_mesh=1,zlo = -0.6,zhi = 0.6,alphalo = 0.2, alphahi = 1) {
  GG$zlo = zlo
  GG$zhi = zhi
  if ((sign(zhi) - sign(zlo)) == 2) {
    g <- LinePlotTheme() + geom_tile(data=GG,aes(x=x,y=y,fill=pmin(pmax(z,zlo),zhi))) +
      scale_fill_gradient2(low=(muted("blue")),mid="light yellow",high=(muted("red")),
                           guide=guide_colourbar(title=leg_title),limits=c(zlo,zhi)) 
  } else {
    g <- LinePlotTheme() + geom_tile(data=GG,aes(x=x,y=y,fill=pmin(pmax(z,zlo),zhi))) +
      scale_fill_gradient(low="light yellow",high=(muted("green")),
                           guide=guide_colourbar(title=leg_title),limits=c(zlo,zhi)) 
  }
    
  
  if (!is.null(GGstd)){
    GGstd$alphalo = alphalo
    GGstd$alphahi = alphahi
    g <- g+ geom_tile(data=GGstd,aes(x=x,y=y, alpha=-pmax(pmin(z,alphahi),alphalo)),fill="white")
    
  }
  #   g <- g + theme(panel.background = element_rect(fill='white', colour='white'),panel.grid=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
  #          axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank()) +
  g <- g + coord_fixed(ratio=1) + xlab("") + ylab("")
  g <- g + theme(legend.position="right") + scale_alpha(guide = 'none')
  
  # theme(legend.key.width=unit(4,"lines")) + theme(legend.key.height=unit(4,"lines")) +
  # theme(text = element_text(size=40))
  return(g)
}

Overlay_data <- function(g,Obs) {
  g <- g + geom_point(data=Obs,aes(x,y),size=5,colour="black") + 
    geom_point(data=Obs,aes(x,y,colour=pmin(pmax(z,-10),10)),size=4) + 
    scale_colour_gradient2(low=muted("blue"),mid="light yellow",high=muted("red"),limits=c(-10,10))
}
Add_coast <- function(g,coastline) {
  g <- g + geom_polygon(data=coastline,aes(x,y,group=id))
}

Triang_to_grid <- function(df,name,detail,wrap=T) {
  
  if(wrap){
   leftoverlap = subset(df, x < (-pi + pi))
   leftoverlap$x = leftoverlap$x + 2*pi
   rightoverlap = subset(df, x > pi - pi)
   rightoverlap$x = rightoverlap$x - 2*pi
   
   df <- rbind(df,leftoverlap,rightoverlap)
   
   detail = detail*2
  }
  
  G <- interp(df$x,df$y,df[name][,1],xo=seq(min(df$x),max(df$x),length=detail),yo=seq(min(df$y),max(df$y),length=detail))
  GG <- cbind(expand.grid(G$x,G$y),as.vector(G$z))
  names(GG) <- c("x","y","z")
  if(wrap) {
   GG <- subset(GG,x > -pi & x < pi) 
  }
  return(GG)
}

LandAttribute <- function(df,coastline) {
  df$in_land =0
  for( i in unique(coastline$id)) {
    my_sub <- subset(coastline,id==i)
    myind <- which(pnt.in.poly(cbind(df$x,df$y),my_sub[c("x","y")])$pip == 1)
    df$in_land[myind] <- i
  }
  return(df)
}

Add_model <- function(Res,Model_df,Sphere_triang,name,r=1) {
  Model3d <- as.data.frame(cbind(Convert_to_xyz(cbind(Model_df$x,Model_df$y),r),val=Model_df$z))
  Res[name]=0
  for (i in 1:nrow(Model3d)) {
    p <- Sphere_triang@pars$p
    dist <- sqrt((p[,1] - Model3d$x[i])^2 + (p[,2] - Model3d$y[i])^2 + (p[,3] - Model3d$z[i])^2)
    Res$Model[which.min(dist)] = Model3d$val[i]
  }
  return(Res)
}

Read_model_data <- function(modelpath,indexpath,model_grid,plot_it=T) {
  Model_data <- cbind(read.table(modelpath,col.names="z"),
                      read.table(indexpath,col.names="n"))
  
  grid_x <- seq(model_grid[1],model_grid[2],model_grid[3])
  grid_y <- seq(model_grid[4],model_grid[5],model_grid[6])
  grid <- data.frame(expand.grid(grid_x,grid_y))
  names(grid) = c("x","y")
  grid$n <- 1:nrow(grid)
  grid$mask <- rep(0,nrow(grid))
  mask_id <- Model_data$n
  grid$mask[mask_id] <- 1
  Model_df <- merge(grid,Model_data,by="n",all.x=T)
  if (plot_it){
     g <- LinePlotTheme() + geom_tile(data=Model_df,aes(x,y,fill=pmin(pmax(z,-5),5))) + 
       scale_fill_gradient2(low=muted("blue"),mid="light yellow",high=muted("red"),limits=c(-5,5)) + coord_fixed()
     print(g)
  }
      
  #Convert to radians
  #Model_df$x <- (Model_df$x)*2*pi/360
  #Model_df$y <- (Model_df$y)*2*pi/360
  
  return(Model_df)
}


Load_sphere_mesh <- function() {
  
  p  <- (as.matrix(read.table('./Code/p.csv',sep=",")))
  tri <- as.matrix(read.table('./Code/tri.csv',sep=","))
  M <-  as.matrix(read.table('./Code/Sphere1_M.csv',sep=","))
  K <-  as.matrix(read.table('./Code/Sphere1_K.csv',sep=","))
  n <- nrow(p)
  M <- sparseMatrix(i=M[,1],j=M[,2],x=M[,3],dims=c(n,n))
  K <- sparseMatrix(i=K[,1],j=K[,2],x=K[,3],dims=c(n,n))
  
  Mesh <- initFEbasis(p=p, t = tri, M = M, K = K)
  
  return(Mesh)
  
}

Load_egg_mesh <- function() {
  
  p  <- (as.matrix(read.table('./Code/p_egg.csv',sep=",")))
  tri <- as.matrix(read.table('./Code/tri_egg.csv',sep=","))
  M <-  as.matrix(read.table('./Code/Egg_M.csv',sep=","))
  K <-  as.matrix(read.table('./Code/Egg_K.csv',sep=","))
  n <- nrow(p)
  M <- sparseMatrix(i=M[,1],j=M[,2],x=M[,3],dims=c(n,n))
  K <- sparseMatrix(i=K[,1],j=K[,2],x=K[,3],dims=c(n,n))
  
  Mesh <- initFEbasis(p=p, t = tri, M = M, K = K)
  
  return(Mesh)
  
}

Read_coastline <- function() {
  
  coastline <- readShapeSpatial("./Code/shapefiles/ne_110m_land.shp")
  coastline_fort <- fortify(coastline)
  coastline_fort$x <- coastline_fort$long*2*pi/360
  coastline_fort$y <- coastline_fort$lat*2*pi/360
  return(coastline_fort)
  
}




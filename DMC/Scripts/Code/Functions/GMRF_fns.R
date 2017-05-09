## Initialise new GMRF fron mean and precision matrix with permutation
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

## Initialise new GMRF as a Random Walk (first order)
GMRF_init.RW <- function(n = 10,order=1,precinc = 1,perm=F) {

    # Create a GMRF of a RW (first order)
    mu <- matrix(0,nrow=n,ncol=1)
    i = c(1:(n-1),1:(n-1))
    j = c(1:(n-1),2:n)
    x <- numeric(length=((n-1)*2))
    x[1:(n-1)] = -1
    x[n:((2*n)-2)] = 1
    Dmat = sparseMatrix(i,j,x=x)
    R = t(Dmat)%*%Dmat
    if (order == 1) {
    	Q = precinc*R
    	intrinsic = 1
    	}
    if (order == 2) {
    	R <- R %*% R
    	R[1,(1:3)] = c(1,-2,1)
    	R[2,(1:4)] = c(-2,5,-4,1)
    	R[(n-1),(n-3):n] = c(1,-4,5,-2)
    	R[(n),(n-2):n] = c(1,-2,1)
    	Q <- precinc*R
    	intrinsic = 2
    }
    obj <- GMRF_init(mu,Q,intrinsic = intrinsic,perm=perm)


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


FindC_EOF <- function(X,Obs,roundobs=1) { # Requires x,y and n (n to maintin order after merge)
  Obs$x <- round(Obs$x/roundobs)*roundobs
  Obs$y <- round(Obs$y/roundobs)*roundobs
  C <- merge(Obs,X,all.x=T)
  C <- arrange(C,desc(-n))
  C <- C[,!(names(C) %in% c("x","y","n"))]
  return(C)
}

## Like FindC_boxaverage2 butfor arbitrary polygons
FindC_polyaverage_EOF  <- function(X,polygons,mulfun=1)  {
  n <- dim(X)[2] - 2
  m <- length(polygons)
  XY <- X[,1:2]
  S <- X[,-(1:2)]
  SS <- S * mulfun
  C <- matrix(0,m,n)
  
  # For each box
  for (i in 1:m) {
    # Find points which lie in box
    pv <- poly_xy(polygons[i])
    xy_in_poly <- pnt.in.poly(XY,pv)$pip
    C[i,] <- apply((SS[which(xy_in_poly==1),]),2,sum)
  }
  return(as(C,"dgCMatrix"))
}



## Find C matrix when observations are isolated points
FindC <- function(p,tri,locs,method="R") {
   if (length(locs[[1]]) > 0)  {
     t_num <- tsearch(p[,1], p[,2], tri, locs[[1]], locs[[2]], bary = FALSE)
     z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
     b <- matrix(0,3,1)
     A <- matrix(0,3,3)
     A[,1] = 1
     
     
     if(method=="C") {
       dyn.load("./Functions/element_interp/R_WINDOWS/element_interp.dll")
       nt <- length(t_num)
       X <- .C("element_interp",as.integer(as.integer(nt)),as.integer(t_num),as.integer(tri[,1]),
              as.integer(tri[,2]),as.integer(tri[,3]),as.double(p[,1]),as.double(p[,2]),
              as.double(locs[[1]]),as.double(locs[[2]]),i_ind1 = double(nt), i_ind2 = double(nt),
              i_ind3 = double(nt), j_ind1 = double(nt), j_ind2 = double(nt),
              j_ind3 = double(nt),z1 = double(nt), z2 = double(nt),z3 = double(nt))
       i_ind <- cbind(X$i_ind1,X$i_ind2,X$i_ind3)
       j_ind <- cbind(X$j_ind1,X$j_ind2,X$j_ind3)       
       z <- cbind(X$z1,X$z2,X$z3)       
       
       dyn.unload("./Functions/element_interp/R_WINDOWS/element_interp.dll")
     } else {
     
         for (i in 1:length(t_num)) {
            t_num_i <- t_num[i]
            this_tri <- tri[t_num_i,]
            this_p <- p[this_tri,]
            A[,2:3] <- this_p
            Ainv <- solve(A) 
    
    #        A <- matrix(c(1,1,1,p[tri[t_num_i,1],1],p[tri[t_num_i,2],1],
    #                 p[tri[t_num_i,3],1],p[tri[t_num_i,1],2],
    #                 p[tri[t_num_i,2],2],p[tri[t_num_i,3],2]),3,3)
           
               for (j in 1:3) {
                  b[,] <- 0
                  b[j] <- 1
                  i_ind[i,j] <- i
                  j_ind[i,j] <- this_tri[j]
                  z[i,j] <- matrix(c(1,locs[[1]][i],locs[[2]][i]),1,3)%*%Ainv%*%b
               }
         }
    }
    
    C <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                    dims = c(length(locs[[1]]),dim(p)[1]))
    return(C)
    }  else {
    return(matrix(1,0,0))
    }

}


 ## Like FindC_boxaverage but instead uses numerical integration to find int_phi(s) more accurately
 FindC_boxaverage2  <- function(p,tri,locs,boxsize,plotit=F,mulfun = 0,method="R")    {
  if (plotit == T) dev.new()
  n <- dim(p)[1]
  m <- length(locs[[1]])
  C_j <- C_i <- C_z <- NULL
  ds = 80
  # For each box
  for (i in 1:m) {
      # Creates fine grid
      x_grid <- linspace(locs[[1]][i]- boxsize[[i]][1]/2,locs[[1]][i] + boxsize[[i]][1]/2,ds)
      y_grid <- linspace(locs[[2]][i]- boxsize[[i]][2]/2,locs[[2]][i] + boxsize[[i]][2]/2,ds)
      dx <- mean(diff(x_grid))
      dy <- mean(diff(y_grid))
      GRID <- meshgrid(x_grid,y_grid)
      x <- as.vector(GRID$x)
      y <- as.vector(GRID$y)

      # Now for each grid point find the triangle it falls in
      t_num <- tsearch(p[,1], p[,2], tri, x, y, bary = FALSE)
      t_num <- t_num[!is.na(t_num)]   # Change this to cater for boundary effects later on!
      z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
      if (class(mulfun) == 'function') {
          mul_vals <- cbind(mulfun(p[tri[t_num,1],]),
                            mulfun(p[tri[t_num,2],]),
                            mulfun(p[tri[t_num,3],]))
      } else {
          mul_vals = 1
      }
       
       if(method=="C") {
         dyn.load("./Functions/element_interp/R_WINDOWS/element_interp.dll")
           nt <- length(t_num)
           X <- .C("element_interp",as.integer(as.integer(nt)),as.integer(t_num),as.integer(tri[,1]),
                  as.integer(tri[,2]),as.integer(tri[,3]),as.double(p[,1]),as.double(p[,2]),
                  as.double(x),as.double(y),i_ind1 = double(nt), i_ind2 = double(nt),
                  i_ind3 = double(nt), j_ind1 = double(nt), j_ind2 = double(nt),
                  j_ind3 = double(nt),z1 = double(nt), z2 = double(nt),z3 = double(nt))
           i_ind <- cbind(X$i_ind1,X$i_ind2,X$i_ind3)
           j_ind <- cbind(X$j_ind1,X$j_ind2,X$j_ind3)       
           z <- cbind(X$z1,X$z2,X$z3) 
       } else { 
       
          # Now see the magnitude for the three basis it touches
          for (j in 1:length(t_num))      {
             A <- matrix(c(1,1,1,p[tri[t_num[j],1],1],p[tri[t_num[j],2],1],
                     p[tri[t_num[j],3],1],p[tri[t_num[j],1],2],
                     p[tri[t_num[j],2],2],p[tri[t_num[j],3],2]),3,3)
    
               for (k in 1:3) {
                  b <- matrix(0,3,1)
                  b[k] <- 1
                  i_ind[j,k] <- j
                  j_ind[j,k] <- tri[t_num[j],k]
                  z[j,k] <- matrix(c(1,x[j],y[j]),1,3)%*%solve(A)%*%b
               }
            }
      }
      z <- z*mul_vals  
      Mmat <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                    dims = c(length(x),n))
      Cm <- colSums(Mmat)*dx*dy
      not_zero <- which(abs(Cm)>0)
      C_j <- c(C_j,not_zero)
      C_i <- c(C_i,rep(i,length(not_zero)))
      C_z <- c(C_z,Cm[not_zero])


  }
  C <- sparseMatrix(C_i,C_j,x = C_z,dims = c(length(locs[[1]]),dim(p)[1]))
  return(C)
 }

 ## Like FindC_boxaverage2 butfor arbitrary polygons. Here ds is the number of points to use for integration
 FindC_polyaverage  <- function(p,tri,polygons,plotit=F,mulfun = 0,method="R",ds=400)  {
  if (plotit == T) dev.new()
  n <- dim(p)[1]
  m <- length(polygons)
  C_j <- C_i <- C_z <- NULL
  
  # Find radius to consider
  max_tri_length <- max(apply(tri,1,function(x) max(rdist(p[x,],p[x,]))))
  
  # For each box
  for (i in 1:m) {
      # Creates fine grid
      pv <- poly_xy(polygons[i])
      pnts_to_consider <- which(apply(rdist(p,pv),1,min) < max_tri_length)
      pnts <- p[pnts_to_consider,]
      tris_to_consider <- which(apply(tri,1,function(x) any(x %in% pnts_to_consider)))
      tris <- tri[tris_to_consider,]
      
      x_grid <- linspace(min(pv[,1]),max(pv[,1]),ds)
      y_grid <- linspace(min(pv[,2]),max(pv[,2]),ds)
      dx <- mean(diff(x_grid))
      dy <- mean(diff(y_grid))
      GRID <- meshgrid(x_grid,y_grid)
      x <- as.vector(GRID$x)
      y <- as.vector(GRID$y)
      xy_in_poly <- pnt.in.poly(cbind(x,y),pv)$pip
      x <- x[which(xy_in_poly == 1)]
      y <- y[which(xy_in_poly == 1)]

      # Now for each grid point find the triangle it falls in
      t_num <- tsearch(p[,1], p[,2], tris, x, y, bary = FALSE)
      t_num <- t_num[!is.na(t_num)]   # Change this to cater for boundary effects later on!
      z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
      # Now see the magnitude for the three basis it touches
       if (class(mulfun) == 'function') {
          mul_vals <- cbind(mulfun(p[tris[t_num,1],]),
                            mulfun(p[tris[t_num,2],]),
                            mulfun(p[tris[t_num,3],]))
      } else {
          mul_vals = 1
      }
     
      
       if(method=="C") {
         dyn.load("./Functions/element_interp/R_WINDOWS/element_interp.dll")
           nt <- length(t_num)
           X <- .C("element_interp",as.integer(as.integer(nt)),as.integer(t_num),as.integer(tris[,1]),
                  as.integer(tris[,2]),as.integer(tris[,3]),as.double(p[,1]),as.double(p[,2]),
                  as.double(x),as.double(y),i_ind1 = double(nt), i_ind2 = double(nt),
                  i_ind3 = double(nt), j_ind1 = double(nt), j_ind2 = double(nt),
                  j_ind3 = double(nt),z1 = double(nt), z2 = double(nt),z3 = double(nt))
           i_ind <- cbind(X$i_ind1,X$i_ind2,X$i_ind3)
           j_ind <- cbind(X$j_ind1,X$j_ind2,X$j_ind3)       
           z <- cbind(X$z1,X$z2,X$z3) 
       } else {
        for (j in 1:length(t_num))      {
             A <- matrix(c(1,1,1,p[tris[t_num[j],1],1],p[tris[t_num[j],2],1],
                     p[tris[t_num[j],3],1],p[tris[t_num[j],1],2],
                     p[tris[t_num[j],2],2],p[tris[t_num[j],3],2]),3,3)
    
               for (k in 1:3) {
                  b <- matrix(0,3,1)
                  b[k] <- 1
                  i_ind[j,k] <- j
                  j_ind[j,k] <- tris[t_num[j],k]
                  z[j,k] <- matrix(c(1,x[j],y[j]),1,3)%*%solve(A)%*%b
                  if (class(mulfun) == 'function') {
                    
                  }
               }
          }   
      }
      z <- z*mul_vals              
      Mmat <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                    dims = c(length(x),n))
      Cm <-   colSums(Mmat)*dx*dy
      not_zero <- which(abs(Cm)>0)
      C_j <- c(C_j,not_zero)
      C_i <- c(C_i,rep(i,length(not_zero)))
      C_z <- c(C_z,Cm[not_zero])
 
  }
  C <- sparseMatrix(C_i,C_j,x = C_z,dims = c(m,n))
  return(C)
 }


## Find C matrix when observations are box averages
FindC_boxaverage  <- function(p,tri,locs,boxsize,plotit=F)  {
  if (plotit == T) dev.new()
  n <- dim(p)[1]
  m <- length(locs[[1]])
  Voronoi <- deldir(p[,1],p[,2],plotit=plotit,sort=F,rw=c(min(p[,1])-0.00001,max(p[,1])+0.00001,min(p[,2])-0.00001,max(p[,2])+.00001))
	Tess_polys <-  vector("list",n)
  for (i in 1:n)  {
    edges <- rbind(as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x1,y1,bp1))),
                   as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x1,y1,bp1))),
                   as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x2,y2,bp2))),
                   as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x2,y2,bp2))))
    X <- unique(edges)
    if(sum(X[,3])>0) {
     X <- rbind(X,c(p[i,],0))
    }
    edges <- X[,(1:2)]
    edges <- edges[chull(edges), ]
    Tess_polys[[i]] <- as(edges,"gpc.poly")
  }
  count <- 0
  z_ind <- j_ind <-i_ind <- rep(0,10*m)

  for (i in 1:m) {
          # Coordinates of box
          box_obs <- matrix(c(locs[[1]][i]-boxsize[[i]][1]/2,
                           locs[[1]][i]+boxsize[[i]][1]/2,
                           locs[[1]][i]+boxsize[[i]][1]/2,
                           locs[[1]][i]-boxsize[[i]][1]/2,
                           locs[[2]][i]-boxsize[[i]][2]/2,
                           locs[[2]][i]-boxsize[[i]][2]/2,
                           locs[[2]][i]+boxsize[[i]][2]/2,
                           locs[[2]][i]+boxsize[[i]][2]/2),4,2)
         poly_obs <- as(box_obs, "gpc.poly")
         # Find states (vertices) which lie within box
         in_points <- pnt.in.poly(p,box_obs)

         # This gives us an indication which states to consider
         if (sum(in_points$pip) == 0)
         {  # Probably mesh is very course so do exhaustive overlap search
             pnts_to_check <- 1:n

          }
          else {
             pnts_to_check <- which(in_points$pip == 1)
             add_points <- NULL
             # Find in which triangles these points lie and get their neighbours
             for (j in pnts_to_check) {
                add_points <- c(add_points,array(tri[which(tri == j,arr.ind=T)[,1],]))
                }
                pnts_to_check <- unique(add_points)
             }

         for (j in sort(pnts_to_check)) {
            area_this_state <- area.poly(intersect(poly_obs,Tess_polys[[j]]))
            if (area_this_state>0) {
                count <- count + 1
                i_ind[count] <- i
                j_ind[count] <- j
                z_ind[count] <- area_this_state
            }
       }
  }
  C <- sparseMatrix(i_ind[1:count],j_ind[1:count],x=z_ind[1:count],dims = c(m,n))
  return(C)
 }

## Test if object exists
testObject <- function(object)
{
   exists(as.character(substitute(object)))
}

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


## Create a precision matrix from neighbourhood list
Prec_from_neighb <- function(neighb,intrinsic=1,precinc=1){

	num_v <- length(neighb)
	num_neighb <- lapply(neighb,length)

	if (intrinsic == 1)  {
     i_list <-  vector("list",num_v)    # row number index
     for (k in 1:num_v) {
    	  i_list[[k]] <- rep(k,num_neighb[[k]])
     }
     i <- unlist(i_list)
	   j <- unlist(neighb)             # column number index

     z <- rep(-1, length(j))  # for intrinsic = 1 all non-zero off-diagonal elements are -1

     # Now add the diagonal elements (= number of neighbours
     i <- c(i,1:num_v)
     j <- c(j,1:num_v)
     zdiag <- unlist(num_neighb)
     #zdiag[zdiag == 0] = 1   # Make lonely elements conditionally independent with finite variance
     z <- c(z,zdiag)
	  }

  if (intrinsic == 2)  {

     # Start off with the diagonal elements
     i1 <- 1:num_v
     j1 <- 1:num_v
     z1 <- rep(0,num_v)
     for (k in 1:num_v) {
        z1[k] <- num_neighb[[k]]^2 + num_neighb[[k]]
     }

     # Now add elements corresponding to the neighbours. Initialise, prob max 10 neighbours per node
     count <- 1
     i2 <- rep(0,num_v*10)
     j2 <- rep(0,num_v*10)
     z2 <- rep(0,num_v*10)
     for (k in 1:num_v) {
         for (l in neighb[[k]]) {
           i2[count] <- k
           j2[count] <- l
           z2[count] <-  -(num_neighb[[k]] +  num_neighb[[l]] -  sum(duplicated(c(neighb[[l]],neighb[[k]]))))
           count <- count + 1
           }
      }
      i2 <- i2[1:count-1]
      j2 <- j2[1:count-1]
      z2 <- z2[1:count-1]

     # Now add elements corresponding to the neighbours of the neighbours. Initialise, prob max 15 neighbours neighbours per node
     count <- 1
     i3 <- rep(0,num_v*15)
     j3 <- rep(0,num_v*15)
     z3 <- rep(0,num_v*15)
     neighb2 <- vector("list",num_v)
     for (k in 1:num_v) {
       # Construct neighbours of neighbours list, with redundancies (the value of the redundancy is then the element in Q)
       for (l in neighb[[k]]) {
         neighb2[[k]] <- c(neighb2[[k]],setdiff(neighb[[l]],c(neighb[[k]],k)))   # return number of elements in neighb[[l]] which are not in neighb[[k]] ( + node under consideration)
       }
       for (l in unique(neighb2[[k]]))  {
          i3[count] <- k
          j3[count] <- l
          z3[count] <- sum(neighb2[[k]] == l)
          count <- count + 1
       }
      }

     i3 <- i3[1:count-1]
     j3 <- j3[1:count-1]
     z3 <- z3[1:count-1]

     i <- c(i1,i2,i3)
     j <- c(j1,j2,j3)
     z <- c(z1,z2,z3)

    }



 z <- precinc*z
 Q <- sparseMatrix(i,j,x=z)
 return(Q)
}

# Create a covariance function from a set of locations
Covariance_from_points <- function(fn,pars,locs)  {
   n <- length(locs[[1]])
   V <- matrix(0,n,n)
   for(i in 1:n)
     for (j in 1:i)  {
        x1 = locs[[1]][i]
        y1 = locs[[2]][i]
        x2 = locs[[1]][j]
        y2 = locs[[2]][j]
        Dx <- x1-x2
        Dy <- y1-y2
        Dist <- sqrt(Dx^2 + Dy^2)
       if (fn == "sqexp") {
          V[j,i] <- V[i,j] <- pars[1]*exp(-Dist/pars[2])
       } else if (fn == "truncsqexp") {
          tau <- sqrt(2*pi/pars[2])
          if (tau*Dist > 2*pi) {
             V[j,i] <- V[i,j] <- 0
          } else {
             V[j,i] <- V[i,j] <- pars[1]/(3*pi)*((2*pi - tau*Dist)*(1 + (cos(tau*Dist))/2)  +
                                 3*sin(tau*Dist)/2)
          }

        }

      }
     return(V)
  }

## Generate a random graph with num_v vertices with a mean of lambda neighbours each
Random_graph <- function(num_v,lambda) {
	neighb <- vector("list",num_v)
	for (i in c(1:num_v)) {
		num_neighb <- rpois(1,lambda = lambda)
		numattempts = 0
		while ((length(neighb[[i]]) <  num_neighb) && (numattempts < 10)) {
			elect_n <- floor(runif(1,min=i,max=num_v+1))
			numattempts <- numattempts + 1
			if (!(elect_n %in% neighb[[i]]) && i != elect_n) {
				neighb[[i]] <- c(neighb[[i]],elect_n)
				neighb[[elect_n]] <- c(neighb[[elect_n]],i)
				}

			}
	}
	return(neighb)
}

## Return neighb list from precision structure
neighb_from_prec <- function(Q) {
 n <- dim(Q)[1]
 return(apply(matrix(1:n),1,function(x) { return(which(Q[x,]>0))}))


}

 

## Create a graph from a triangulation
Graph_from_tri <- function(p,tri) {

  neighb <- vector("list",dim(p)[1])
  for (i in 1:dim(tri)[1])  {
    this_tri = tri[i,]
    for (j in 1:3) {
        for (k in 1:3) {
              if (!(this_tri[k] %in% neighb[[this_tri[j]]]) && j != k)
				      neighb[[this_tri[j]]] <- c(neighb[[this_tri[j]]],this_tri[k])
				}
    }
  }
  return(neighb)
  }

 ## Create a graph from pixels
Graph_from_grid <- function(x,y,d) {
       n = length(x)
       X <- vector("list",n)
       for (i in 1:n) {
        neighb1 = intersect(which(x == x[i]),which(y == (y[i]+d)))
        neighb2 = intersect(which(x == x[i]),which(y == (y[i]-d)))
        neighb3 = intersect(which(x == (x[i]+d)),which(y == (y[i])))
        neighb4 = intersect(which(x == (x[i]-d)),which(y == (y[i])))
        X[[i]] <- c(neighb1,neighb2,neighb3,neighb4)

       }
       return(X)
      }


## Explore a posterior function
explore_theta <- function(Covariance,log_theta_post,max_log_theta,dz=0.5,prune_grid = T) {
   # Standardize
   X <- eigen(Covariance,symmetric=T)
   ind <- sort(sort(diag(Covariance),decreasing=T,index.return=T)$ix,index.return=T)$ix
   X$values <- X$values[ind]
   X$vectors <- X$vectors[,ind]

   n <- dim(Covariance)[1]
   zexplore = matrix(rep(0,n))
   z_list <- vector("list",n)
   log_theta_post_max <- log_theta_post(max_par)
   if (n > 1) {
   VsqrtLambda <- X$vectors%*%sqrt(diag(X$values))
   } else {
   VsqrtLambda <- X$vectors%*%sqrt(X$values)
   }

   cat("Finding axial points...",sep="\n")
   flush.console()
   for (i in 1:n)    {
     while((log_theta_post_max - log_theta_post(c(max_par + VsqrtLambda%*%zexplore))) < 2.5) {
      z_list[[i]] <- c(z_list[[i]],zexplore[i])
      zexplore[i] <- zexplore[i] + dz
      }
      zexplore[i] = -dz
      while(((log_theta_post_max - log_theta_post(c(max_par + VsqrtLambda%*%zexplore))) < 2.5)
          &&(max_par[i] + (VsqrtLambda%*%zexplore)[i] > 0)) {
      z_list[[i]] <- c(z_list[[i]],zexplore[i])
      zexplore[i] <- zexplore[i] - dz
     }
     zexplore[i] <- 0
   }
   for (i in 1:n)    {
     z_list[[i]] <- sort(z_list[[i]])
    }

   cat("Forming grid...",sep="\n")
   flush.console()
   if (n == 1) {
       z = matrix(unlist(z_list))
   } else if (n == 2) {
      Z <- meshgrid(z_list[[1]],z_list[[2]])
      z <- matrix(c(c(Z$x),c(Z$y)),length(Z$x),2)
  # } else if (n == 3) {   #3d meshgrid NOT WORKING!!
  #    Z <- meshgrid(z_list[[1]],z_list[[2]],z_list[[3]])
  #    z <- matrix(c(c(Z$x),c(Z$y),c(Z$z)),length(Z$x),3)
   } else {
     nlist <- unlist(lapply(z_list,length))
     z <- matrix(0,prod(nlist),length(nlist))
     z[,length(nlist)] <- rep(z_list[[length(nlist)]],prod(nlist[1:(length(nlist)-1)]))
     for(i in seq(length(nlist),1,-1)) {
       if (i==1) {
           repmul = 1
        } else {
           repmul = prod(nlist[1:(i-1)])
        }
        temp <- rep(z_list[[i]],repmul)

        if ((i+1) > length(nlist)) {
           repmul = 1
        } else {
           repmul = prod(nlist[(i+1):length(nlist)])
        }

        whole_vec <- rep(temp,repmul)
        A <- reshape(as.matrix(whole_vec),length(temp),repmul)
        z[,i] = reshape(t(A),prod(nlist),1)
     }
   }

   theta_explore = t(max_par +  VsqrtLambda%*%t(z))
   d_theta <- VsqrtLambda%*%matrix(rep(dz,n))

   # Now we have a list of points on a grid and we can remove those which are not important
   cat("Computing log-posterior at all points and pruning...",sep="\n")
   flush.console()
   ind_keep <- NULL
   log_theta_post_vals <- rep(0,dim(theta_explore)[1])
   for (i in 1:length(theta_explore[,1])) {
      log_theta_post_vals[i] <- log_theta_post(theta_explore[i,])
      if ((log_theta_post_max - log_theta_post_vals[i]) < 2.5)
        ind_keep <- c(ind_keep,i)
      cat(paste(i/length(theta_explore[,1])*(100),"% complete"),sep="\n")
      flush.console()
   }


   if (prune_grid == T) {
      theta_explore <- theta_explore[ind_keep,]
      log_theta_post_vals <- log_theta_post_vals[ind_keep]
   }

   # Unnormalised posterior
   theta_post_unnorm <- exp(log_theta_post_vals - max(log_theta_post_vals))

   # Normalised posterior in theta scale
   vol <- abs(sum(theta_post_unnorm)*abs(prod(d_theta)))
   theta_post_norm <- theta_post_unnorm/vol

   # Plot marginals
   # Normalised posterior in z scals
   vol <- abs(sum(theta_post_unnorm)*abs(dz^n))
   theta_post_norm_z <- theta_post_unnorm/vol

   #for(i in 1:n) {
   # dev.new()
   # zaxis <- z[,i]
   # zmarg <- unlist(lapply(split(theta_post_norm_z,zaxis),function(x){sum(x*dz^(n-1)) }),use.names='F')
   # thetaaxis <- max_par[i] + VsqrtLambda[i,i]*unique(zaxis)
   # vol <- abs(sum(zmarg)*d_theta[i])
   # theta_marg_norm <- zmarg/vol
   # plot(thetaaxis,theta_marg_norm,type="o",xlab="theta_i",ylab="p(theta_i)")
   #}

  #theta_axis_all <- matrix(max_par,3,dim(z)[1]) + VsqrtLambda%*%t(z)
  thetaaxis <- vector("list",n)
  theta_marg_norm <- vector("list",n)
  for(i in 1:n) {
    #dev.new()
    zaxis <- z[,i]
    zmarg <- unlist(lapply(split(theta_post_norm_z,zaxis),function(x){sum(x*dz^(n-1)) }),use.names='F')

    zmat <- matrix(0,length(unique(zaxis)),n)
    zmat[,i] <- unique(zaxis)
    thetaaxis[[i]] <- (max_par + VsqrtLambda%*%t(zmat))[i,]
    vol <- abs(sum(zmarg)*d_theta[i])
    theta_marg_norm[[i]] <- zmarg/vol
    plot(thetaaxis[[i]],theta_marg_norm[[i]],type="o",xlab="theta_i",ylab="p(theta_i)")
   }


   return(list(post_unnorm = theta_post_unnorm,post_norm = theta_post_norm,
               theta_explore = theta_explore,d_theta = d_theta,
               thetaaxis = thetaaxis, theta_marg_norm = theta_marg_norm))

}



 optimset <- function()
 return (list(maxit = 100L, tol=1e-7))

 myoptim <- function(init,fn,gr,hessian,control=list(),aggressive=100,pars=NULL) {
    
   if (!is.null(pars)) {
    fn_wrap <- function(theta){ 
       return(fn(theta,pars))
    }
    gr_wrap <- function(theta){ 
      return(gr(theta,pars))
    }
   }
   
     # Establish control parameter
      con <- optimset()
      nmsC <- names(con)
      con[(namc <- names(control))] <- control
      if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

  

     delta <- Inf
     theta <- init
     gradient <- gr_wrap(theta)
     thetanew <- init
     alpha = rep(0,length(theta))
     stepsize <- rep(0.01,length(theta))
     itercount <- 0
     while((abs(delta) > con$tol) && (stepsize > 1e-11) && (itercount < con$maxit)) {
         gradient <- gr_wrap(theta)
         for (i in 1:length(theta))   {
            alpha[i] <- gradient[i]
            # Find appropriate stepsize
           if (abs(alpha[i]) < stepsize[i])  { while (abs(alpha[i]) < stepsize[i]) alpha[i] <- alpha[i]*10
              } else if (abs(alpha[i])>(10*stepsize[i])) {while (abs(alpha[i])> (10* stepsize[i])) alpha[i] <- alpha[i]/10 }
          #thetanew[i] <- theta[i] - alpha[i]*abs(theta[i])
          thetanew[i] <- theta[i] - alpha[i]*stepsize[i]*aggressive
          }

          delta = as.vector(fn_wrap(thetanew) - fn_wrap(theta))
          if (delta < 0) {
              theta = thetanew
              stepsize <- rep(0.01,length(theta))
          } else {            # Overshoot detected
          stepsize <- stepsize/10
          }

          itercount <- itercount + 1
          cat(paste("Step ",itercount,": delta = ",delta,"theta = ",do.call(paste,as.list(theta))),sep="\n")
          flush.console()
     }

     if (abs(delta) < con$tol) cat("Reached local minimum",sep="\n")
     if (stepsize[1] < 1e-11) cat("Reached minimum stepsize (local minimum)",sep="\n")
     if (itercount == con$maxit) cat("Maximum iteration count reached",sep="\n")

     if (hessian == T) {
        H <- matrix(0,length(theta),length(theta))

        # get decent stepsizes in alpha (as above)
        stepsize <- 0.001
        for (i in 1:length(gradient)) {
        alpha[i] <- gradient[i]
        if (abs(alpha[i]) < stepsize)  { while (abs(alpha[i]) < stepsize) alpha[i] <- alpha[i]*10
        } else if (abs(alpha[i])>stepsize) {while (abs(alpha[i])> stepsize) alpha[i] <- alpha[i]/10 }
        }

        stepsizes <- abs(alpha*theta)


        for (i in 1:length(theta))
          for (j in 1:(i)) {
           maskvecj <- maskveci <- rep(0,length(theta))
           maskvecj[j] = stepsizes[j]
           maskveci[i] = stepsizes[i]
           gradient = 0.5*((gr_wrap(theta + .5*maskvecj)[i] - gr_wrap(theta - .5*maskvecj)[i])/(stepsizes[j]) +
                          (gr_wrap(theta + .5*maskveci)[j] - gr_wrap(theta - .5*maskveci)[j])/(stepsizes[i]) )
           H[j,i] <- H[i,j] <- gradient

          }

     } else { H = NA}
     return(list(par = theta,hessian = H))
    }

marg_prec_from_kappa <- function(kappa_l,nu) {
  return(1/(gamma(nu)/(gamma(nu+1)*4*pi*kappa_l^(2*nu))))
}

  Prec_from_SPDE <- function(M,K,tau,kappa,alpha=1)  {
    if (!(alpha %in% c(1,2,3,4))) {
     stop("alpha > 2 not implemented yet")
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


my_Matern <- function(r=0:100,nu=3/2,var=1,kappa=0.1) {
     K <- var/((2^(nu-1))*gamma(nu))*(kappa*abs(r))^nu*besselK(kappa*abs(r),nu=nu)
     if (class(K) == "matrix") {
       diag(K) = var
     }
     return(K)
}

my_RBF <- function(r,mu=matrix(0,1,2),sigma2=1,A=1) {
     y <- A*exp(-r^2/(2*sigma2))
     return(y)
}

Build_AQA <- function(Qx,A,T) {
  # For now T >= 3
  n <- nrow(Qx)
  QA <- Qx %*% A
  AQA <- A %*% Qx %*% A + Qx
  AQ <- t(A) %*% Qx
  for ( i in 0:(T-3)) {
    if (i == 0) {
      Q <- cBind(-QA, AQA, -AQ , Zeromat(n,((T-3)-i)*n))
    } else if (i == (T-3)) {
      Q <- rBind(Q,cBind(Zeromat(n,n*i),-QA, AQA, -AQ))
    } else {
      Q <- rBind(Q,cBind(Zeromat(n,n*i),-QA, AQA, -AQ , Zeromat(n,((T-3)-i)*n)))
    }
  }
  Q <- rBind(cBind(AQA, -AQ, Zeromat(n,(T-2)*n)),
             Q,
             cBind(Zeromat(n,n*(T-2)),-QA, Qx))
  return(Q)
}


# Distance function for polygon - taken from p_poly_dist (MATLAB)
  mesh.dpoly <-  function(p,pv) {

  p1 <- p[,1]
  p2 <- p[,2]
  xv <- pv[,1]
  yv <- pv[,2]

  # If (xv,yv) is not closed, close it.
  xv <- c(xv);
  yv <- c(yv);
  Nv = length(xv);
  if ((xv[1] != xv[Nv]) || (yv[1] != yv[Nv])) {
      xv <- c(xv,xv[1])
      yv <- c(yv,yv[1])
      Nv <- Nv + 1;
  }

  # linear parameters of segments that connect the vertices
  A <- -diff(yv);
  B <-  diff(xv);
  C <- yv[-1]*xv[-Nv] - xv[-1]*yv[-Nv];

  dvec <- rep(0,length(p1))

  for (i in 1:length(p1)) {
	  x <- p1[i]
	  y <- p2[i]

	  # find the projection of point (x,y) on each rib
	  AB <- 1/(A^2 + B^2);
	  vv <- (A*x+B*y+C);
	  xp <- x - (A*AB)*vv;
	  yp <- y - (B*AB)*vv;

	  # find all cases where projected point is inside the segment
	  idx_x <- (((xp>=xv[-Nv]) & (xp<=xv[-1])) | ((xp>=xv[-1]) & (xp<=xv[-Nv])));
	  idx_y <- (((yp>=yv[-Nv]) & (yp<=yv[-1])) | ((yp>=yv[-1]) & (yp<=yv[-Nv])));
	  idx <- idx_x & idx_y;

	  # distance from point (x,y) to the vertices
	  dv <- sqrt((xv[-Nv]-x)^2 + (yv[-Nv]-y)^2);

	  if(!any(idx)) {# all projections are outside of polygon ribs
	     d <- min(dv);
	  } else {
	     # distance from point (x,y) to the projection on ribs
	     dp <- sqrt((xp[idx]-x)^2 + (yp[idx]-y)^2);
	     d <- min(min(dv), min(dp));
	  }


	  if (in.poly(matrix(c(x,y),1,2),cbind(xv,yv))) {
	     d = -d
	    }
	dvec[i] <- d
      }
     return(matrix(dvec))

}


PolygonfromVoronoi <- function(Voronoi,p) {
      n = dim(p)[1]
      polygons <- vector("list",n)
      for (i in 1:n)  {
          #edges <- rbind(as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x1,y1,bp1))),
          #         as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x1,y1,bp1))),
          #         as.matrix(subset(Voronoi$dirsgs,ind1 ==i,select=c(x2,y2,bp2))),
          #         as.matrix(subset(Voronoi$dirsgs,ind2 ==i,select=c(x2,y2,bp2))))
          #X <- unique(edges)
          #if(sum(X[,3])>0) {
          #  X <- rbind(X,c(p[i,],0))
          #}
          #edges <- X[,(1:2)]
          
          X <- subset(Voronoi$dirsgs,(ind1 ==i | ind2 == i))
          X <- matrix(c(X$x1,X$x2,X$y1,X$y2,X$bp1,X$bp2),ncol=3)
          X <- unique(X)
          if(sum(X[,3])>0) {
            X <- rbind(X,c(p[i,],0))
          }
          edges <- X[,(1:2)]
          
          edges <- edges[chull(edges), ]
          polygons[[i]] <- as(edges,"gpc.poly")
      }
      return(polygons)

    }
    
    
inpolygon_MATLAB <- function(x,y,xv,yv) {
# vectorize the computation.
# Translated from MATLAB
inputSize = c(length(x),1);
Nv <- length(xv)
if (!any(is.nan(xv) | is.nan(yv))) {
    if ((xv[1] != xv[Nv]) || (yv[1] != yv[Nv])) {
        xv = c(xv,xv[1])
        yv = c(yv,yv[1])
        Nv = Nv + 1;
    }
}

mask = (x >= min(xv)) & (x <= max(xv)) & (y>=min(yv)) & (y<=max(yv));
inbounds = which(mask);
x = x[mask];
y = y[mask];

Np <- length(x)
x <- t(matrix(x,Np,Nv))
y <- t(matrix(y,Np,Nv))


# Compute scale factors for eps that are based on the original vertex
# locations. This ensures that the test points that lie on the boundary
# will be evaluated using an appropriately scaled tolerance.
# (m and mp1 will be reused for setting up adjacent vertices later on.)
m = 1:(Nv-1)
mp1 = 2:Nv
avx = abs(0.5*(  xv[m] + xv[mp1]))
avy = abs(0.5*(yv[m]+yv[mp1]))
scaleFactor = pmax(avx[m], avy[m])
scaleFactor = pmax(scaleFactor, avx[m]*avy[m])
# Translate the vertices so that the test points are
# at the origin.
xv = matrix(xv,Nv,Np) - x
yv = matrix(yv,Nv,Np) - y

# Compute the quadrant number for the vertices relative
# to the test points.
posX = xv > 0
posY = yv > 0
negX = !posX
negY = !posY
quad = (negX & posY) + 2*(negX & negY) + 3*(posX & negY)

# Ignore crossings between distinct edge loops that are separated by NaNs
nanidx = which(is.nan(xv) | is.nan(yv));
quad[nanidx] = NaN
# Compute the sign() of the cross product and dot product
# of adjacent vertices.
theCrossProd = xv[m,]* yv[mp1,] - xv[mp1,]* yv[m,]
signCrossProduct = sign(theCrossProd)


# Adjust values that are within epsilon of the polygon boundary.
# Making epsilon larger will treat points close to the boundary as
# being "on" the boundary. A factor of 3 was found from experiment to be
# a good margin to hedge against roundoff.
scaledEps = scaleFactor*.Machine$double.eps*3
idx = abs(theCrossProd) < scaledEps
signCrossProduct[idx] = 0
dotProduct = xv[m,]* xv[mp1,] + yv[m,]* yv[mp1,]

# Compute the vertex quadrant changes for each test point.
diffQuad = diff(quad)

# Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
# Any quadrant difference with an absolute value of 2 should have
# the same sign as the cross product.
idx = (abs(diffQuad) == 3)
diffQuad[idx] = -diffQuad[idx]/3
idx = (abs(diffQuad) == 2)
diffQuad[idx] = 2*signCrossProduct[idx]

# Find the inside points.
# Ignore crossings between distinct loops that are separated by NaNs
nanidx = which(is.nan(diffQuad))
diffQuad[nanidx] = 0
inpoly = (apply(diffQuad,2,sum) != 0)

# Find the points on the polygon.  If the cross product is 0 and
# the dot product is nonpositive anywhere, then the corresponding
# point must be on the contour.
on = apply(((signCrossProduct == 0) & (dotProduct <= 0)),2,any)
inpoly = inpoly | on

mask[inbounds[!inpoly]] = 0
return(c(reshape(matrix(mask),inputSize)))
}


    
inpolygon_MATLAB_big <- function(x,y,xv,yv) {
# vectorize the computation.
# Translated from MATLAB

x_full <- x
y_full <- y
xv_orig <- xv
yv_orig <- yv
concat_dist <- NULL
for (i in 1: ceil(length(x_full)/1000)) {
  this_batch <- ((i-1)*1000+1):min(c(i*1000,length(x_full)))
  x <- x_full[this_batch]
  y <- y_full[this_batch]
  xv <- xv_orig
  yv <- yv_orig
  
  inputSize = c(length(x),1);
  Nv <- length(xv)
  if (!any(is.nan(xv) | is.nan(yv))) {
      if ((xv[1] != xv[Nv]) || (yv[1] != yv[Nv])) {
          xv = c(xv,xv[1])
          yv = c(yv,yv[1])
          Nv = Nv + 1;
      }
  }
  
  mask = (x >= min(xv)) & (x <= max(xv)) & (y>=min(yv)) & (y<=max(yv));
  inbounds = which(mask);
  if( length(inbounds) == 0) {
    concat_dist <- rbind(concat_dist,matrix(0,length(this_batch),1))
  } else {
    x = x[mask];
    y = y[mask];
    
    Np <- length(x)
    x <- t(matrix(x,Np,Nv))
    y <- t(matrix(y,Np,Nv))
    
    
    # Compute scale factors for eps that are based on the original vertex
    # locations. This ensures that the test points that lie on the boundary
    # will be evaluated using an appropriately scaled tolerance.
    # (m and mp1 will be reused for setting up adjacent vertices later on.)
    m = 1:(Nv-1)
    mp1 = 2:Nv
    avx = abs(0.5*(  xv[m] + xv[mp1]))
    avy = abs(0.5*(yv[m]+yv[mp1]))
    scaleFactor = pmax(avx[m], avy[m])
    scaleFactor = pmax(scaleFactor, avx[m]*avy[m])
    # Translate the vertices so that the test points are
    # at the origin.
    xv = matrix(xv,Nv,Np) - x
    yv = matrix(yv,Nv,Np) - y
    
    # Compute the quadrant number for the vertices relative
    # to the test points.
    posX = xv > 0
    posY = yv > 0
    negX = !posX
    negY = !posY
    quad = (negX & posY) + 2*(negX & negY) + 3*(posX & negY)
    
    # Ignore crossings between distinct edge loops that are separated by NaNs
    nanidx = which(is.nan(xv) | is.nan(yv));
    quad[nanidx] = NaN
    # Compute the sign() of the cross product and dot product
    # of adjacent vertices.
    theCrossProd = xv[m,]* yv[mp1,] - xv[mp1,]* yv[m,]
    signCrossProduct = sign(theCrossProd)
    
    
    # Adjust values that are within epsilon of the polygon boundary.
    # Making epsilon larger will treat points close to the boundary as
    # being "on" the boundary. A factor of 3 was found from experiment to be
    # a good margin to hedge against roundoff.
    scaledEps = scaleFactor*.Machine$double.eps*3
    idx = abs(theCrossProd) < scaledEps
    signCrossProduct[idx] = 0
    dotProduct = xv[m,]* xv[mp1,] + yv[m,]* yv[mp1,]
    
    # Compute the vertex quadrant changes for each test point.
    diffQuad = diff(quad)
    
    # Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
    # Any quadrant difference with an absolute value of 2 should have
    # the same sign as the cross product.
    idx = (abs(diffQuad) == 3)
    diffQuad[idx] = -diffQuad[idx]/3
    idx = (abs(diffQuad) == 2)
    diffQuad[idx] = 2*signCrossProduct[idx]
    
    # Find the inside points.
    # Ignore crossings between distinct loops that are separated by NaNs
    nanidx = which(is.nan(diffQuad))
    diffQuad[nanidx] = 0
    inpoly = (apply(diffQuad,2,sum) != 0)
    
    # Find the points on the polygon.  If the cross product is 0 and
    # the dot product is nonpositive anywhere, then the corresponding
    # point must be on the contour.
    on = apply(((signCrossProduct == 0) & (dotProduct <= 0)),2,any)
    inpoly = inpoly | on
    
    mask[inbounds[!inpoly]] = 0
    concat_dist <- rbind(concat_dist,reshape(matrix(mask),inputSize))
    }
  }
return(c(concat_dist))
}

# Finds nearest neighbour of a continuous point from a set of discrete data points
interp_nn <- function(s,X) {
          s1 <- s[,1]
          s2 <- s[,2]
          x_all <- X[,1]
          y_all <- X[,2]
          z <- X[,3]

          # cater for nodes outside observation box
          outbox <- ((s1 < min(x_all)) | (s1 > max(x_all)) | (s2 < min(y_all)) |  (s2 > max(y_all)))
          if (length(s1)[1] == 1) {
            x <- z[which((abs(s1-x_all) == min(abs(s1-x_all)))&(abs(s2-y_all) == min(abs(s2-y_all))))[1]]  # Just in case it's in middle choose one corner
           } else {
            x <- apply(s,1,function(x) z[which((abs(x[1]-x_all) == min(abs(x[1]-x_all)))&
                                          (abs(x[2]-y_all) == min(abs(x[2]-y_all))))][1])
           }
           x[is.na(x)] <- 0
           return(x*!outbox)  
}

# interpolates a continuous data set onto a grid
nn_grid_interp <- function(s,df,delta=10,miss_value=NA) {
    s_rnd <- data.frame(round(s/delta)*delta)
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
  
# Regress gridded data on a triangulation
Regress_triang <- function(p,tri,z) {
    C_mean <- FindC(p,tri,z[,c(1,2)],method="C")
    EmptyC <-  which(colSums(C_mean) == 0)
    C_mean2 <- C_mean[,-EmptyC]
    n <- dim(p)[1]
    x_prior <- matrix(rep(0,n))
    x_prior[setdiff(1:n,EmptyC),] <- matrix(solve(t(C_mean2)%*%C_mean2)%*%t(C_mean2)%*%z[,3]) 
    return(x_prior)
    
    #C_mean <- FindC(p,tri,z[,c(1,2)],method="C")
    #EmptyC <-  which( apply(C_mean,2,sum) == 0)
    #z <- rbind(z,cbind(p[EmptyC,],0))
    #C_mean <- FindC(p,tri,z[,c(1,2)],method="C")
    #n <- dim(p)[1]
    #x_prior <- matrix(rep(0,n))
    #x_prior <- matrix(solve(t(C_mean)%*%C_mean)%*%t(C_mean)%*%z[,3]) 
    #return(x_prior)
}
  
  
  # Check if there is at least one corner of a square inside a polygon
  corner_in_poly <- function(cx,cy,lx,ly,poly) {
      return(which(!(inpolygon_MATLAB(cx+lx/2,cy+ly/2,poly[,1],poly[,2]) |  # Remove mascons out at sea
                          inpolygon_MATLAB(cx+lx/2,cy-ly/2,poly[,1],poly[,2]) |
                          inpolygon_MATLAB(cx-lx/2,cy+ly/2,poly[,1],poly[,2]) |
                          inpolygon_MATLAB(cx-lx/2,cy-ly/2,poly[,1],poly[,2]))))
  }

  # Check if polygons intersect
  poly_in_poly <- function(pol_list,poly) {
      x <- rep(0,length(pol_list))
      for (i in 1:length(pol_list)) {
        x[i] <- area.poly(intersect(pol_list[[i]],poly)) - area.poly(pol_list[[i]])
      }
      return(which(x<0))
  }  
  
#Order vertices
Order_vertices <- function(p) {
p2 <- p
for (i in 1:(nrow(p)-1)) {



    dist_mat = as.matrix(dist(p2));
    if ((i > 3) && (dist_mat[1,i] < min(dist_mat[i,-c(1:i)])))  {
         break
        }
    else {
        ind = which(dist_mat[i,-(1:i)] == min(dist_mat[i,-(1:i)]))+i;
        temp =  p2[i+1,];
        p2[i+1,] = p2[ind[1],];
        p2[ind[1],] = temp;
        }
}
#p2 <- p2[1:(i-1),]
return(p2)
}


poly_xy <- function(poly) {
  return(cbind(poly[[1]]@pts[[1]]$x,poly[[1]]@pts[[1]]$y))
}
  
  
  
Average_observations <- function(df,box_size=20) {
    cat("Averaging altimetry over boxes",sep="\n")
    breaksx <- seq(min(df$x)-1,max(df$x)+1,by=box_size)
    breaksy <- seq(min(df$y)-1,max(df$y)+1,by=box_size)
    df$box_x <- cut(df$x,breaksx)
    df$box_y <- cut(df$y,breaksy) 
    df <- ddply(df,.(box_x,box_y,t), function(x) {
      X <- x[1,]
      if (nrow(x) > 4) {
        X$z <- mean(x$z)
        X$std <- std(x$z)
      }
      return(X)
    } )
    df <- arrange(df,t)
    df$n <- 1:nrow(df)
    return(df)
  }
  
# Carried out piecewise interpolation of a function using tent functions. "first" denotes the first knot
Piecewise_interp <- function(X,ds=1,first=NULL) {

    library(plyr)
    library(Matrix)
    
    if (is.null(first)) {             # Automatically place knots
      first_val <- ceiling(min(X$x)/ds)*ds - ds
      last_val <- floor(max(X$x)/ds)*ds + ds
      knots <- seq(first_val,last_val,by=ds)  # Form notes
    } else {
      first_val <- first
      knots <- first_val
      while( tail(knots,n=1) < max(X$x)) {
        knots <- c(knots,tail(knots,n=1)+ds)
      }
      last_val <- tail(knots,n=1)
      X <- subset(X,x > first_val)
    }
    # For each point, find basis which they belong to
    X$section <- cut(X$x,knots)
    count <- 1
    X <- ddply(X,"section",function(df) {
      df$sectionnum = count
      count <<- count+1
      return(df)
    }   )

    # For each section calculate values of phi
    X <- ddply(X,"section",function(df) {
        lo_x = knots[df$sectionnum]
        hi_x = knots[df$sectionnum+1]
        c1 = 1 - (-1/ds)*lo_x
        c2 = 1 - (1/ds)*hi_x

        df$phi1 <- (-1/ds)*df$x + c1
        df$phi1_ind <- df$sectionnum
        df$phi2 <- (1/ds)*df$x + c2
        df$phi2_ind <- df$sectionnum  + 1

        return(df)
     } )

    # Build the 'observation matrix'
    C <- sparseMatrix(i = c(1:nrow(X),1:nrow(X)), j = c(X$phi1_ind, X$phi2_ind), x = c(X$phi1, X$phi2))

    # Inference with EM algorithm for observation precision
    gamma_prec <- rep(1/var(X$y),10) #Initialise precision
    maxiter=100
    
    for (m in 2:maxiter) {
        # Estimate x
        xML <- data.frame(x=knots,y=as.vector(chol2inv(chol(t(C) %*% C)) %*% t(C) %*% X$y))
        Qobs <- sparseMatrix(i=1:nrow(X),j=1:nrow(X),x=gamma_prec[m-1])
        Qtot <- t(C)%*%Qobs%*%C
        Varx <- solve(Qtot)
        
        # Update gamma_prec
        Gamma = as.vector(t(X$y) %*% X$y - 2*t(xML$y)%*% t(C) %*% X$y + sum(diag(t(C) %*% C %*% (Varx + xML$y%*%t(xML$y)))))
        gamma_prec[m] <- nrow(X)/Gamma
        if (abs(gamma_prec[m] - gamma_prec[m-1]) < 1e-5) break
    }
    cat(paste("EM converged after",m,"iterations",sep=" "),sep="\n")

    Vartrends <- NULL
    for (i in 1:(nrow(Varx)-1)) {
        Vartrends[i] <- (sum(diag(Varx)[i:(i+1)]) -2*(Varx[i,i+1]))/ds^2
    }
    trends = data.frame(year = head(knots,n=(length(knots)-1)),
         trend = diff(xML$y)/ds,
         std = sqrt(Vartrends)
         )
    return(list (x=xML, trends=trends))
}

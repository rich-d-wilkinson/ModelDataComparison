
FreqAnal2D <- function(locs,x,theta.grid=0,smoothness = 5/2,d = seq(-5,5,0.1),dither = T,plotit=F) {

  if (class(locs) == "list") {
    locs <- cbind(locs[[1]],locs[[2]])
  }

  if (class(x) == "numeric") x <- matrix(x)

   if ((dim(x)[1] > 200)&(dither==T)) {
      cat("Dithering: too much data",sep="\n")
      while(dim(x)[1] > 200) {
        x <- matrix(x[seq(1,length(x),2)])
        locs <- locs[seq(1,dim(locs)[1],2),]
        }
      }



  cat('Using a Matern field to model data',sep="\n")
  dd <- mean(diff(d))   # In km
  fs_rho = 1/dd         # In cycles per km
  d1 <- d2 <- d
  N_rho1 = length(d1)
  N_rho2 = length(d2)
  f_rho1 = fs_rho * (0 : (N_rho1-1)) / N_rho1        # Frequency axis for correlation function
  f_rho2 = fs_rho * (0 : (N_rho2-1)) / N_rho2


  D <- meshgrid(d1,d2)
  D <- cbind(c(D$x),c(D$y))
  r <- apply(D,1,function(x) {sqrt(x[1]^2 + x[2]^2) } )     # find radial distance from origin of every point

  x <- apply(x,2,function(x) {x - mean(x)})   # detrend signal       
  if (length(theta.grid) == 0) {
  fit <- MLE.Matern.fast(locs,x,smoothness=smoothness)
  
  } else {
    fit <- MLE.Matern.fast(locs,x,smoothness=smoothness,theta.grid=theta.grid)        # estimate Matern parameters
    # In MLE.Matern.fast theta is the range, rho is the sill and sigma2 is the nugget
  }
  
  rho <- matrix(Matern(abs(r),range=fit$par['theta'],nu=smoothness,phi=fit$par['rho']),N_rho2,N_rho1)    
  # In Matern range is theta and phi is the marginal variance (rho)
  RHO <- fft(rho*dd)
  PS = abs(RHO)

  prediction <- Krig(locs,x,Covariance="Matern",smoothness=smoothness,theta=fit$par['theta'],rho=fit$par['rho'],sigma2=fit$par['sigma']^2)
  out<- predict.surface(prediction)

  if (plotit) {
      image(d1,d2,t(rho),xlab='s1 (km)',ylab='s2 (km)')
      title('Matern Kernel')
    
      image(f_rho1,f_rho2,t(PS),xlab='f1 (cycles per unit)',ylab='f2 (cycles per unit)')
      title('Power Spectrum')
      
      surface( out, type="C") # option "C" our favorite
      title("Kriged Field")
  }
  
  
  fit$cutoff <- f_rho1[min(which(abs(RHO)[1,] < 0.1*abs(RHO)[1,1]))]
  fit$rho_0_1 <- as.numeric(fit$pars["theta"])*sqrt(8*smoothness)        #Lindgren
  
   out2=0
   theta_axis =  seq(1,50,0.5)
   info<- list( x=locs,y=x,smoothness=smoothness, ngrid=1)
   for(i in theta_axis) out2[i] <-  MLE.objective.fn2(log(i), info, value=T,lambda.grid=1e3)
   theta = theta_axis[which.min(out2)]
   fit$rho_0_1 <- as.numeric(theta)*sqrt(8*smoothness)        #Lindgren 
  
  
  return(fit)
}


## Find 2D Power Spectrum
PS2D <- function(s1,s2,x) {
  N1 <- length(s1)
  N2 <- length(s2)
  ds <- mean(diff(s1))
  fs <- 1/ds

  S = meshgrid(s1,s2)
  X <- fft(x*ds)
  f1 = fs * (0 : (N1-1)) / N1
  f2 = fs * (0 : (N2-1)) / N2
  f_grid <- meshgrid(f1,f2)
  PS = abs(X)^2
  return(list(f1 = f1,f2 = f2, PS = PS))
}




 MLE.objective.fn2 <- function (ltheta, info, value = TRUE,lambda.grid=100)
{
    y <- as.matrix(info$y)
    x <- info$x
    smoothness <- info$smoothness
    ngrid <- info$ngrid
    M <- ncol(y)
    Tmatrix <- fields.mkpoly(x, 2)
    qr.T <- qr(Tmatrix)
    N <- nrow(y)
    Q2 <- qr.yq2(qr.T, diag(1, N))
    ys <- t(Q2) %*% y
    N2 <- length(ys)
    theta <- exp(ltheta)
    K <- Matern(rdist(x, x)/theta, smoothness = smoothness)
    Ke <- eigen(t(Q2) %*% K %*% Q2, symmetric = TRUE)
    u2 <- t(Ke$vectors) %*% ys
    u2.MS <- c(rowMeans(u2^2))
    D2 <- Ke$values
    N2 <- length(D2)
    ngrid <- min(ngrid, N2)
    # lambda.grid <- 100
    trA <- minus.pflike <- rep(NA, ngrid)
    temp.fn <- function(llam, info) {
        lam.temp <- exp(llam)
        u2 <- info$u2.MS
        D2 <- info$D2
        N2 <- length(u2.MS)
        rho.MLE <- (sum((u2.MS)/(lam.temp + D2)))/N2
        lnDetCov <- sum(log(lam.temp + D2))
        -1 * M * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(rho.MLE) -
            (1/2) * lnDetCov)
    }
    info <- list(D2 = D2, u2 = u2.MS, M = M)
    out <- golden.section.search(f = temp.fn, f.extra = info,
        gridx = log(lambda.grid), tol = 1e-07)
    minus.LogProfileLike <- out$fmin
    lambda.MLE <- exp(out$x)
    rho.MLE <- (sum((u2.MS)/(lambda.MLE + D2)))/N2
    sigma.MLE <- sqrt(lambda.MLE * rho.MLE)
    trA <- sum(D2/(lambda.MLE + D2))
    pars <- c(rho.MLE, theta, sigma.MLE, trA)
    names(pars) <- c("rho", "theta", "sigma", "trA")
    if (value) {
        return(minus.LogProfileLike)
    }
    else {
        return(list(minus.lPlike = minus.LogProfileLike, lambda.MLE = lambda.MLE,
            pars = pars, mle.grid = out$coarse.search))
    }
}


 lscale_from_variogram  <- function(data,lim=500,plotit=T) {
      dist <- rdist(data[,1:2],data[,1:2])
      diff <- rdist(data[,3],data[,3])
      D <- data.frame(d <- c(dist),var <- c(diff)^2)
      D <- subset(D, d < lim)
      break_diff = lim/50
      breaks <- seq(0,lim*2,break_diff)
      D$group <- cut(D$d,breaks)
      D_mean <- ddply(D,"group",function(x) mean = mean(x$var))
      if (plotit == T)  {
         dev.new()
         plot(D_mean[,1],-D_mean[,2])
      }
      th <- min(D_mean[,2])+diff(range(D_mean[,2]))*0.9
      return(which(D_mean[,2]>th)[1]*break_diff)
  }

  
lscale_from_Matern  <- function(data,rho=100,nu=3/2,var=1) {
      marg <- matrix(0,length(rho),length(var))
      result_frame <- NULL
      if (ncol(data) == 3) {
        dist <- rdist(data[,1:2],data[,1:2])
      } else {
        dist <- rdist(data[,1],data[,1])
      }  
      
      for(k in 1:length(nu)) {
        nu_sel <- nu[k]
        for(i in 1:length(rho)) {
          kappa = sqrt(8*nu_sel)/rho[i]
          K_norm <- my_Matern(r=dist,nu=nu_sel,var=1,kappa=kappa)
          diag(K_norm) <- 1
          L_norm <- chol(K_norm)
          Kinv_norm <- chol2inv(L_norm) 
          marg1_norm <-  t(data$z)%*%Kinv_norm%*%data$z
          marg2_norm <-  logdet(L_norm)
          for(j in 1:length(var)) {
            #Kinv <- Kinv_norm/var[j]    
            #marg[i,j] <- -0.5*t(data$z)%*%Kinv%*%data$z - 0.5*logdet(sqrt(var[j])*L_norm)
            marg[i,j] <- -0.5*marg1_norm/var[j] - 0.5*(nrow(L_norm)*log(var[j]) + marg2_norm)
          }
        }
        rho_max_ind <- which.max(apply(marg,1,max))
        var_max_ind <- which.max(apply(marg,2,max))
        if((rho_max_ind %in% c(1,length(rho)) )| (var_max_ind %in% c(1,length(var))) ) {
          cat(paste("Warning: reached search boundary for nu = ",nu_sel,sep=""),sep="\n"); flush.console()
        }
        result_frame <- rbind(result_frame,data.frame(lscale = rho[rho_max_ind],var = var[var_max_ind], nu = nu_sel, lk = max(marg)))
      }
      
      
      return(result_frame)

    }
 
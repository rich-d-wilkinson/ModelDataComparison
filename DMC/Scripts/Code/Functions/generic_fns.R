## Find a weighted sum of normals on an axis "x", with weights "w" and mean and variance "mu" and "variance"
weighted_sum_of_normals <- function(x,w,mu,variance,method="R") {

   if (method == "R") {
     f <- rep(0,length(x))
     for (i in 1:length(w)) {
       f <- f + w[i]/sqrt(2*pi*variance[i])*exp(-(x - mu[i])^2/(2*variance[i]))
    }
   } else {
    if (.Platform$OS.type == "windows") {
      dyn.load("./Functions/gauss_mixture/R_WINDOWS/gauss_mixture.dll")
    } else {
      dyn.load("./Functions/gauss_mixture/R_UNIX/gauss_mixture.so")
    }
    X <- .C("gauss_mixture",as.double(x),as.double(w),as.double(mu),as.double(sqrt(variance)),
            as.integer(length(w)),as.integer(length(x)),f = double(length(x)))
    f <- X$f
    }
    return(f)
   }


## return the nth element from a vector (for use with lapply)
nth_element <- function(x,i)
return(x[i])


## Accept-Reject sampler
AcceptReject <- function(f,propose = 'uniform',pars = c(0.5,4),M=1)  {

 accept = F
 while (accept == F) {
    if (propose == 'uniform') {
         u <- runif(1)
         g <- runif(1)*diff(pars)+pars[1]
         if (u < f(g)/(M*1)) accept = T
    }
  }
  return(g)
 }
 
 
 # Find which rows of Tab2 are in Tab1
common_rows <- function(Tab1,Tab2) {
  nrows <- length(Tab1[,1])
  myind <-which(duplicated(rbind(Tab1, Tab2))) - nrows
  return(myind)
}

# Find difference between the p5 and p95 of the equivalent distribution of the standard devisaiton pars=pars
Findalphabeta_invgamma <- function(pars,p5,p95) {
     #return( sum(abs(1/sqrt(qgamma(c(0.05,0.95),shape=pars[1],rate=pars[2])) - c(p95,p5))))
     return( sum(abs(qinvgamma(c(0.05,0.95),shape=pars[1],scale=pars[2]) - c(p5^2,p95^2))))
  }

Findalphabeta_gamma <- function(pars,p5,p95) {
  #return( sum(abs(1/sqrt(qgamma(c(0.05,0.95),shape=pars[1],rate=pars[2])) - c(p95,p5))))
  return( sum(abs(qgamma(c(0.05,0.95),shape=pars[1],rate=pars[2]) - c(1/p95^2,1/p5^2))))
}



# Memory consumption of individual variables  
object.sizes <- function()
{
    return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name)
        object.size(get(object.name))))))
}
 
KLnorm <- function(mu1,S1,mu2,S2) {
  KL <- 0.5*(log(det(S2)/det(S1)) + sum(diag(solve(S2)%*%S1)) + t(mu1 - mu2)%*%solve(S2)%*%(mu1 - mu2) - length(mu1))
  return(KL)
  
}
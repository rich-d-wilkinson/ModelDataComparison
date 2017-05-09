

GRBF <- function(c,std,s) {
  c_ext <- matrix(c,nrow(s),2,byrow=T)
  dist_sq <- (rowSums((s - c_ext)^2))
  return(exp(-0.5* dist_sq/(std^2) ))
}

setClass("Basis",
         representation(pars="list",
                        n = "numeric",
                        fn="list"))

setClass("GRBFBasis",
         contains="Basis")

setClass("ConstBasis",
         contains="Basis")

setClass("FEBasis",
        contains="Basis")

setGeneric("basisinterp",
           function(G,s,weights)
             standardGeneric("basisinterp"))

setMethod("basisinterp",signature(G = "Basis"),  # GRBF basis with mean offset as last weight
          function(G,s,weights=NULL) {
            y <- matrix(0,nrow(s),1)
            n <- length(weights)
            if (is.null(weights)) {
              weights <- rep(1,n)
            }
            for (i in 1:n) {
              y <- y + weights[i]*(G@fn[[i]](G@pars[[i]],s))
            }
            return(as.vector(y))
          })

concatBasis <- function(G1,G2){
            this_basis <- new("Basis", pars=c(G1@pars,G2@pars), n=G1@n + G2@n, fn=c(G1@fn,G2@fn))
            return(this_basis)
          }

initGRBFbasis = function(x,y,std,nx,ny) {
  knots_x <- seq(x[1],x[2],length=(nx+2))
  knots_y <- seq(y[1],y[2],length=(ny+2))
  centres <- expand.grid(knots_x[-c(1,(nx+2))],knots_y[-c(1,(ny+2))])
  n <- nrow(centres)
  stds <- rep(std,n)
 
  fn <- pars <- list()
  for (i in 1:n) {
    fn[[i]] <-  function(pars,s) {
      return(GRBF(matrix(as.numeric(pars$centres),1,2),pars$stds,s))
    }
    pars[[i]] <- list(centres = as.matrix(centres[i,]), stds=stds[i])
  }
  this_basis <- new("GRBFBasis", pars=pars, n=nrow(centres), fn=fn)
  return(this_basis)
}

initConstbasis = function(c) {
  fn <- pars <- list()
  fn[[1]] <- function(pars,s) {
     return(pars$const) 
  }
  pars[[1]] <- list(const = c)
  this_basis <- new("ConstBasis", pars=pars, n=1, fn=fn)
  return(this_basis)
}

initFEbasis = function(p,t,M,K) {
  fn <- pars <- list()
  pars$p <- p
  pars$t <- t
  pars$M <- M
  pars$K <- K
  
  this_basis <- new("FEBasis", pars=pars, n=nrow(p), fn=fn)
  return(this_basis)
}


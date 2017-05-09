### Covariance matrices for dim(x)>1

X <- cbind(rnorm(150), rnorm(150))
X2 <- cbind(rnorm(150), rnorm(150))

create_k <- function(l=1, sigma2=1){
  # creates a covariance function with required length scale and variance.
  k <- function(r){
    # sigma2*(1+sqrt(3)*r/l)*exp(-sqrt(3)*r/l) # Matern 3/2
    sigma2 * exp(-(r/l)^2)  # squared exponential
  }
  return(k)
}

dimx=2
klist <- lapply(1:dimx,FUN=create_k)

system.time({
out <- abs(outer(X,X2, FUN="-"))
tmp<-lapply(1:dimx, FUN=function(i) klist[[i]](out[,i,,i])     )
K = tmp[[1]]*tmp[[2]]
})

system.time({
### Or by loop
K2 <- matrix(nr=dim(X)[1], nc=dim(X2)[1])
for(i in 1:dim(X)[1]){
  for(j in 1:dim(X2)[1]){
    K2[i,j] <- klist[[1]](X[i,1]-X2[j,1])*klist[[2]](X[i,2]-X2[j,2])
  }
}

})


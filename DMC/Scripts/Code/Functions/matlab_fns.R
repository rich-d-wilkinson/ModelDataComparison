##R equivalent of repmat (matlab)
repmat <- function(X,m,n){
mx = dim(X)[1]
nx = dim(X)[2]
return(matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}

spy <- function(X) {
 Y <- which(t(X)!=0,arr.ind=TRUE)
 plot(Y[,1],-Y[,2],pch=".")
}


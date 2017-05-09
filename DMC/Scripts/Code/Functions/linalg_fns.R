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


# Create a sparse diagonal matrix
sparsediag <- function(xx) {
  n <- length(xx)
  return(sparseMatrix(i=1:n,j=1:n,x=xx))
}

## Takahashi diag: Find marginal variance from Cholesky factor
# For now convert to full matrices for indexing. We need to do this intelligently in the future for matrices > 10000x10000 which would fill up memory
Takahashidiag <- function(L) {
 n <- dim(L)[1]
 X <- which(L!=0,arr.ind=T)
 i_ind <- X[,1]
 j_ind <- X[,2]
 Sigma <- Sigma + t(Sigma) - sparseMatrix(i=1:n,j=1:n,x=diag(Sigma))
 Sigma <- as(Sigma,"dgTMatrix") # coerce to i,j format
 Sigma[n,n] <- 1/L[n,n]^2
 numcomp = 0

# Sigma <- as.matrix(Sigma)
# L <- as.matrix(L)

 tic()
 for (i in seq(n-1,1,-1)) {
     Lii <- L[i,i]
     nz_indices <- intersect(i_ind[j_ind == i],(i+1):n)
     Lip1n_i<-L[nz_indices,i]
     for (j in intersect(seq(n,i,-1),which(abs(L[,i])>0))) {
          Sigma[i,j] <-  Sigma[j,i] <- (i==j)/(Lii^2) - 1/Lii*sum(Lip1n_i*Sigma[j,nz_indices ])
          numcomp = numcomp + 1
          }
     }
 toc()



 return(diag(Sigma))
}

Takahashidiag_Cseke <- function(L) {

n <- dim(L)[1]
invdiagL2 <- 1/diag(L)^2
S <- L + t(L)
S[S >0] = 1
S[n,n] =invdiagL2[n]

for (i in seq(n-1,1,-1)) {
    I   = i+which(abs(L[seq(i+1,n,1),i]) > 0)
    S[I,i] = -(S[I,I]%*%L[I,i])/L[i,i]
 	  S[i,I] = t(S[I,i])
 	  S[i,i] = invdiagL2[i] - (S[i,I]%*%L[I,i])/L[i,i];

    }
return(diag(S))
}

 ## Find the log-determinant of Q from its Cholesky L
 logdet <- function(L)  {
 diagL <- diag(L)
  return(2*sum(log(diagL)))
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


#Solve X = AQ^{-1}t(A)
cholsolveAQinvAT <- function(Q,A,Lp,P) {
  W <- t(solve(Lp,t(P)%*%t(A)))
  return(W %*% t(W))
  
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
      dyn.load("./Functions/sparseinv/R_WINDOWS/sparseinvR.dll")
    } else {
      dyn.load("./Functions/sparseinv/R_UNIX/sparseinvR.so")
    }
    X <- .C("sparseinv",as.integer(n),as.integer(Lp),as.integer(Li),as.double(Lx),as.double(d),as.integer(Up),as.integer(Uj),as.double(Ux),as.integer(Zpatp),as.integer(Zpati),result = double(znz))

    Lower <- L + sparseMatrix(i=1:n,j=1:n,x = d)
    Zp <- Lower@p
    Zi <- Lower@i
    Zx <- X$result[-which(X$result == 0)]
    Z <- sparseMatrix(p = Zpatp, i =Zpati, x = X$result,index1=F)
    return(Z)
}


#### R and C implementation of Takahashi diagonal equations by Jonty
Takahashi <- function(L, diag = TRUE, method = c("R", "C")) {

  method <- match.arg(method)
  stopifnot(inherits(L, "CsparseMatrix")) # has @i and @p


  n <- nrow(L)
  ii <- L@i + 1 # in {1,...,n}
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp) # in {1,...,n}, non-decreasing
  N = length(ii)

  stopifnot(ii >= jj,              # lower triangular
            1:n %in% ii[ii == jj]) # full diagonal


 if (method == "C") {
    tic();
    dyn.load("Test.so")
    X <- .C("TakahashiC",as.integer(n),as.integer(N),as.integer(ii),as.integer(jj),as.double(L@x),results = double(N))
    toc();
    return(X$results[ii==jj])
  } else if (method == "R") {

    if (diag) { # this to speed up the calculation
      tic()
      S <- L; S@x[] <- -999

      S@x[ii == n & jj == n] <- 1 / L[n, n]^2

      if (n > 1)
        for (i in (n-1):1) {

          k <- ii[ii > i & jj == i]      # Find row numbers with non zero indices at this column
          if (length(k) == 0) {
              S@x[ii == i & jj == i] <- 1 / L[i, i]^2
          } else {
            Lii <- L[i, i]
            Lki <- L[k, i]

            js <- rev(jj[ii %in% k & jj >= i]) # going backwards
            #for (j in js) {
            #  skj <- S@x[ii == pmax(k, j) & jj == pmin(k, j)] # select from lower triangle
            #  S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            #}
            js = unique(js)
	    js <- c(k,i)
            for(j in js) {
               skj <- apply(matrix(k),1,function(ind){ S@x[ii == max(ind,j) & jj == min(j,ind)] } )
               S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            }

          }
        }
      toc()
      return(diag(S))

    } else { # diag = FALSE : use full-size S but only fill lower triangle

      S <- matrix(NA, n, n)

      S[n, n] <- 1 / L[n, n]^2

      if (n > 1)
        for (i in (n-1):1) {

          k <- ii[ii > i & jj == i]
          if (length(k) == 0) {
              S@x[ii == i & jj == i] <- 1 / L[i, i]^2
          } else {
            Lii <- L[i, i]
            Lki <- L[k, i]

            js <- n:i # going backwards

            for (j in js) {
              skj <- S[pmax(k, j), pmin(k, j)] # select from lower triangle
              S[j, i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            }
          }
        }

      return(ifelse(is.na(S), t(S), S))
    }

  } else stop("Never get here!")
}


Zeromat <- function(ni,nj) {
  if ((ni == 0) | (nj == 0)) {
    return(NULL)
  } else {
    return(sparseMatrix(i=1, j=1,x=0,dims=c(ni,nj)))
  }
}


Imat <- function(n) {
  sparseMatrix(i=1:n, j=1:n,x=1)
}

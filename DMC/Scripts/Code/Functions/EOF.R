#Created by Pretty R at inside-R.org

#This function performs either a EOF ('Empirical Orthogonal Function' analysis) 
#on a single field (F1) or an MCA ('Maximum Covariance Analysis') on two fields
#(F1 and F2).
#The most typical setup is to have fields where columns are
#spatial locations and rows are temporal measurement. For an MCA, the row 
#dimension of F1 and F2 must be identical, but the columns can be different.

#Arguments 'F1_cols_incl' and 'F2_cols_incl' allow for the specification of
#using a subset of columns.

#In both cases, a covariance matrix is computed between columns using the 
#function 'cov4gappy()'. This function gives comparable results to: 
#'cov(F1,y=F2,use="pairwise.complete.obs")'
#whereby each covariance value is divided by n number of shared values (as opposed
#to n-1 in the case of cov(). Futhermore, the function will return a 0 (zero) in
#cases where no shared values exist between columns; the advantage being that a
#covariance matrix will still be calculated in cases of very gappy data, or when
#spatial locations have accidentally been included without observations (i.e. land
#in fields of aquatic-related parameters). 
#EOF analysis could be conducted through a direct decomposition of the 
#field F1 using SVD, although, for consistency with MCA, the the decomposition is applied
#to the covariance matrix. This is also an essential step for gappy data.

#The default settings include the centering of columns (centered=TRUE), which is usually
#recommended for EOF. Scaling (scaled=TRUE) of columns can also be set. This may
#be desired in MCA when unit ranges differ between fields. There are many different
#viewpoints on these approaches, so take caution as to the settings you choose.

#'nu' and 'nv' generally allow one to specify the return of a smaller number
#of EOF or MCA modes in order to create a smaller sized results object. 
#Additionally, two 'methods' are possible for the decomposition of the covariance
#matrix - 'svd' ('singular value decomposition') and 'irlba' ('implicitly restarted
#Lanczos bi-diagonalization') from the 'irlba' package. irlba can be a great 
#time saver in the decomposition of very large matrices, whereby the nu and nv
#settings define a stopping point. To the contrary, the svd method computes the
#full set of modes regardless, and nu and nv settings will only be used to trim
#the returned result.

#Finally, a lag may be introduced in order to explore changes in the associated
#statistics. This really only makes sense in cases where the row dimension
#represents temporal measurements taken at regular intervals.

#Main statistics returned:
#1. expl_var - 'explained variance' for each mode (not calcuable with irlba method
#as it requires the sum of the full singular value vector. 
#2. sq_cov_frac - 'squared covariance fraction' for each mode
#3. n_sig - 'number of significant' modes. This is calculated accrding to 'North's
#Rule of Thumb' which compares the difference between neighboring singular values
#(Lambda) with their associated error (Lambda_err)

#Additional information is returned which is helpful in the reconstruction of fields
#using a truncated number of modes and their associated coefficients.

eof.mca <- function(F1, F2=NULL,
                    centered=TRUE, scaled=FALSE,
                    F1_cols_incl=1:length(F1[1,]), 
                    F2_cols_incl=1:length(F2[1,]),
                    nu=NULL, nv=NULL, method="svd", F2_lag=NULL
){
  
  if(method=="irlba"){
    print("Only squared variance, not variance, is calculated using irlba method")
  }
  
  ####################
  ###eof vectors######
  ####################
  if(is.null(F2)){ # EOF
    F1 <- scale(F1, center=centered, scale=scaled)
    
    F1_center=attr(F1,"scaled:center")
    F1_scale=attr(F1,"scaled:scale")
    F2_center=NULL
    F2_scale=NULL
    
    F1_ts <- rownames(F1)
    F2_ts <- NULL
    
    F1_dim <- dim(F1)
    F2_dim <- NULL
    
    C <- cov4gappy(F1[,F1_cols_incl])
    F2_cols_incl=NULL
    
    if(method=="svd") {   
      if(is.null(nu)){
        nu=F1_dim[2]
      }
      if(is.null(nv)){
        nv=F1_dim[2]
      }
      L <- svd(C)
    }
    
    if(method=="irlba") {
      if(is.null(nu)){
        nu=5
      }
      if(is.null(nv)){
        nv=5
      }
      L <- irlba(C, nu=nu, nv=nv)
    }
    
    L$v<-NULL
  }
  
  if(!is.null(F2)){ # MCA
    F1 <- scale(F1, center=centered, scale=scaled)
    F2 <- scale(F2, center=centered, scale=scaled)
    
    F1_center=attr(F1,"scaled:center")
    F1_scale=attr(F1,"scaled:scale")
    F2_center=attr(F2,"scaled:center")
    F2_scale=attr(F2,"scaled:scale")
    
    if(!is.null(F2_lag)){
      if(sign(F2_lag)==1){
        F1 <- F1[(1+F2_lag):length(F1[,1]),]
        F2 <- F2[1:(length(F2[,1])-F2_lag),]
      }
      if(sign(F2_lag)==-1){
        F1 <- F1[1:(length(F1[,1])-F2_lag),]
        F2 <- F2[(1+F2_lag):length(F2[,1]),]
      }
    }
    
    F1_ts <- rownames(F1)
    F2_ts <- rownames(F2)
    
    F1_dim <- dim(F1)
    F2_dim <- dim(F2)
    
    C <- cov4gappy(F1[,F1_cols_incl], F2=F2[,F2_cols_incl])
    
    if(method=="svd") {
      if(is.null(nu)){
        nu=min(F1_dim[2], F2_dim[2])
      }
      if(is.null(nv)){
        nv=min(F1_dim[2], F2_dim[2])
      }
      L <- svd(C)
    }
    
    if(method=="irlba") {
      if(is.null(nu)){
        nu=5
      }
      if(is.null(nv)){
        nv=5
      }           
      L <- irlba(C, nu=nu, nv=nv)
    }
    
  }
  
  ###################
  ###eof mca stats###
  ###################
  if(method=="svd"){
    expl_var=L$d/sum(L$d) #explained variance
    sq_cov_frac=L$d^2/sum(L$d^2) #squared covariance fraction
  }
  if(method=="irlba"){
    expl_var=NULL
    if(dim(C)[1] <= dim(C)[2]){
      sq_cov_frac=L$d^2/sum(diag(C%*%t(C)))
    } else {
      sq_cov_frac=L$d^2/sum(diag(t(C)%*%C))
    }
  }
  
  if(length(L$d)>1){
    Lambda_err <- sqrt(2/min(F1_dim[2], F2_dim[2]))*L$d
    upper.lim <- L$d+Lambda_err
    lower.lim <- L$d-Lambda_err
    NORTHok=0*L$d
    for(i in seq(L$d)){
      Lambdas <- L$d
      Lambdas[i] <- NaN
      nearest <- which.min(abs(L$d[i]-Lambdas))
      if(nearest > i){
        if(lower.lim[i] > upper.lim[nearest]) NORTHok[i] <- 1
      }
      if(nearest < i){
        if(upper.lim[i] < lower.lim[nearest]) NORTHok[i] <- 1
      }
    }
    n_sig <- min(which(NORTHok==0))-1
  }
  
  ##########################################################
  ###expansion of eof coefficients "principle components"###
  ##########################################################
  
  A_coeff = NULL
  A_norm = NULL
  A = NULL
  B_coeff = NULL
  B_norm = NULL
  B = NULL
  
  #trim columns of original data
  F1 <- as.matrix(F1[,F1_cols_incl])
  
  #setup for norm
  F1_val<-replace(F1, which(!is.na(F1)), 1)
  F1_val<-replace(F1_val, which(is.na(F1_val)), 0)
  
  #calc of expansion coefficient and scaling norm
  A_coeff <- replace(F1, which(is.na(F1)), 0)%*%L$u[,1:nu]
  A_norm <- F1_val%*%(L$u[,1:nu]^2)
  A=A_coeff/A_norm
  
  if(!is.null(F2)){
    #trim columns of original data then center then scale
    F2 <- F2[,F2_cols_incl]       
    
    #setup for norm
    F2_val<-replace(F2, which(!is.na(F2)), 1)
    F2_val<-replace(F2_val, which(is.na(F2_val)), 0)
    
    #calc of expansion coefficient and scaling norm
    B_coeff <- replace(F2, which(is.na(F2)), 0)%*%L$v[,1:nv]
    B_norm <- F2_val%*%(L$v[,1:nv]^2)
    B=B_coeff/B_norm
  }
  
  ############
  ###result###
  ############
  RESULT <- list(
    u=L$u[,1:nu], v=L$v[,1:nv], 
    Lambda=L$d, Lambda_err=Lambda_err,
    expl_var=expl_var, sq_cov_frac=sq_cov_frac, 
    NORTHok=NORTHok, n_sig=n_sig, nu=nu, nv=nv,
    F1_dim=F1_dim, F2_dim=F2_dim,
    F1_cols_incl=F1_cols_incl,  F2_cols_incl=F2_cols_incl,
    F1_center=F1_center, F1_scale=F1_scale, 
    F2_center=F2_center, F2_scale=F2_scale,
    A=A, B=B,
    F1_ts=F1_ts, F2_ts=F2_ts, F2_lag=F2_lag
  )
  
  RESULT
} 
cov4gappy <- function(F1, F2=NULL){
  if(is.null(F2)){
    F1 <- as.matrix(F1)
    F1_val<-replace(F1, which(!is.na(F1)), 1)
    F1_val<-replace(F1_val, which(is.na(F1_val)), 0) 
    n_pairs=(t(F1_val)%*%F1_val)
    
    F1<-replace(F1, which(is.na(F1)), 0)
    cov_mat <- (t(F1)%*%F1)/n_pairs
    cov_mat <- replace(cov_mat, which(is.na(cov_mat)), 0)
  }
  
  if(!is.null(F2)){
    if(dim(F1)[1] == dim(F2)[1]){
      F1 <- as.matrix(F1)
      F2 <- as.matrix(F2)
      
      F1_val<-replace(F1, which(!is.na(F1)), 1)
      F1_val<-replace(F1_val, which(is.na(F1_val)), 0) 
      F2_val<-replace(F2, which(!is.na(F2)), 1)
      F2_val<-replace(F2_val, which(is.na(F2_val)), 0) 
      n_pairs=(t(F1_val)%*%F2_val)
      
      F1<-replace(F1, which(is.na(F1)), 0)
      F2<-replace(F2, which(is.na(F2)), 0)
      cov_mat <- (t(F1)%*%F2)/n_pairs
      cov_mat <- replace(cov_mat, which(is.na(cov_mat)), 0)
      
    } else {
      print("ERROR; matrices columns not of the same lengths")
    }
  }
  
  cov_mat
}
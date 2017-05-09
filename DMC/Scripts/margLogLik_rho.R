### Plot the marginal likelihood over rho:
#
#
# 1. Fix other parameters:
fix.param <- c(0.5,1,0.01)
rho.seq <- seq(100,1500,10)
n.rho <- length(rho.seq)
n.res <- length(residual_Obs)
negLogLik <- array(0, dim=c(n.rho,1))
#
# 2. calculate log-likelihood at the parameters
for (i in 1:n.rho){
  params <- c(fix.param[1], rho.seq[i],fix.param[2:3])
  negLogLik[i] <- covariance.neg.lik.calc(params, residual_Obs,n.res,meas.error, dist.mat, "Matern")
}

plot(rho.seq, negLogLik)


#####
# using the maaxmised parameters to check if it does give the minimum:
fix.param <- mle.output.matern.nug$par 
rho.seq <- seq(1000,3500,10)
n.rho <- length(rho.seq)
n.res <- length(residual_Obs)
negLogLik <- array(0, dim=c(n.rho,1))
#
# 2. calculate log-likelihood at the parameters
for (i in 1:n.rho){
  params <- c(fix.param[1], rho.seq[i],fix.param[3:4]) #fix.param[2] is the MLE of rho
  negLogLik[i] <- covariance.neg.lik.calc(params, residual_Obs,n.res,meas.error, dist.mat, "Matern")
}

plot(rho.seq, negLogLik)
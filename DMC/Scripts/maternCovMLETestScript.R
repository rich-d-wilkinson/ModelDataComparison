## R script to test the matern covariance function parameter estimation.
#We generate data using a provided variogram function and then try and recover the set parameter values. 
#
# 1. Choose variogram model:
#
rho <- 35
sigma <- 50
nu <- 0.5
matern.Exp <- c(0.5,rho,sigma) # Exponential variogram model (valid on the sphere) (psill = sigma and range = rho)

# 2. Simulate data on the sphere:
#
n <- 80
lon_ <- runif(n,-25,-20)
lat_ <- runif(n,-45,-40)


#sph.coord <- SpatialPoints(cbind(lon_,lat_), CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
#model.choice <- vgm(sigma,"Exp",rho) #exponential covariance function (valid on the sphere)
#g.dummy <- gstat(formula = val~1, dummy = TRUE, beta = 0,
#                 model = model.choice, nmax = 50) # for speed -- 10 is too small!!
#The predict code uses is.projected to check for lon,lat coordinates (I think in degrees)
#val <- predict(g.dummy, sph.coord, nsim = 1)

# simulate data from GP:
dist.mat <- geoDistance(cbind(lon_,lat_), cbind(lon_,lat_))
cov.param <- c(nu, rho, sigma)
matern.Cov <- maternCov(dist.mat, cov.param)
data <- mvrnorm(1, rep(0,n), matern.Cov)

# 3. Recover parameters:
#
ini.true <- c(rho,sigma)
ini.off <- c(45,55)
param.lower <- c(0,0) #all parameters are strictly positive (should we have c(0,0) then?

#check we don't start at non-positive definite matrix:
cov.param <- c(nu, ini.true)
Cov <- maternCov(dist.mat, cov.param)
det(Cov)
cov.param <- c(nu, ini.off)
Cov <- maternCov(dist.mat, cov.param)
det(Cov)
##

#param.upper <- c(inf,inf,inf)
mle.output <- optim(ini.true, sub.matern.neg.lik.calc.constrain, nu = 0.5, data = data, n = n, dist.mat = dist.mat, method="BFGS")
mle.output.NM <- optim(ini.true, sub.matern.neg.lik.calc.constrain, nu = 0.5, data = data, n = n, dist.mat = dist.mat, method="Nelder-Mead")
mle.output.NM
mle.output

mle.output.off <- optim(ini.off, sub.matern.neg.lik.calc.constrain, nu = 0.5, data = data, n = n, dist.mat = dist.mat, method="BFGS") #, control = list(trace=4))
mle.output.off.NM <- optim(ini.off, sub.matern.neg.lik.calc.constrain, nu = 0.5, data = data, n = n, dist.mat = dist.mat, method="Nelder-Mead") #, control = list(trace=4))
mle.output.off.NM
mle.output.off

mle.output.off.2 <- optim(mle.output.off.NM$par + c(0.01,0.01), sub.matern.neg.lik.calc, nu = 0.5, data = data, n = n, dist.mat = dist.mat, method="BFGS") #, control = list(trace=4))
mle.output.off.2

#constrained lik:
mle.output.uncon <- optim(ini.true, sub.matern.neg.lik.calc, nu = 0.5, data = data, n = n, dist.mat = dist.mat, lower = param.lower, method="L-BFGS-B")
mle.output.uncon

#Estimating nu as well:
ini.true <- c(nu, rho, sigma)
ini.off <- c(0.3, rho, sigma)

cov.param <- ini.true
Cov <- maternCov(dist.mat, cov.param)
det(Cov)
cov.param <- ini.off
Cov <- maternCov(dist.mat, cov.param)
det(Cov)

param.lower <- c(0.01,0.01,0.01)
#needs upper limit for nu <= 1/2 too.
param.upper <- c(1/2,10000,100000)
mle.output.uncon <- optim(ini.true, matern.neg.lik.calc, data = data, n = n, dist.mat = dist.mat, lower = param.lower, upper = param.upper, method="L-BFGS-B",control = list(trace=5))
mle.output.uncon


####################
# Why does rho not change?
#This was because the cov matrix was too close to a diagonal because rho and sigma were too small. 

rho_vec <- seq(0.2,0.6,by = 0.01)
n1 <- length(rho_vec)
lik <- array(0, dim = c(1,n1))
det.cov <- array(0, dim=c(1,n1))
for(i in 1:n1){
  cov.param <- c(nu, rho_vec[i],sigma)
  Cov <- maternCov(dist.mat, cov.param)
  det.cov[i] <- det(Cov)
  lik[i] <-  1/2*log(det(Cov)) + 1/2*t(data)%*%solve(Cov)%*%data
}






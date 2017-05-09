dyn.load("Test.so")

hello2 <- function(n) {
.C("hello", as.integer(n))
}

ksmooth <- function(x, xpts, h) {
n <- length(x)
D <- outer(x, xpts, "-")
K <- dnorm(D / h)
dens <- colSums(K) / (h * n)
}

ksmoothC <- function(x, xpts, h) {
n <- length(x)
nxpts <- length(xpts)
dens <- .C("kernel_smooth", as.double(x), as.integer(n),
as.double(xpts), as.integer(nxpts), as.double(h),
result = double(length(xpts)))
dens[['result']]
}

library(dpca)


## random spectrum
set.seed(123)
dim <- 500L
m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)

spec <- m %*% t(Conj(m))
spec[1:5, 1:5]


res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 3L)

edec <- eigen(spec)

res$values
edec$values[1:3]

res$vectors[1:10,1:3]
edec$vectors[1:10,1:3]


plot(x = Re(res$vectors[,1]), y = Re(edec$vectors[,1]))


## spectrum of x_it = u_{t-i} at pi/4
spec <- diag(1 + 0i, nrow = dim, ncol = dim)
for (i in 1:dim) {
  for (j in 1:dim) {
    spec[i,j] <- exp((i-j) * 1i * pi/4)
  }
}
spec[1:5, 1:5]

res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 3L)

edec <- eigen(spec)

res$values
edec$values[1:3]

res$vectors[1:4,1:3]
edec$vectors[1:4,1:3]

unique(round(Re(res$vectors[,1]), 4))
unique(round(Re(edec$vectors[,1]), 4))

unique(round(log(edec$vectors[,1] / res$vectors[,1]), 4))




## spectrum of x_it = u_{t-i} at 0
spec <- matrix(1 + 0i, nrow = dim, ncol = dim)
res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 2L, 1e-16)

edec <- eigen(spec)

res$values
res$vectors[1:10,1]
edec$values[1:3]

res$vectors[1:4,1]
edec$vectors[1:4,1:3]

plot(x = Re(res$vectors[,1]), y = Re(edec$vectors[,1]))
plot(x = Re(res$vectors[,2]), y = Re(edec$vectors[,2]))


## spectrum of x_it = u_{t-i} at 0
set.seed(123)
dim <- 10L
k <- 5
m <- matrix(cos(runif(dim*k, -pi, pi)) + 1i * sin(runif(dim*k, -pi, pi)), ncol = dim)

spec <- t(m) %*% Conj(m)
res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 1L, 1e-26)

edec <- eigen(spec)

res$values
res$vectors[1:10,1]
res$avectors[1:10,1]
##res$vectors[1:10,3]
edec$values[1:3]

res$vectors[1:4,1]
edec$vectors[1:4,1]

round(log(edec$vectors[,1] / res$vectors[,1]), 4)
round(edec$vectors[,1] / res$vectors[,1], 4)

plot(x = Re(res$vectors[,1]), y = Re(edec$vectors[,1]))
plot(x = Re(res$vectors[,2]), y = Re(edec$vectors[,2]))


spec %*% res$vectors[,1] - res$values[1] * res$vectors[,1]
spec %*% edec$vectors[,1] - edec$values[1] * edec$vectors[,1]

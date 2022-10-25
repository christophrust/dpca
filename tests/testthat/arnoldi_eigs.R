library(dpca)


## random spectrum
test_that("Eigenvalue and -vector equality", {
  dim <- 150L
  m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
  spec <- m %*% t(Conj(m))

  r1 <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 4L, .Machine$double.eps, 1L, 0L)
  r2 <- eigen(spec)

  diff_cp <- abs(Re(t(Conj(r1$vectors)) %*% r2$vectors[,1:4])) - diag(4)

  expect_lte(sum(Re(r1$values - r2$values[1:4])^2 ), 1e-12)
  expect_lte(sum(diff_cp^2 ), 1e-12)

})








c_vec <- runif(dim) + runif(dim) * 1i
system.time(a <- .Call("R_zMatVec", spec, c_vec, as.integer(length(c_vec)), 1L))
system.time(b <- .Call("R_zMatVec", spec, c_vec, as.integer(length(c_vec)), 2L))
spec %*% c_vec
a-b
plot(Re(a)~Re(b))
sum(abs(a-b))


system.time(res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 4L, .Machine$double.eps))

system.time(edec <- eigen(spec))

res$values
edec$values[1:4]

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
res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 1L, 1e-16)

edec <- eigen(spec)

res$values
res$vectors[1:10,1]
res$vectors[1:10,3]
edec$values[1:3]

res$vectors[1:4,1]
edec$vectors[1:4,1]
round(log(edec$vectors[,1] / res$vectors[,1]), 4)
round(edec$vectors[,1] / res$vectors[,1], 4)

plot(x = Re(res$vectors[,1]), y = Re(edec$vectors[,1]))
plot(x = Re(res$vectors[,2]), y = Re(edec$vectors[,2]))


spec %*% res$vectors[,1] - res$values[1] * res$vectors[,1]
spec %*% edec$vectors[,1] - edec$values[1] * edec$vectors[,1]


spec <- matrix(1, 5,5) + 1i
res <- .Call("R_arnoldi_eigs", mat = spec, dim = 5L, 1L, 1e-16)




## spectrum of x_it = u_{t-i} at 0
set.seed(123)
dim <- 10L
k <- 5
m <- matrix(cos(runif(dim*k, -pi, pi)) + 1i * sin(runif(dim*k, -pi, pi)), ncol = dim)

spec <- round(t(m) %*% Conj(m), digits = 4)

res <- .Call("R_arnoldi_eigs", mat = spec, dim = dim, 1L, .Machine$double.eps)

spec %*% res$vectors[,1] - res$values[1] * res$vectors[,1]

paste0("m = [", paste0(vapply(1:nrow(spec), function(x){
  paste0(spec[x,], collapse = ",")
}, ""), collapse= ";"), "];")
##  [V,D] = eigs(m, 1)





res$vectors/ c(-0.001307 + 0.032186i,
  -0.089224 - 0.197479i,
  -0.440288 + 0.321553i,
  -0.164862 - 0.303179i,
   0.010640 + 0.033460i,
  -0.127271 - 0.099264i,
  -0.118197 - 0.381526i,
   0.136564 - 0.232786i,
  -0.285395 - 0.338237i,
  -0.206347 - 0.193822i)

spec %*% c(0.715292-0.909374i
, 0.250232-0.983513i
, 0.800279-0.551713i
, 0.613106-0.150973i
, -0.191170+0.920083i
, 0.377532+0.864621i
, -0.940258-0.202890i
, 0.022660+0.377654i
, 0.550207-0.990628i
, -0.345216+0.577269i)

v <- c(0.715292-0.909374i
, 0.250232-0.983513i
, 0.800279-0.551713i
, 0.613106-0.150973i
, -0.191170+0.920083i
, 0.377532+0.864621i
, -0.940258-0.202890i
, 0.022660+0.377654i
, 0.550207-0.990628i
, -0.345216+0.577269i)
.Call("R_zMatVec", spec, v, as.integer(length(v)))


spec %*% c(0.040282+0.161256i
, 0.814374+0.437483i
, 0.780500+0.919456i
, 0.676322+0.031354i
, 0.395854+0.130229i
, 0.145070+0.713476i
, 0.651700+0.363842i
, 0.844054+0.733210i
, 0.014318+0.421981i
, 0.108302+0.530101i)



mm <- matrix(runif(9) + runif(9)*1i, ncol = 3)
v <- runif(3) + runif(3)* 1i

mm
v
mm %*% v
mm[1,] * v
cumsum(mm[1,] * v)
.Call("R_zMatVec", mm, v, as.integer(length(v)), 1L)
.Call("R_zMatVec", mm, v, as.integer(length(v)), 2L)
mm %*% v

dim <- 150L
m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
spec <- m %*% t(Conj(m))

cat("Integer size:", .Machine$sizeof.long, "\n")
cat("sizeof(a_int) in C should match\n")

cat("calling eigen\n")
r2 <- eigen(spec)

cat("called eigen\n")
gc()

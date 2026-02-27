# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)

dim <- 150L
m <- matrix(cos(runif(dim^2, -pi, pi)) + 1i * sin(runif(dim^2, -pi, pi)), ncol = dim)
spec <- m %*% t(Conj(m))

cat("TESTTEST:", .Machine$sizeof.long, "\n")
cat("Integer size:", .Machine$sizeof.long, "\n")
cat("sizeof(a_int) in C should match\n")

cat("calling eigen\n")
r2 <- eigen(spec)

cat("called eigen\n")
gc()

cat("TESTTEST finished:", .Machine$sizeof.long, "\n")

library(dpca)

test_check("dpca")

# Dynamic Principal Component Analysis and Generalized Dynamic Factor Model Estimation

<!-- badges: start -->
[![R-CMD-check](https://github.com/christophrust/dpca/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/christophrust/dpca/actions/workflows/check-standard.yaml)
[![License](https://img.shields.io/github/license/christophrust/dpca)](./LICENSE)
<!-- badges: end -->

This package provides a fast and reliable implementation of dynamic principal component analysis ala Brillinger. The main use case is estimation of (real) dynamic factor models (ala Forni, Lippi, Reichlin, Hallin,...).

We are aware of the R package `freqdom`, developed by Siegried Hörmann and Lukas Kidzinsiki, which is a pure R implementation and therefore quite slow and not really useful for the use in simulations. We draw on many implementation designs used by them. This package to some extent is a `C` port of `freqdom` but deviates at some points. For instance, the convoluted filter which computes the dynamic common component from the output in freqdom is obtained by filtering the output twice: first to get the inputs \(what in `freqdom` is called "scores" in analogy to their FDA context\) and, second these inputs are filtered again to get the dynamic common component \("KLexpansion"\). `dpca` computes the convolution in the frequency domain. The advantage of this approach is that this filter is invariant with respect to multiplications of dynamic eigenvectors by a unit-length complex number.

This package currently is not in a very mature state and will hopefully soon be available on CRAN in a more mature version.

## Things implemented so far:

- [x] Estimation of multivariate spectral density using lag-window technique.
- [x] Fast computation of dynamic eigenvalues/eigenvectors of spectrum using [`ARPACK`](https://en.wikipedia.org/wiki/ARPACK).
- [x] Discrete fourier transformation to obtain filters/transfer functions.
- [x] The selection criterion of [`Hallin & Liska (2007)`](https://doi.org/10.1198/016214506000001275) to determine the number of dynamic factors.
- [ ] One-Sided representation of the the dynamic common component using the approach of [`Forni, Hallin, Lippi, Zaffaroni (2015)`](http://dx.doi.org/10.1016/j.jeconom.2013.10.017).
- [ ] Forecasting methods...
- [ ] Model assessment...
- [ ] ...
- [x] Ship ARPACK



## Installation

`dpca` requires [`ARPACK`](https://en.wikipedia.org/wiki/ARPACK) to compute dynamic eigenvalues/eigenvectors using the Implicitly Restarted Arnoldi Method which is much faster than `base::eigen()` which does full spectral decomposition which is in most times not required. `ARPACK` is shipped with `dpca`.

`dpca`can be installed with `devtools`:

```r
devtools::install_github("https://github.com/christophrust/dpca.git")
```


## Examples

```r
set.seed(123456)
## simulate some data
nrx <- 100L
ncx <- 1000L
q <- 4L

epsilon <- matrix(rnorm(ncx * q), nrow = q)

b_filter <- vapply(1:10, function(l) {
  matrix(rnorm(q * nrx, sd = 1/l), q, nrx)
}, matrix(0, q, nrx))

chi <- .Call("R_filter_process", b_filter, epsilon, as.integer(1:10),
             nrx, q, q, ncx, 10L, 1L, 0L, 0L, PACKAGE = "dpca")

x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))
bw <- as.integer(floor(ncol(x)^(1/3)))



dpc <- dpca(x, q = q, bandwidth = bw)
str(dpc)
```

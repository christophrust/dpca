# Dynamic Principal Component Analysis and Generalized Dynamic Factor Model Estimation

<!-- badges: start -->
[![R-CMD-check](https://github.com/christophrust/dpca/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/christophrust/dpca/actions/workflows/check-standard.yaml)
[![Codecov test coverage](https://codecov.io/gh/christophrust/dpca/branch/main/graph/badge.svg)](https://app.codecov.io/gh/christophrust/dpca?branch=main)
[![License](https://img.shields.io/github/license/christophrust/dpca)](./LICENSE)
<!-- badges: end -->

This package provides a fast and reliable implementation of Brillinger's dynamic principal component analysis. The main use case is estimation of Generalized Dynamic Factor Models ([Forni & Lippi, 2001](https://doi.org/10.1017/S0266466601176048), [Forni, Hallin, Lippi and Reichlin, 2000](https://www.jstor.org/stable/2646650)).

At its core, this package features multivariate spectral density estimation based on the lag-window estimator. Based on an estimated spectrum, dynamic principal components are computed using an efficient numerical procedure for eigenvalue/eigenvector decomposition (Implicitly Restarted Arnoldi Method).

Moreover, `dpca` implements the method to select the number of dynamic principal components of [Hallin & Liska (2007)](https://doi.org/10.1198/016214506000001275). Also the information criteria by [Bai & Ng, 2003](https://doi.org/10.1111/1468-0262.00273) are implemented, with the additional feature that the constant multplier in the penalty term is chosen in a data-driven way analogously to [Hallin & Liska (2007)](https://doi.org/10.1198/016214506000001275).

## Things implemented so far:

- [x] Estimation of multivariate spectral density using lag-window technique.
- [x] Fast computation of dynamic eigenvalues/eigenvectors of spectrum using [`ARPACK`](https://en.wikipedia.org/wiki/ARPACK).
- [x] Discrete fourier transformation to obtain filters/transfer functions.
- [x] The selection criterion of [Hallin & Liska (2007)](https://doi.org/10.1198/016214506000001275) to determine the number of dynamic factors.
- [x] Ship ARPACK
- [x] Refactor C code
- [x] Interface to common time series data formats (`zoo`, `ts`).
- [ ] One-Sided representation of the the dynamic common component using the approach of [`Forni, Hallin, Lippi, Zaffaroni (2015)`](http://dx.doi.org/10.1016/j.jeconom.2013.10.017).
- [ ] Forecasting methods...
- [ ] Model assessment...
- [ ] ...

We are aware of the R package [`freqdom`](https://cran.r-project.org/web/packages/freqdom/index.html), developed by Siegfried Hörmann and Lukas Kidzinsiki which is a pure `R` implementation. `dpca` is written mainly in `C`. Although providing a similiar interface to that of `freqdom`, `dpca` has some unique features apart from being much faster.

For instance, the convoluted filter which computes the dynamic common component from the output in `freqdom` is obtained by filtering the output twice: first to get the inputs \(what in `freqdom` is called "scores" in analogy to their FDA context\) and, second these inputs are filtered again to get the dynamic common component \("KLexpansion"\). `dpca` computes the convolution in the frequency domain. The advantage of this approach is that this filter is invariant with respect to multiplications of dynamic eigenvectors by a unit-length complex number.

## Installation

`dpca` depends on [`ARPACK`](https://en.wikipedia.org/wiki/ARPACK) to compute dynamic eigenvalues/eigenvectors using the Implicitly Restarted Arnoldi Method which is much faster than R's `base::eigen()` whenever a truncated instead of the full spectral decomposition is required. `ARPACK` is shipped with `dpca`.

To install `dpca` using `devtools`:

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

epsilon <- matrix(rnorm((ncx + 10) * q), nrow = q)

b_filter <- vapply(1:10, function(l) {
  matrix(rnorm(q * nrx, sd = 1/l), nrx, q)
}, matrix(0, nrx, q))

chi1 <- multivariate_filter(epsilon, b_filter, as.integer(1:10))

x <- chi + rnorm(nrx * ncx, sd = 0.1 * sd(chi))
bw <- as.integer(floor(ncol(x)^(1/3)))

dpc <- dpca(x, q = q, bandwidth = bw)
str(dpc)
```

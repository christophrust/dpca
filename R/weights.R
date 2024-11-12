weights_bartlett <- function(x) {
  1 - abs(x)
}
weights_trunc <- function(x) {
  x[] <- 1
  x
}
weights_tukey <- function(x) {
  0.54 + 0.46 * cos(pi * x)
}
weights_parzen <- function(x) {
  res <- x
  yes <- abs(x) <= 0.5
  res[yes] <- 1 - 6 * x[yes]^2 + 6 * abs(x[yes])^3
  res[!yes] <- 2 * (1 - abs(x[!yes]))^3
  res
}
weights_bohman <- function(x) {
  ((1 - abs(x))^3) * cos(pi * x) + sin(pi * abs(x)) / pi
}
weights_daniell <- function(x) {
  res <- sin(pi * x) / (pi * x)
  res[res > 1] <- 1
  res[is.nan(res)] <- 1
  res
}

weights_barlett_priestley <- function(x) {
  min(abs(3 * (min(sin(pi * x) / pi * x, 1) - cos(pi * x)) / (pi * pi * x * x)), 1)
}

weights_parzen_cogburn_davis <- function(x, r = 1) {
  1 / (1 + abs(x)^(2 * r))
}

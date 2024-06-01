##
c_ts <- function(x, info_ts){             # convert ts
  start = extract_start_date(info_ts)
  return(ts(x, start = start, freq = attributes(info_ts)$tsp[3]))
}

extract_start_date <- function(x){

 freq <- attributes(x)$tsp[3]
 start_date <- time(x)[1]
 year <- trunc(start_date)

 if(freq == 4){
  rest_qu <- c(0, 0.25, 0.5, 0.75)
  Q <- match((start_date - year), rest_qu)
  return(c("Year" = year, "Quarter" = Q))
 }

 if(freq == 12){
  rest_md <- c(0, sapply(1:11, function(i) i / 12))
  M <- which.min(abs(rest_md - (start_date - year)))
  return(c("Year" = year, "Month" = M))
 }

}

dta_md_l <- read.csv("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv")     # l for level

start_md <- c(1959, 1)    # 1st month 1959

tcode_md <- dta_md_l[1, - 1]     # tranformation code

dta_md_l <- ts(as.matrix(dta_md_l[- 1, - 1]), start = start_md, freq = 12)

## define a function in order to transform the data to stationary data
## according McCracken and Ng (2015) (appendix, p.26)

transform <- function(xx, tcode) {
  if(tcode == 1) return(xx)
  if(tcode == 2) return(c(NA, diff(xx)))
  if(tcode == 3) return(c(NA, NA, diff(xx, differences = 2)))
  if(tcode == 4) return(log(xx))
  if(tcode == 5) return(c(NA, diff(log(xx))))
  if(tcode == 6) return(c(NA, NA, diff(log(xx), differences = 2)))
  if(tcode == 7) return(c(NA, NA, diff(xx / lag(xx, k = - 1) - 1)))
}


N_md <- ncol(dta_md_l)

fredmd <- c_ts(
  sapply(
    seq_len(N_md),
    function(i) transform(dta_md_l[, i], tcode_md[[i]])
  ),
  dta_md_l
)
colnames(fredmd) <- colnames(dta_md_l)

## remove missings
fmd <- fredmd[, colSums(is.na(fredmd)) < 10]
which(rowSums(is.na(fmd)) > 0)

fmd <- window(fmd, start = start(fmd) + c(0, 2), end = end(fmd) - c(0, 2))
which(rowSums(is.na(fmd)) > 0)
which(colSums(is.na(fmd)) > 0)
fmd <- fmd[, colSums(is.na(fmd)) < 1]
any(is.na(fmd))
fredmd <- fmd

save(fredmd, file = "data/fredmd.rda")

data(fredmd)

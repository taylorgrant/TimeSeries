
arfimaSim <- function(n, d, p=NULL, q=NULL) {
  require(signal)
  # creates a fractionally integrated (p,d,q) process of length n (extends to (1,d,1))
  # d: fractional differencing parameter (e.g. d=.4, produces series that FI(.4)
  # p: value of AR parameter
  # q: value of MA
  #
  x <- rnorm(n)
  d <- -d
  n <- length(x)
  k <- seq(1:(n-1))
  b <- (k-d-1)/k 
  b <- c(1,cumprod(b))
  if (is.null(p) && is.null(q)) {
    df <- filter(b,1,x)
  } else if (!is.null(p) && is.null(q)) {
    df <- filter(b,1,x)
    arx <- filter(1,c(1,-p), df)
  } else if (is.null(p) && !is.null(q)) {
    df <- filter(b,1,x)
    max <- filter(c(1,q), 1, df)
  } else if (!is.null(p) && !is.null(q)) {
    df <- filter(b, 1, x)
    marx <- filter(c(1,q),c(1,-p), df)
  }
}

############ Examples ############
### ARFIMA (0,d,0) 
# df <- arfimaSim(500, .75)
#
### ARFIMA (1,d,0) 
# ar.df <- arfimaSim(500, .75, .4)
#
### ARFIMA (0,d,1)
# ma.df <- arfimaSim(500, .75, 0, .4)
# 
### ARFIMA (1,d,1)
# arma.df <- arfimaSim(500, .75, .4, .4)

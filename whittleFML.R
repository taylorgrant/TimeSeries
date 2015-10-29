# ------------------------------------------ #
# Frequency Domain Whittle Maximum Likelihood 
# Estimator of Fractional Integration
# 
# Main Function: whittleFML
# Requires: specDens; whitFunc

# ------------------------------------------ #
# Estimating the spectral density 
# ------------------------------------------ #
specDens <- function(params, p,q, n) {
  # At Fourier frequencies 2*pi*j/m (j = floor(m/2))
  # n = length of time series
  
  j <- (n-1) %/% 2
  d <- params[1]
  
  ## Estimate the Fourier frequencies
  fs <- 2*pi*(1:j)/n
  
  ## Calculate the spectral density function
  if(p > 0) {
    phi <- params[2:(p+1)]
    px <- outer(fs, 1:p)
    Rar <- cos(px) %*% phi # matrix product
    Iar <- sin(px) %*% phi
    
    far <- (1-Rar)^2 + Iar^2
  } else {
    phi <- numeric()
    far <- 1
  }
  
  if(q > 0) {
    theta <- params[(p+2):(p+q+1)]
    px <- outer(fs, 1:q)
    Rma <- cos(px) %*% theta
    Ima <- sin(px) %*% theta
    
    fma <- (1+Rma)^2 + Ima^2
  } else {
    theta <- numeric()
    fma <- 1
  }
  
  ## Spectral Density: 
  sif <- (2 - 2*cos(fs))^(-d)*fma/far  
  
  r <- list(fs = fs, sif = sif,
            pq = c(p,q), params = c(d = d, phi = phi, theta = theta))
}

# ------------------------------------------ #
# whitFunc - function to be minimized by MLE
# ------------------------------------------ #
whitFunc <- function(params, n, In, pq=pq, give.w.only = FALSE) {
  
  p <- pq[1]
  q <- pq[2]
  
  spD <- specDens(params,p,q,n)
  sif <- spD$sif
  fy <- In/sif
  
  ## w: criterion function from Hauser (1999)
  K <- length(sif)
  w <- K*(log(2*pi) -1 + log(sum(fy)/K)) + sum(log(sif)) 
  
  if(give.w.only)
    return(w)
  
  list(n = n, d = params[1], params = params, w = w, sif = sif)
}


# ------------------------------------------ #
# Frequency Domain MLE - Whittle Estimator #
# ------------------------------------------ #
whittleFML <- function(x, p, q, n=length(x), 
                       inits = list(d=0, AR = numeric(), MA = numeric())) {
  # check initial call
  if(missing(p))  p <- length(inits$AR)
  else {
    stopifnot(length(inits$AR) == p)
    if(0 > (p <- as.integer(p))) stop("must have integer p >= 0")
  }
  if(missing(q))  q <- length(inits$MA)
  else {
    stopifnot(length(inits$MA) == q)
    if(0 > (q <- as.integer(q))) stop("must have integer q >= 0")
  }
  
  ## lining up parameters 
  pq <- c(p,q)
  params <- unlist(inits) # c(H, ar[1:p], ma[1:q]) # where p=0, q=0 is possible
  d0 <- params[1]
  n.params <- length(params)
  
  ## generate periodogram using fast fourier transform, drop 1st
  In <- (Mod(fft(x))^2/(2*pi*n)) [2:((n+1) %/% 2)]
  
  ## define function to be minimized by ML estimator 
  FMLwrap <- function(params) {
    r <- whitFunc(params, n = n, In = In, pq = pq, give.w.only = TRUE)
  }
  
  ## MLE method depends on number of parameters to estimate
  if(n.params == 1) { # one dimensional - ARFIMA (0,d,0) 
    result <- optim(par = params, fn = FMLwrap, lower = -.99, upper = 0.99,
                    method = "Brent", hessian = TRUE)
    params <- c(d = result$par)
  } 
  else { # ARFIMA (p,d,q)
    result <- optim(par = params, fn = FMLwrap, lower = -.99, upper = .99, 
                    method = "L-BFGS-B", hessian = TRUE)
    params <- result$par
  }
  
  ## Estimated Results
  whittle <- whitFunc(params, n, In, pq)
  
  Vcov <- solve(result$hessian)
  ssd <- sqrt(diag(Vcov))
  coef <- cbind(Estimate = round(params,4), "Std.Error" = round(ssd,4),
                "z value" = round(params/ssd,4),
                "Pr(>|z|)" = round(2 * pnorm(-abs(params/ssd)),4))
  dimnames(Vcov) <- list(names(params), names(params))
  
  ## Information Criteria (AIC, BIC values from Davidson (2015))
  LL <- whittle$w
  AIC <- -LL - n.params
  BIC <- -LL - (n.params* log(n))/2
  
  IC <- data.frame(rbind("Log-Likelihood" = -LL, 
                         "AIC" = AIC,
                         "BIC" = BIC))
  colnames(IC)[1] <- "IC"
  
  ## output results
# print(list(coef, IC))
  
  ## save results as list 
  whit.result <- list(n = n, p=p, q=q,
                      coefficients = coef, vcov = Vcov,
                      IC = IC, In = In, sif = whittle$sif)
  
}

############# Examples ############
# using arfimaSim to generate series
# 
# Note: limits on MLE are between -0.99 and 0.99 while 
# limit of FML is (-0.5, 0.5). First estimate can be 
# considered exploratory, then difference if necessary. 
#
# df <- arfimaSim(200, d=.4)
# whittleFML(df)
# same as whittleFML(df, inits = list(d=0))
# 
# ardf <- arfimaSim(200, d=.75, p = .2)
# whittleFML(ardf, p=1, inits = list(d=0, AR=0))
#
# madf <- arfimaSim(200, d=.3, q = .4)
# whittleFML(madf, q=1, inits = list(d=0, MA=0))

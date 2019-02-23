# ------------------------------------------ #
# Frequency Domain Whittle Maximum Likelihood 
# Estimator of Fractional Integration
# 
# Main Function: whittleFML
# ARFIMA estimation requires: specDens; whitFunc
#
# Autocovariances and residuals require:
# tacvfARFIMA, tacvfARMA, tacvfFD, symtacvf, 
# mix, DLResiduals

# ------------------------------------------ #
## Autocovariances: calculate theoretical auto-
## covariance functions of FI and ARMA series 
## and mix. 
## 
## Borrowed (and adapted) from Veenstra and McLeod: 
## Veenstra, J. and McLeod, A. I. (Working Paper):
## The arfima R package: Exact Methods for Hyperbolic 
## Decay Time Series
# ------------------------------------------ # 
tacvfARFIMA <-function(phi, theta, d, maxlag) {
  x <- tacvfFD(d = d, maxlag = maxlag)
  y  <-  tacvfARMA(phi, theta, maxlag)
  z <- mix(x, y)
  return(z[1:(maxlag+1)])
}

## autocovariance of ARMA part
tacvfARMA <- function(phi, theta, maxlag) {
  p <- length(phi)
  q <- length(theta)
  maxlagp1 <- maxlag + 1
  
  if(max(p, q) == 0) {
    return(c(numeric(maxlag)))
  }
  r <- max(p, q) + 1
  b <- numeric(r)
  C <- numeric(q + 1)
  
  C[1] <- 1
  theta2 <- c(-1, theta)
  phi2 <- numeric(3 * r)
  phi2[r] <- -1
  if(p > 0) {
    phi2[r + 1:p] <- phi
  }
  if(q > 0) {
    for(k in 1:q) {
      C[k + 1] <- - theta[k]
      if(p > 0) {
        for(i in 1:min(p, k)) {
          C[k + 1] <- C[k + 1] + phi[i] * C[k + 1 - i]
        }
      }
    }
  }
  
  for(k in 0:q) {
    for(i in k:q) {
      b[k + 1] <- b[k + 1] - theta2[i + 1] * C[i - k + 1]
    }
  }
  
  if(p == 0) {
    g   <-   c(b, numeric(maxlagp1))[1:maxlagp1]
    return(g)
  }
  else if(p > 0) {
    a <- matrix(numeric(r^2), ncol = r)
    for(i in 1:r) {
      for(j in 1:r) {
        if(j == 1) {
          a[i, j] <- phi2[r + i - 1]
        }
        else if(j != 1) {
          a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 2]
        }
      }
    }
    
    g   <-   solve(a,  - b)
    
    if(length(g) <= maxlag) {
      g <- c(g, numeric(maxlagp1 - r)) 
      
      for(i in (r + 1):maxlagp1) {
        g[i] <- phi %*% g[i - 1:p]
      }
      return(g[1:maxlagp1])
    }
    else if(length(g) >= maxlagp1) {
      return(g[1:maxlagp1])
    }
  }
}

## autocovariance of ARFIMA (0,d,0)
tacvfFD <- function(d, maxlag) {
  x <- numeric(maxlag + 1)
  x[1] <- gamma(1 - 2 * d)/gamma(1 - d)^2
  for(i in 1:maxlag) 
    x[i + 1] <- ((i - 1 + d)/(i - d)) * x[i]
  x
}

## convolution
symtacvf <- function(x) {
  c(rev(x[-1])[-1], x)
}
mix <- function(x, y) {
  n <- 2*length(x)-2
  rev(Re(fft(fft(symtacvf(x)) * fft(symtacvf(y)), inverse = TRUE)/n)[(n/2 - 1):(n - 1)])
}

## calculate residuals
DLResiduals <- function(r,z) {
  n <- length(z)
  error <- numeric(n)
  sigmasq <- numeric(n)
  error[1] <- z[1]
  sigmasq[1] <- r[1]
  phi <- r[2]/r[1]
  error[2] <- z[2] - phi * z[1]
  sigmasqkm1 <- r[1] * (1 - phi^2)
  sigmasq[2] <- sigmasqkm1
  for(k in 2:(n - 1)) {
    phikk <- (r[k + 1] - phi %*% rev(r[2:k]))/sigmasqkm1
    sigmasqk <- sigmasqkm1 * (1 - phikk^2)
    phinew <- phi - phikk * rev(phi)
    phi <- c(phinew, phikk)
    sigmasqkm1 <- sigmasqk
    error[k + 1] <- z[k + 1] - crossprod(phi, rev(z[1:k]))
    sigmasq[k + 1] <- sigmasqk
  }
  res<-error
  res 
}


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
                       inits = list(d=0, AR = numeric(0), MA = numeric(0))) {
  
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
  params <- unlist(inits) # c(d, ar[1:p], ma[1:q]) # where p=0, q=0 is possible
  d <- params[1]
  n.params <- length(params)
  
  ## generate periodogram using fast fourier transform, drop 1st
  In <- (Mod(fft(x))^2/(2*pi*n)) [2:((n+1) %/% 2)]
  
  ## define function to be minimized by ML estimator 
  FMLwrap <- function(params) {
    r <- whitFunc(params, n = n, In = In, pq = pq, give.w.only = TRUE)
  }
  
  ## MLE method depends on number of parameters to estimate
  if (n.params == 1) { # one dimensional - ARFIMA (0,d,0) 
    result <- optim(par = params, fn = FMLwrap, lower = -.99, upper = .99,
                    method = "Brent", hessian = TRUE)
    params <- c(d = result$par)
  } 
  else { # ARFIMA (p,d,q)
    result <- optim(par = params, fn = FMLwrap, lower = -.99, upper = .99, 
                    method = "L-BFGS-B", hessian = TRUE)
    if (is.null(names(result$par) == "MA")) {
      params <- result$par
    } else {
      params <- result$par
      params[which(names(params) == "MA")] <- params[which(names(params) == "MA")]*-1
    }
    
  }
  
  ## Estimated Results
  whittle <- whitFunc(params, n, In, pq)
  
  if  (pq[1] == 0 && pq[2] == 0) {
    phi <- 0     
    theta <- 0
  }
  else if (pq[1] == 0 && pq[2] != 0) { 
    theta = params[2:(sum(pq)+1)]
    phi = 0
  }
  else if (pq[1] != 0 && pq[2] != 0) {
    phi = params[2:(2+pq[1]-1)]
    theta = params[(2+pq[1]):(sum(pq)+1)]
  }
  else if (pq[1] != 0 && pq[2] == 0) {
    phi = params[2:2+pq[1]-1]
    theta = 0;       
  }
  d <- params[1]
  
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
  
  ## Theoretical autocovariances and residuals
  rr <- tacvfARFIMA(phi = phi, theta = theta, d = d, maxlag = n-1)
  res <- DLResiduals(rr, x)
  
  ## output results
  #print(list(coef, IC))
  
  ## save results as list 
  whit.result <- list(n = n, p=p, q=q,
                      coefficients = coef, vcov = Vcov,
                      IC = IC, In = In, sif = whittle$sif,
                      rr = rr, residuals = res) 
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

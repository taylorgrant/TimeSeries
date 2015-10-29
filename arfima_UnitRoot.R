# ------------------------------------- # 
# Monte Carlo simulation of Dickey Fuller 
# and KPSS tests. Returns mean rate 
# of rejecting the null hypothesis. 
# DF null: unit root
# KPSS null: stationarity
# NOTE: 1000 MC is slow b/c of for loops
# ------------------------------------- #

source("arfimaSim.R")
require(tseries)

ns <- seq(100,500,100);  # number of observations
d1s = seq(.05, 1, .05)  # d -- after fractional differencing, series is at this level of d
rep = 1000 
nd <- length(d1s)

# initialize empty matrix
# dftest - null of unit root; kpss - null of stationarity
dftest <- replicate(5, diag(0,nd,rep), simplify=F)
kptest = replicate(5, diag(0,nd,rep), simplify=F)

ptm <- proc.time()

for (j in 1:length(ns)) {
  n <- ns[j]

  for (i in 1:nd) {
    d1 <- d1s[i]

      for (m in 1:rep) {
        
        x1 <- arfimaSim(n, d1)
        adf <- suppressWarnings(adf.test(x1, alternative = c("stationary"), 
                                          k = trunc((length(x1)-1)^(1/3)))) 
        
        # k = lags of ADF: 100=4;200=5;300=6;400=7;500=7
        dftest[[j]][i,m] <- ifelse(adf[4] < .05, 1, 0)
        rownames(dftest[[j]]) <- d1s
        x <- sapply(dftest, rowMeans) 
        
        # KPSS lag parameter: 100=2; 200=3; 300=3; 400=4; 500=5
        kpp <- suppressWarnings(kpss.test(x1))
        kptest[[j]][i,m] <- ifelse(kpp[3][1] < .05, 1, 0) # why is kpps named?
        rownames(kptest[[j]]) <- d1s
        y <- sapply(kptest, rowMeans)
      }
    }
}
proc.time() - ptm

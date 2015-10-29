MC_whittleFML <- function(ts, d, n, AR=NULL, MA=NULL) {
  # --------------------------------------- #
  ## Dependencies: MC_whittleFML.R; arfimaSim.R
  ##
  ## Inputs: 
  ## ts: length(s) of time series
  ## d: order(s) of fractional integration, after fractional
  ##    differencing, series at order of I(d) [-0.499, 0.499]
  ## n: number of simulations
  ## AR,MA: additional noise parameter [-.99, .99]
  ##        AR or MA, not both at same time
  ##        one AR/MA parameter value at a time
  # --------------------------------------- #  
  set.seed(1234)
  ts <- ts 
  d1s = d  
  rep = n 
  nd <- length(d1s)
  nts <- length(ts)

  # initialize empty matrix
  mat1 <- replicate(4, diag(0,rep, nd), simplify=FALSE)
  # mat1: each list is n; each row of list is d
  
    for (j in 1:length(ts)) {
      n <- ts[j]
    
      for (i in 1:nd) {
        d1 <- d1s[i]
      
        for (m in 1:rep) {
          if(!is.null(AR)) {
            x1 <- arfimaSim(n, d=d1, p=AR)
            dw <- whittleFML(x1, p=length(AR), inits=list(d0=0, AR=rep(0,times =length(AR))))
            mat1[[j]][m,i] <- dw[[4]][1] # place d estimates in matrix (dhat)
            colnames(mat1[[j]]) <- paste(d1s, AR, sep="-")
          } 
          else if(!is.null(MA)) {
            x1 <- arfimaSim(n, d=d1, q=MA)
            dw <- whittleFML(x1, q=length(MA), inits = list(d0=0, MA = rep(0, times=length(MA))))
            mat1[[j]][m,i] <- dw[[4]][1] # place d estimates in matrix (dhat)
            colnames(mat1[[j]]) <- d1s
          } 
          else {
          x1 <- arfimaSim(n, d1) # simulate ARFIMA series
          dw <- whittleFML(x1)   # estimate ARFIMA parameters
          mat1[[j]][m,i] <- dw[[4]][1] # place d estimates in matrix (dhat)
          colnames(mat1[[j]]) <- d1s
          }
        }
      }
    }
      
    # --------------------------------#
    # Results
    mat2 <- do.call(cbind.data.frame, mat1)
    mu_dhat <- apply(mat2, 2, mean)
    sd <- apply(mat2, 2, sd)
    d1 <- rep((d1s), times=nd)
    n <- rep((ts), each=nts)
    bias <- mean - d1
    mse <- bias^2 + sd^2
    rmse <- sqrt(mse)
    var <- sd^2
    AR <- ifelse(!is.null(AR), AR, NA)
    MA <- ifelse(!is.null(MA), MA, NA)
    
    full <- cbind(n, d1, mu_dhat, bias, sd, mse, rmse, var, AR,MA)
    full <- list(suppressWarnings(data.frame(round(full,3))))
}

# ------------------------------------------- #
## Examples
# 
## ARFIMA (0,d,0)
## z0_d_z0 <- MC_whittleFML(ts = seq(60,220,40), n=1000, 
##            d = c(-.25,0,.25,.4))
#
## ARFIMA (1,d,0)
## AR_d_z0 <- MC_whittleFML(ts = seq(60,220,40), n=1000, 
##            d = c(-.25,0,.25,.4), AR = .4)
#
## ARFIMA (0,d,1)
## z0_d_MA <- MC_whittleFML(ts = seq(60,220,40), n=1000, 
##            d = c(-.25,0,.25,.4), MA = .4)

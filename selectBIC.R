selectBIC <- function(x) {
## BIC selection for the whittleFML estimator
## note selectBIC differences the variable
## by default. whittleFML does not...
  n_params <- expand.grid(p=0:1, q=0:1)
  sel_bic <- rep(0, nrow(n_params))
  name <- rep(0, nrow(n_params))
  
  for (i in seq_along(sel_bic)) {
    sel_bic[i] <- whittleFML(diff(x), p=n_params[i,1], q=n_params[i,2], 
                             inits = list(d=0, AR=c(rep(0,n_params[i,1])),
                                          MA = c(rep(0,n_params[i,2]))))$IC[3,]
    
    name[i] <- paste("(", n_params[i,1],",", "d",",", n_params[i,2], ")", sep="")
    }
  
  ans <- which.max(sel_bic)
  out <- data.frame(BIC = round(sel_bic, 4))
  row.names(out) <- name
  print(out)
  print(paste("BIC optimized by Model", ans, name[ans], sep= " " ))

}

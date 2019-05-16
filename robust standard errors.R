#Function cse is useful to apply robust standard errors in case of heteroscedasticity
#Applied with OLS
cse = function(reg) {
  rob = sqrt(diag(vcovHC(reg, type = "HC1")))
  return(rob)
}


#clse function for correct standard errors
#Applied in case of panel data (with FE method)
#clustered SEs, clustered on "group"... could also cluster on "time" 
clse = function(reg) { 
  # index(reg, "id") returns the id or entity variable vector 
  G = length(unique(index(reg,"id")))
  N = length(index(reg,"id"))
  dfa = (G/(G - 1))   # note Bluhm multiplies this by finite-sample df adjustment
  rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                             cluster = "group")))
  return(rob)
}

#ivse function is applied for robust standard errors with TSLS method
ivse = function(reg) {
  rob = robust.se(reg)[,2]
  return(rob)
}
#density function for lifetime
f_density = function(u, mu){
  
  size = (omega - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  return( dbinom(x, size, mu) )
  
}

#link function
mu = function(eta){
  
  return( exp( eta ) / ( 1 + exp(eta) ) )
  
}

#covariates
eta = function(xi, beta){
  
  return( t(xi) %*% beta )
  
}

g_density = function(v, gv){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v -Delta
    return( gv[vec_idx] )
  }
  
}

X_support = c( (Delta + 1) : (omega) )
Y_support = c( (Delta + 1) : (Delta + m) )


alpha = function( g_parameters, xi, beta ){
  
  X_probs = sapply(X_support, f_density, mu = mu(eta(xi, beta)) )
  Y_probs = sapply(Y_support, g_density, gv = g_parameters )
  
  return(
    sum( (Y_probs %*% t(X_probs)) * outer(Y_support, X_support, "<=") )
  )
  
}

#inputs


#needs density functions from likelihood-testing.R
sim_sample.seed = function(n, omega, Delta, m, beta, gv, x.cov, seed.start = 1){
  
  #global input to function
  X_support = c( (Delta + 1) : (omega) )
  Y_support = c( (Delta + 1) : (Delta + m) )
  
  Xi = c()
  Yi = c()
  #begin function
  for(i in c(1:n)){
    
    #build current version of h-star
    xi = as.numeric(x.cov[i,])
    
    xpmf <- sapply(X_support, f_density, mu = mu(eta(xi, beta)) )
    ypmf <- sapply(Y_support, g_density, gv = gv)
    
    marginal <- list(cumsum(xpmf)[-length(X_support)],
                     cumsum(ypmf)[-length(Y_support)])
    support <- list(unique(X_support), unique(Y_support))
    
    ind <- outer(Y_support, X_support, "<=")
    joint <- outer(ypmf, xpmf)
    
    hInv <- data.frame(x = X_support[col(ind)[ind]],
                       y = Y_support[row(ind)[ind]],
                       u = cumsum(joint[ind] / alpha(gv, xi, beta)))
    
    set.seed(i + n *(seed.start - 1))
    smp <- findInterval(runif(1), hInv$u) + 1
    Xi <- append(Xi, hInv$x[smp])
    Yi <- append(Yi, hInv$y[smp])
    
  }
  
  xy_dat = data.frame("Xi" = Xi, "Yi" = Yi)
  
  return(cbind(xy_dat, x.cov))
  
}
#inputs

#support of Y
#number of betas
#sample data



#helper functions

################################################################################
#second derivatives
# B_function.star = function(v, v.star, theta, yi, xi){
#   
#   g.parameters = theta[1:(m - 1)]
#   beta.parameters = theta[m:length(theta)]
#   
#   if(v == v.star){
#     idx = v - (Delta + 1) + 1
#     a1 = -(yi == v) * (1 / g.parameters[idx]^2)
#   }
#   else{
#     a1 = 0
#   }
#   
#   a2 = -(yi == (Delta + m)) / ((1 - sum(g.parameters))^2 )
#   
#   a3.1 = surv.v(v, xi, beta.parameters) - surv.v(Delta + m, xi, beta.parameters)
#   a3.2 = surv.v(v.star, xi, beta.parameters) - surv.v(Delta + m, xi, beta.parameters)
#   a3 = a3.1 * a3.2
#   
#   a4 = (1 / (alpha.star(g.parameters, xi, beta.parameters)))^2
#   
#   return(a1 + a2 + a3 * a4)
#   
# }
B_function.star = function(theta, yi, zi){
  
  g.param = theta[1:(m - 1)]
  beta = theta[m:length(theta)]
  
  S1 = diag( as.numeric(
    sapply(c((Delta + 1):(Delta + m - 1)), S_function, mu = mu(eta(zi, beta))) -
      S_function(Delta + m, mu = mu(eta(zi, beta)))
  ))
  
  S2 = matrix(rep( sapply(c((Delta + 1):(Delta + m - 1)), S_function, mu = mu(eta(zi, beta))),
                   length(Y_support) - 1),
              nrow = length(Y_support) - 1,
              ncol = length(Y_support) - 1,
              byrow = TRUE)
  
  S3 = S_function( Delta + m, mu = mu(eta(zi, beta)) ) * matrix(1,
                                                                nrow = length(Y_support) - 1,
                                                                ncol = length(Y_support) - 1)
  
  return(
    diag( as.numeric(
      -(1 * (c((Delta + 1):(Delta + m - 1)) == yi)) / g.param^2
    )) -
      (yi == (Delta + m)) / ((1 - sum(g.param))^2 ) * matrix(1,
                                                             nrow = length(Y_support) - 1,
                                                             ncol = length(Y_support) - 1) +
      (1 / alpha.star(g.param, zi, beta))^2 * (S1 %*% (S2 - S3))
  )
  
}

r_function.star = function(zi, theta){
  
  g.parameters = theta[1:(m - 1)]
  beta.parameters = theta[m:length(theta)]
  
  a2i = c()
  for(v in c( (Delta + 1):(Delta + m - 1)) ){
    
    a2.1 = sum(sapply(c(v:omega), q_function, zi = zi, beta = beta.parameters))
    a2.2 = theta[(v - Delta)]
    a2i = append(a2i, a2.2 * a2.1)
    
  }
  
  a2 = sum(a2i)
  
  a3.1 = sum(sapply(c((m + Delta):omega), q_function, zi = zi, beta = beta.parameters))
  a3 = (1 - sum(theta[1:(m - 1)])) * a3.1
  
  return(a2 + a3)
  
}


A_function.star = function(u, zi, theta){
  
  g.parameters = theta[1:(m - 1)]
  beta.parameters = theta[m:length(theta)]
  
  a1 = s_function(u, zi, beta.parameters) / f_density(u, mu(eta(zi, beta.parameters)) )
  a2 = (q_function(u, zi, beta.parameters) / f_density(u, mu(eta(zi, beta.parameters)) ))^2
  
  a3 = (r_function.star(zi, theta) / alpha.star(g.parameters, zi, beta.parameters) )^2
  
  a4i = c()
  for(v in c( (Delta + 1):(Delta + m - 1)) ){
    
    a4.1 = sum(sapply(c(v:omega), s_function, zi = zi, beta = beta.parameters))
    a4.2 = theta[(v - Delta)]
    a4i = append(a4i, a4.2 * a4.1)
    
  }
  
  a4 = sum(a4i)
  
  a5.1 = sum(sapply(c((m + Delta):omega), s_function, zi = zi, beta = beta.parameters))
  a5 = (1 - sum(theta[1:(m - 1)])) * a5.1
  
  a6 = (1 / alpha.star(g.parameters, zi, beta.parameters))
  
  return( a1 - a2 + a3 - a6 * (a4 + a5) )
  
}

#off-diagonals
# A_function.star = function(v, theta, xi){
#   
#   g.parameters = theta[1:(m - 1)]
#   beta.parameters = theta[m:length(theta)]
#   
#   a1 = r_function.star(xi, theta)
#   a2 = surv.v(v, xi, beta.parameters) - surv.v(Delta + m, xi, beta.parameters)
#   a3 = alpha.star(g.parameters, xi, beta.parameters)
#   
#   a4.1 = sum(sapply(c((m + Delta):omega), q_function, xi = xi, beta = beta.parameters))
#   a4.2 = sum(sapply(c(v:omega), q_function, xi = xi, beta = beta.parameters))
#   a4 = a4.2 - a4.1
#   
#   return( (a1 * a2) / (a3^2) - a4 / a3 )
#   
# }
D_function.star = function(theta, zi){
  
  g.param = theta[1:(m - 1)]
  beta = theta[m:length(theta)]
  
  S =
    matrix( data = rep( sapply( c((Delta + 1):(omega)),
                                q_function,
                                zi,
                                beta),
                        length(Y_support) - 1),
            ncol = length(Y_support) - 1,
            nrow = length(X_support)
    )
  S = colSums(S * outer(X_support, Y_support[-length(Y_support)], ">="))
  
  return(
    (r_function.star(zi, theta) / alpha.star(g.param, zi, beta)^2) *
      (
        sapply( c((Delta + 1):(Delta + m - 1)),
                S_function,
                mu = mu(eta(zi, beta)) ) -
          S_function( (m + Delta), mu(eta(zi, beta)) )
      ) -
      (1 / alpha.star(g.param, zi, beta)) *
      S +
      (1 / alpha.star(g.param, zi, beta)) *
      sum(sapply(c((m + Delta):omega), q_function, zi, beta))
  )
  
}

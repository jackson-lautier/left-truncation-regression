#inputs

#support of Y
#number of betas
#sample data

#omega, Delta, m, and data: Xi, Yi, x.cov
#tolerance?
#NR iterations limit?
#step 4 iterations limit?

# X.df = cbind("x0" = rep(1,n), sim.data[,c("x1", "x2")])
# X = as.matrix(X.df)
# J = rep(1, n)

tol = 1/10e7 #indv NR methods
global.tol = 1/10e6 #for iteration step 4


################################################################################
#helper functions
################################################################################


#Econometrics and Statistics: initial estimation of g, beta parameters
g_est = function(g){
  
  idx = g - Delta
  
  if( g == (Delta + m)){
    return( rev_haz_est[idx] )
  }
  if( ((Delta + 1) <= g) & (g < Delta + m ) ){
    return(rev_haz_est[idx] * prod(1 - rev_haz_est[(idx + 1):(m)]))
  }
  
}

Cnx <- function(x) {
  ind_x <- ifelse( (sim.data$Yi <= x) & (x <= sim.data$Xi),1,0 )
  return( (1/nrow(sim.data)) * (sum(ind_x)) )
}

bnx <- function(y) {
  return( (1/nrow(sim.data)) * (1/Cnx(y)) * (sum(sim.data$Yi == y)) )
}

lnx <- function(x) {
  return( (1/nrow(sim.data)) * (1/Cnx(x)) * (sum(sim.data$Xi == x)) )
}

fnx <- function(x) {
  return(sum(sim.data$Xi == x) / nrow(sim.data))
}

f_est = function(u){
  
  idx = u - Delta
  
  if( u == (Delta + 1)){
    return( haz_est[idx] )
  }
  if( ((Delta + 1) < u) & (u <= omega ) ){
    return(haz_est[idx] * prod(1 - haz_est[1:(idx - 1)]))
  }
  
}

#NR algorithm for the betas
#derivatives
df_dmu = function(u, mu){
  
  size = (omega - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  a = choose(size, x)
  
  num = ((1 - mu)^(size - x)) * (mu^(x-1)) * (size * mu - x)
  den = mu - 1
  
  return( a * (num / den) )
  
}

dmu_deta = function(eta){
  
  num = exp( eta )
  den = ( exp(eta) + 1 )^2
  
  return( num / den)
  
}

q_function = function(u, zi, beta){
  
  a = df_dmu(u, mu(eta(zi, beta)) )
  b = dmu_deta( eta(zi, beta) )
  
  return(a * b)
  
}

W_function = function(u, g_parameters, zi, beta){
  
  A = q_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) )
  
  Y_probs = sapply(Y_support, g_density, gv = g_parameters )
  q_results = sapply(X_support, q_function, zi = zi, beta = beta)
  
  c1 = sum( (Y_probs %*% t(q_results)) * outer(Y_support, X_support, "<=") )
  
  B = c1 / alpha(g_parameters, zi, beta)
  
  return(A - B)
  
}

lB_function = function(beta){
  
  w_data = c()
  for(i in c(1:n)){
    w_data = append(w_data,
                    W_function(u = sim.data$Xi[i],
                               g_parameters = gv,
                               zi = as.numeric(Z.df[i,]),
                               beta = beta))
  }
  
  W = diag(w_data)
  
  return( t(Z) %*% W %*% J )
  
  
}

d2mu_deta2 = function(eta){
  
  num = (exp( eta ) - 1) * exp( eta )
  den = ( exp(eta) + 1 )^3
  
  return( num / den)
  
}

d2f_dmu2 = function(u, mu){
  
  size = (omega - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  a = choose(size, x)
  
  v1 = a * ((1 - mu)^(size - x)) * (mu^(x - 2))
  v2 = (size^2 - size) * mu^2 + (2 - 2*size) * x * mu + x^2 - x
  num = v1 * v2
  den = (mu - 1)^2
  
  return(num / den)
  
}

r_function = function( g_parameters, zi, beta){
  
  Y_probs = sapply(Y_support, g_density, gv = g_parameters )
  q_results = sapply(X_support, q_function, zi = zi, beta = beta)
  
  c1 = sum( (Y_probs %*% t(q_results)) * outer(Y_support, X_support, "<=") )
  
  return(c1)
  
}

s_function = function(u, zi, beta){
  
  v1 = d2f_dmu2(u, mu(eta(zi, beta)) )
  v2 = (dmu_deta( eta(zi, beta) ))^2
  
  v3 = df_dmu(u, mu(eta(zi, beta)) )
  v4 = d2mu_deta2( eta(zi, beta) )
  
  return( v1 * v2 - v3 * v4 )
  
}


A_function = function(u, g_parameters, zi, beta){
  
  A = s_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) )
  B = (q_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) ))^2
  
  Y_probs = sapply(Y_support, g_density, gv = g_parameters )
  s_results = sapply(X_support, s_function, zi = zi, beta = beta)
  c1 = sum( (Y_probs %*% t(s_results)) * outer(Y_support, X_support, "<=") )
  
  C = c1 / alpha(g_parameters, zi, beta)
  D = (r_function(g_parameters, zi, beta) / alpha(g_parameters, zi, beta))^2
  
  return(A - B - C + D)
  
}

l2B_function = function(beta){
  
  a_data = c()
  for(i in c(1:n)){
    a_data = append(a_data,
                    A_function(u = sim.data$Xi[i],
                               g_parameters = gv,
                               zi = as.numeric(Z.df[i,]),
                               beta = beta))
  }
  
  A = diag(a_data)
  
  return(t(Z) %*% A %*% Z)
  
  
}

#NR method for g parameters
alpha.star = function( g_parameters, zi, beta ){
  
  g.parameters.star = c(g_parameters, 1-sum(g_parameters))
  
  X_probs = sapply(X_support, f_density, mu = mu(eta(zi, beta)) )
  Y_probs = sapply(Y_support, g_density, gv = g.parameters.star )
  
  return(
    sum( (Y_probs %*% t(X_probs)) * outer(Y_support, X_support, "<=") )
  )
  
}

#g-parameters derivatives
# m_function = function(g_parameters, yi, xi, beta){
#   
#   M = matrix(NA, nrow = 1, ncol = length(Y_support) - 1)
#   
#   for(v in (Y_support[-length(Y_support)])){
#     
#     idx = v - (Delta + 1) + 1
#     a1 = (yi == v)
#     a2 = g_parameters[idx] / alpha.star(g_parameters, xi, beta)
#     a3 = sum(sapply(c(v:omega), f_density, mu = mu(eta(xi, beta))))
#     
#     M[1, idx] = (a1 - a2 * a3)
#     
#   }
#   
#   return(M)
#   
# }
# 
# B_function = function(v, v.star, g_parameters, yi, xi, beta){
#   
#   if(v == v.star){
#     idx = v - (Delta + 1) + 1
#     a1 = sum(sapply(c(v:omega), f_density, mu = mu(eta(xi, beta))))
#     a2 = g_parameters[idx]
#     a3 = alpha.star(g_parameters, xi, beta)
#     a4 = sum(sapply(c((m+Delta):omega), f_density, mu = mu(eta(xi, beta))))
#     
#     return( -(a1 * a3 - a2 * a1 * (a1 - a4)) / (a3^2) )
#   }
#   else{
#     idx = v - (Delta + 1) + 1
#     a1 = sum(sapply(c(v:omega), f_density, mu = mu(eta(xi, beta))))
#     a2 = g_parameters[idx]
#     a3 = alpha.star(g_parameters, xi, beta)
#     a4 = sum(sapply(c(v.star:omega), f_density, mu = mu(eta(xi, beta))))
#     a5 = sum(sapply(c((m+Delta):omega), f_density, mu = mu(eta(xi, beta))))
#     
#     return( (a2 * a1 * (a4 - a5)) / (a3^2) )
#   }
#   
# }
# 
# F_function = function( g.param ){
#   
#   M = matrix(NA, nrow = nrow(sim.data), ncol = length(Y_support) - 1)
#   
#   #g_parameters = c(g.param, 1 - sum(g.param))
#   
#   for(i in c(1:nrow(sim.data))){
#     
#     M[i,] = m_function(g_parameters = g.param,
#                        yi = sim.data$Yi[i],
#                        xi = as.numeric(X.df[i,]),
#                        beta = beta.est)
#     
#   }
#   
#   return( (t(M) %*% J) )
#   
# }
# 
# F2g_function = function( g.param ){
#   
#   B = matrix(NA, nrow = length(Y_support) - 1, ncol = length(Y_support) - 1)
#   
#   #g_parameters = c(g.param, 1 - sum(g.param))
#   
#   for(v in (Y_support[-length(Y_support)])){
#     for(v.star in (Y_support[-length(Y_support)])){
#       b = c()
#       for(i in c(1:nrow(sim.data))){
#         
#         b = append(b,
#                    B_function(v, v.star,
#                               g_parameters = g.param,
#                               yi = sim.data$Yi[i],
#                               xi = as.numeric(X.df[i,]),
#                               beta = beta.est))
#         
#       }
#       B[(v - Delta),(v.star - Delta)] = sum(b)
#     }
#   }
#   return(B)
#   
# }

B_function = function(g.param, yi, zi, beta){
  
  S1 = diag( as.numeric(
    (g.param *
       sapply(c((Delta + 1):(Delta + m - 1)), S_function, mu = mu(eta(zi, beta))))
  ))
  
  S2 = matrix(rep( sapply(c((Delta + 1):(Delta + m - 1)), S_function, mu = mu(eta(zi, beta))),
                   length(Y_support) - 1),
              nrow = length(Y_support) - 1,
              ncol = length(Y_support) - 1,
              byrow = TRUE)
  
  S3 = S_function( Delta + m, mu = mu(eta(zi, beta)) ) * matrix(1,
                                                                nrow = length(Y_support) - 1,
                                                                ncol = length(Y_support) - 1)
  
  S4 = diag (as.numeric(
    ( alpha.star(g.param, zi, beta) *
        sapply(c((Delta + 1):(Delta + m - 1)), S_function, mu = mu(eta(zi, beta))))
  ))
  
  return(
    (1 / alpha.star(g.param, zi, beta)^2) *
      (S1 %*% (S2 - S3) - S4)
  )
  
  
}

F_function = function( g.param ){
  
  S.mat = matrix(NA, nrow = nrow(sim.data),
                 ncol = length(Y_support) - 1)
  
  for(i in c(1:nrow(sim.data))){
    
    S.mat[i, ] = sapply(c((Delta + 1):(Delta + m - 1)),
                        S_function, mu = mu(eta(zi = as.numeric(Z.df[i,]),
                                                beta = beta.est))) /
      alpha.star(g.param,
                 zi = as.numeric(Z.df[i,]),
                 beta = beta.est)
  }
  
  Ind.mat = matrix(0, nrow = length(sim.data$Yi), ncol = max(sim.data$Yi) - Delta)
  Ind.mat[cbind(1:length(sim.data$Yi), (sim.data$Yi - Delta))] <- 1
  Ind.mat = Ind.mat[,-(max(sim.data$Yi) - Delta)]
  
  G.mat = diag(as.numeric(g.param))
  
  return( t(Ind.mat - S.mat %*% G.mat) %*% J )
  
}

F2g_function = function( g.parameters ){
  
  B = matrix(0, nrow = length(Y_support) - 1, ncol = length(Y_support) - 1)
  
  for(i in c(1:nrow(sim.data))){
    
    B = B + B_function( g.param = g.parameters,
                        yi = sim.data$Yi[i],
                        zi = as.numeric(Z.df[i,]),
                        beta = beta.est)
    
  }
  
  return(B)
  
}

#for step 4 iteration

#h-star function
h_star = function(u , v, zi, gv, beta){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + m)) &
      (v <= u)){
    
    return(
      (f_density(u, mu = mu(eta(zi, beta))) * g_density(v, gv)) / alpha(gv, zi, beta)
    )
  }
  else{
    return(0)
  }
  
}

#needs to be generalized if sim study works
log.like.function = function(g1, g2, beta0, beta1, beta2){
  
  THETA = c(g1, g2, 1 - g1 - g2, beta0, beta1, beta2)
  data = sim.data
  
  gv = THETA[1:3]
  beta = THETA[4:6]
  
  L = c()
  
  for(i in c(1:n)){
    
    cXi = data[i, "Xi"]
    cYi = data[i, "Yi"]
    cxi = c(1, data[i, "x1"], data[i, "x2"])
    
    L = append(L, log( h_star(cXi, cYi, cxi, gv, beta) ) )
    
  }
  
  return(sum(L))
  
}

#derivative functions
surv.v = function(v, zi, beta){
  
  sum(sapply(c(v:omega), f_density, mu = mu(eta(zi, beta))))
  
}

dl.dgv = function(v, theta.est){
  
  gv = theta.est[1:(m - 1)]
  beta = theta.est[m:length(theta.est)]
  
  idx = v - (Delta + 1) + 1
  dl.cur = c()
  for(i in c(1:nrow(sim.data))){
    
    zi = as.numeric(Z.df[i,])
    
    a1 = (sim.data$Yi[i] == v) / (theta.est[idx])
    a2 = (sim.data$Yi[i] == (Delta + m)) / (1 - sum(theta.est[1:(m - 1)]))
    a3 = surv.v(v, zi, beta) - surv.v(Delta + m, zi, beta)
    a4 = alpha.star(gv, zi, beta)
    
    dl.cur = append( dl.cur, a1 - a2 - a3 / a4 )
    
  }
  
  return(sum(dl.cur))
  
}

W_function.star = function(u, zi, theta){
  
  beta = theta[m:length(theta)] 
  
  a1 = q_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) )
  
  a2i = c()
  for(v in c( (Delta + 1):(Delta + m - 1)) ){
    
    a2.1 = sum(sapply(c(v:omega), q_function, zi = zi, beta = beta))
    a2.2 = theta[(v - Delta)]
    a2i = append(a2i, a2.2 * a2.1)
    
  }
  
  a2 = sum(a2i)
  
  a3.1 = sum(sapply(c((m + Delta):omega), q_function, zi = zi, beta = beta))
  a3 = (1 - sum(theta[1:(m - 1)])) * a3.1
  
  a4 = alpha.star(theta[1:(m - 1)], zi, beta)
  
  return( a1 - (1/a4) * (a2 + a3) )
  
}
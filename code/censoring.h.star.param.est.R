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

#aoas paper
gnv.aoas = function(v){
  return((1/nrow(sim.data)) * sum(sim.data$Yi == v))
}

f_X.aoas = function(u, THETA){
  
  p = THETA[1]
  
  size = (xi - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  return( dbinom(x, size, p) )
  
}

df_dp.aoas = function(u, THETA){
  
  p = THETA[1]
  
  A = f_X.aoas(u, p)
  B = (u - (Delta + 1)) / p
  C = (xi - u) / (1 - p)
  
  return(A * (B - C))
  
}

dS_dp.aoas = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return( sum(sapply(c(u:xi), df_dp.aoas, THETA)) )
  }
  else{
    return(0)
  }
  
}

P_constraint = function(p_input){
  
  n = nrow(sim.data)
  
  LHS = c()
  for(k in c((Delta + 1):(Delta + m))){
    A = gnv.aoas(k)
    B = sum(mapply(f_X.aoas, c(k:xi), p_input))
    C = sum(mapply(df_dp.aoas, c(k:xi), p_input))
    LHS = append(LHS, (A / B) * C )
  }
  
  RHS = c()
  for(i in c(1:n)){
    if(sim.data$Di[i] == 1){
      A = df_dp.aoas(sim.data$Ti[i], p_input)
      B = f_X.aoas(sim.data$Ti[i], p_input)
    }
    if(sim.data$Di[i] == 0){
      A = dS_dp.aoas(sim.data$Ti[i] + 1, p_input)
      B = S_X.aoas(sim.data$Ti[i] + 1, p_input)
    }
    
    RHS = append(RHS, A/B)
    
  }
  
  return( ((sum(RHS) / n) - sum(LHS))^2 )
  
}

S_X.aoas = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return( sum(sapply(c(u:xi), f_X.aoas, THETA)) )
  }
  else{
    return(0)
  }
  
}

g_tau_hat = function(v, p_input){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = gnv.aoas(v) / S_X.aoas(v, p_input)
  B = mapply(gnv.aoas, c(v_min:v_max))
  C = mapply(S_X.aoas, c(v_min:v_max), p_input)
  return( A * (sum( B/C ))^(-1) )
  
}

#derivatives
df_dmu = function(u, mu){
  
  size = (xi - (Delta + 1) + 1) - 1
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


W_tau_function = function(u, d, g_parameters, zi, beta){
  
  A =
    ifelse(d == 1,
           q_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) ),
           sum(sapply(c((u+1):xi), q_function, zi = zi, beta = beta)) /
             S_function(u + 1, mu(eta(zi, beta))))
  
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
                    W_tau_function(u = sim.data$Ti[i],
                                   d = sim.data$Di[i],
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
  
  size = (xi - (Delta + 1) + 1) - 1
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


A_tau_function = function(u, d, g_parameters, zi, beta){
  
  if(d == 1){
    A = s_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) )
    B = (q_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) ))^2
  }
  
  if(d == 0){
    A = sum(sapply(c((u+1):xi), s_function, zi = zi, beta = beta)) /
      S_function(u + 1, mu(eta(zi, beta)))
    B = (sum(sapply(c((u+1):xi), q_function, zi = zi, beta = beta)) /
           S_function(u + 1, mu(eta(zi, beta))))^2
  }
  
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
                    A_tau_function(u = sim.data$Ti[i],
                                   d = sim.data$Di[i],
                                   g_parameters = gv,
                                   zi = as.numeric(Z.df[i,]),
                                   beta = beta))
  }
  
  A = diag(a_data)
  
  return(t(Z) %*% A %*% Z)
  
  
}

alpha.star = function( g_parameters, zi, beta ){
  
  g.parameters.star = c(g_parameters, 1-sum(g_parameters))
  
  X_probs = sapply(X_support, f_density, mu = mu(eta(zi, beta)) )
  Y_probs = sapply(Y_support, g_density, gv = g.parameters.star )
  
  return(
    sum( (Y_probs %*% t(X_probs)) * outer(Y_support, X_support, "<=") )
  )
  
}

#g-parameters derivatives
# m_function = function(g_parameters, yi, zi, beta){
#   
#   M = matrix(NA, nrow = 1, ncol = length(Y_support) - 1)
#   
#   for(v in (Y_support[-length(Y_support)])){
#     
#     idx = v - (Delta + 1) + 1
#     a1 = (yi == v)
#     a2 = g_parameters[idx] / alpha.star(g_parameters, zi, beta)
#     a3 = sum(sapply(c(v:xi), f_density, mu = mu(eta(zi, beta))))
#     
#     M[1, idx] = (a1 - a2 * a3)
#     
#   }
#   
#   return(M)
#   
# }

dF_dgv = function(v, g_parameters){
  
  idx = v - (Delta + 1) + 1
  
  ret = c()
  for(i in c(1:n)){
    
    a1 = sim.data$Yi[i] == v
    a2 = (g_parameters[idx] / alpha.star(g_parameters, zi = as.numeric(Z.df[i,]), beta.est))
    a3 = sum(sapply(c(v:xi), f_density, mu = mu(eta(zi = as.numeric(Z.df[i,]), beta.est))))
    
    ret = append(ret, a1 - a2 * a3)  
  }
  
  return(sum(ret))
  
}

# B_function = function(v, v.star, g_parameters, yi, zi, beta){
#   
#   if(v == v.star){
#     idx = v - (Delta + 1) + 1
#     a1 = sum(sapply(c(v:xi), f_density, mu = mu(eta(zi, beta))))
#     a2 = g_parameters[idx]
#     a3 = alpha.star(g_parameters, zi, beta)
#     a4 = sum(sapply(c((m+Delta):xi), f_density, mu = mu(eta(zi, beta))))
#     
#     return( -(a1 * a3 - a2 * a1 * (a1 - a4)) / (a3^2) )
#   }
#   else{
#     idx = v - (Delta + 1) + 1
#     a1 = sum(sapply(c(v:xi), f_density, mu = mu(eta(zi, beta))))
#     a2 = g_parameters[idx]
#     a3 = alpha.star(g_parameters, zi, beta)
#     a4 = sum(sapply(c(v.star:xi), f_density, mu = mu(eta(zi, beta))))
#     a5 = sum(sapply(c((m+Delta):xi), f_density, mu = mu(eta(zi, beta))))
#     
#     return( (a2 * a1 * (a4 - a5)) / (a3^2) )
#   }
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

#numeric derivative checks
dF_dg1 = function(g1, g2){
  
  g_parameters = c(g1, g2)
  v = 1
  
  idx = v - (Delta + 1) + 1
  
  ret = c()
  for(i in c(1:n)){
    
    a1 = sim.data$Yi[i] == v
    a2 = (g_parameters[idx] / alpha.star(g_parameters, zi = as.numeric(Z.df[i,]), beta.est))
    a3 = sum(sapply(c(v:xi), f_density, mu = mu(eta(zi = as.numeric(Z.df[i,]), beta.est))))
    
    ret = append(ret, a1 - a2 * a3)  
  }
  
  return(sum(ret))
  
}

dF_dg2 = function(g1, g2){
  
  g_parameters = c(g1, g2)
  v = 2
  
  idx = v - (Delta + 1) + 1
  
  ret = c()
  for(i in c(1:n)){
    
    a1 = sim.data$Yi[i] == v
    a2 = (g_parameters[idx] / alpha.star(g_parameters, zi = as.numeric(Z.df[i,]), beta.est))
    a3 = sum(sapply(c(v:xi), f_density, mu = mu(eta(zi = as.numeric(Z.df[i,]), beta.est))))
    
    ret = append(ret, a1 - a2 * a3)  
  }
  
  return(sum(ret))
  
}

#functions for NR method
# F_function = function( g.param ){
# 
#   M = matrix(NA, nrow = nrow(sim.data), ncol = length(Y_support) - 1)
# 
#   for(i in c(1:nrow(sim.data))){
# 
#     M[i,] = m_function(g_parameters = g.param,
#                        yi = sim.data$Yi[i],
#                        zi = as.numeric(Z.df[i,]),
#                        beta = beta.est)
# 
#   }
# 
#   return( (t(M) %*% J) )
# 
# }

#new
# F_function = function( g.param ){
#   
#   process_row.f <- function(row) {
#     
#     m_function(g_parameters = g.param,
#                yi = row["Yi"],
#                zi = c(row["z0"],
#                       row["z1"],
#                       row["z2"],
#                       row["z3"],
#                       row["z4"]),
#                beta = beta.est)
#     
#   }
# 
#   M = apply(sim.data, 1, process_row.f)
# 
#   return( M %*% J )
# 
# }

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


# F2g_function = function( g.param ){
# 
#   B = matrix(NA, nrow = length(Y_support) - 1, ncol = length(Y_support) - 1)
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
#                               zi = as.numeric(Z.df[i,]),
#                               beta = beta.est))
# 
#       }
#       B[(v - Delta),(v.star - Delta)] = sum(b)
#     }
#   }
#   return(B)
# 
# }
# F2g_function = function( g.param ){
# 
#   B = matrix(NA, nrow = length(Y_support) - 1, ncol = length(Y_support) - 1)
# 
#   for(v in (Y_support[-length(Y_support)])){
#     for(v.star in (Y_support[-length(Y_support)])){
# 
#       process_row <- function(row) {
# 
#         B_function(v = v, v.star = v.star,
#                    g_parameters = g.param,
#                    yi = row["Yi"],
#                    zi = c(row["z0"],
#                           row["z1"],
#                           row["z2"],
#                           row["z3"],
#                           row["z4"]),
#                    beta = beta.est)
#       }
# 
#       B[(v - Delta),(v.star - Delta)] = sum(apply(sim.data, 1, process_row))
#     }
#   }
#   return(B)
# 
# }
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
      (u <= xi) &
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

#h-bar function
h_bar = function(u, v, zi, gv, beta){
  
  a = S_function(u + 1, mu = mu(eta(zi, beta)))
  b = g_density(v, gv)
  return( (a*b) / alpha(gv, zi, beta))
  
}


log.like.function = function(g1, g2, beta0, beta1, beta2){
  
  THETA = c(g1, g2, 1 - g1 - g2, beta0, beta1, beta2)
  data = sim.data
  
  gv = THETA[1:3]
  beta = THETA[4:6]
  
  L = c()
  
  D.0 = data[data$Di == 0,]
  D.1 = data[data$Di == 1,]
  
  for(i in c(1:(nrow(D.0)))){
    
    cTi = D.0[i, "Ti"]
    cYi = D.0[i, "Yi"]
    czi = c(1, D.0[i, "z1"], D.0[i, "z2"])
    
    L = append(L, log( h_bar(cTi, cYi, czi, gv, beta) ) )
    
  }
  
  for(i in c(1:(nrow(D.1)))){
    
    cTi = D.1[i, "Ti"]
    cYi = D.1[i, "Yi"]
    czi = c(1, D.1[i, "z1"], D.1[i, "z2"])
    
    L = append(L, log( h_star(cTi, cYi, czi, gv, beta) ) )
    
  }
  
  return(sum(L))
  
}

#derivative functions
surv.v = function(v, zi, beta){
  
  sum(sapply(c(v:xi), f_density, mu = mu(eta(zi, beta))))
  
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

W_tau_function.star = function(u, d, zi, theta){
  
  beta = theta[m:length(theta)]
  
  a1 = ifelse(d == 1,
              q_function(u, zi, beta) / f_density(u, mu(eta(zi, beta)) ),
              sum(sapply(c((u+1):xi), q_function, zi = zi, beta = beta)) /
                S_function(u + 1, mu(eta(zi, beta))))
  
  #a1 = q_function(u, xi, beta) / f_density(u, mu(eta(xi, beta)) )
  
  a2i = c()
  for(v in c( (Delta + 1):(Delta + m - 1)) ){
    
    a2.1 = sum(sapply(c(v:xi), q_function, zi = zi, beta = beta))
    a2.2 = theta[(v - Delta)]
    a2i = append(a2i, a2.2 * a2.1)
    
  }
  
  a2 = sum(a2i)
  
  a3.1 = sum(sapply(c((m + Delta):xi), q_function, zi = zi, beta = beta))
  a3 = (1 - sum(theta[1:(m - 1)])) * a3.1
  
  a4 = alpha.star(theta[1:(m - 1)], zi, beta)
  
  return( a1 - (1/a4) * (a2 + a3) )
  
}

g_Y.aoas = function(v, THETA){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v + 1 - Delta
    return( THETA[vec_idx] )
  }
  
}

alpha.aoas = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(xi))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y.aoas, THETA = THETA)
    res = append(res,
                 f_X.aoas(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}
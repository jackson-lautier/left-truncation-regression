h_u_v = function(u,v){
  
  n = nrow(obs_data)
  val = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
  return(sum(val) / n)
  
}

h_dot_v = function(v){
  
  n = nrow(obs_data)
  val = c()
  for(u in c(v:omega)){
    cur = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
    val = append(val, cur/n)
  }
  
  return(sum(val))
}

f_X.aoas = function(u, THETA){
  
  p = THETA[1]
  
  size = (omega - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  return( dbinom(x, size, p) )
  
}

#note THETA can be one dimensional 
df_dp = function(u, THETA){
  
  p = THETA[1]
  
  A = f_X.aoas(u, p)
  B = (u - (Delta + 1)) / p
  C = (omega - u) / (1 - p)
  
  return(A * (B - C))
  
}

P_constraint = function(p_input){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  LHS = c()
  for(k in c((Delta + 1):(Delta + m))){
    A = h_dot_v(k)
    B = sum(mapply(f_X.aoas, c(k:omega), p_input))
    C = sum(mapply(df_dp, c(k:omega), p_input))
    LHS = append(LHS, (A / B) * C )
  }
  
  RHS = c()
  for(k in c((Delta + 1):(Delta + m))){
    for(j in c(k:omega)){
      A = h_u_v(j,k)
      B = f_X.aoas(j, p_input)
      C = df_dp(j, p_input)
      RHS = append(RHS, (A / B) * C)
    }
  }
  
  return( (sum(RHS) - sum(LHS))^2 )
  
}

S_X.aoas = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( sum(sapply(c(u:omega), f_X.aoas, THETA)) )
  }
  else{
    return(0)
  }
  
}

g_tau_hat = function(v, p_input){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = h_dot_v(v) / S_X.aoas(v, p_input)
  B = mapply(h_dot_v, c(v_min:v_max))
  C = mapply(S_X.aoas, c(v_min:v_max), p_input)
  return( A * (sum( B/C ))^(-1) )
  
}

g_Y.aoas = function(v, THETA){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v + 1 - Delta
    return( THETA[vec_idx] )
  }
  
}

alpha.aoas = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y.aoas, THETA = THETA)
    res = append(res,
                 f_X.aoas(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}
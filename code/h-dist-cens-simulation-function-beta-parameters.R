#density function for lifetime
f_density = function(u, mu){
  
  size = (omega - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  return( dbinom(x, size, mu) )
  
}

S_function = function(u, mu){
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return(sum(sapply(c(u:xi),
                      f_density,
                      mu)))
  } else {
    return( 0 )
  }
  
}

#link function
mu = function(eta){
  
  return( exp( eta ) / ( 1 + exp(eta) ) )
  
}

#covariates
eta = function(zi, beta){
  
  return( t(zi) %*% beta )
  
}

g_density = function(v, gv){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v - Delta
    return( gv[vec_idx] )
  }
  
}

tau = e - (m + Delta + 1)
xi = min(omega, e - 1)
X_support = c( (Delta + 1) : (xi) )
Y_support = c( (Delta + 1) : (Delta + m) )


alpha = function( g_parameters, zi, beta ){
  
  X_probs = sapply(X_support, f_density, mu = mu(eta(zi, beta)) )
  Y_probs = sapply(Y_support, g_density, gv = g_parameters )
  
  return(
    sum( (Y_probs %*% t(X_probs)) * outer(Y_support, X_support, "<=") )
  )
  
}

#inputs
# n = 1000 #sample size
# 
# #trapezoid support
# omega = 4
# m = 3
# Delta = 0
# e = 5

#parameters
# gv = c(0.3, 0.5, 0.2) #left-trun RV
# beta = c(0.5, 2, 1)
# 
# #covariates
# set.seed(1000)
# z.cov = data.frame(z0 = rep(1,n),
#                    z1 = rnorm(n, 0, 1),
#                    z2 = rnorm(n, 0, 1))



#generate support
# tau = e - (m + Delta + 1)
# xi = min(omega, e - 1)
# X_support = c( (Delta + 1) : (xi) )
# Y_support = c( (Delta + 1) : (Delta + m) )

#create h-bar

#all possible outcomes
# ind <- outer(Y_support, X_support, "<=")
# 
# h.bar = data.frame(Y.i = Y_support[row(ind)[ind]],
#                    T.i = X_support[col(ind)[ind]],
#                    D.i = 0,
#                    R.i = NA,
#                    g.Yi = sapply(Y.i, g_density, gv = gv),
#                    S.Ti1 = sapply((T.i + 1),
#                                   S_function,
#                                   mu = mu(eta(zi, beta)),
#                                   pi = pi,
#                                   f.02 = f.02),
#                    valid.obs = 1 * (T.i == Y.i + tau))
# h.bar$prob = h.bar$valid.obs *
#   ((h.bar$g.Yi * h.bar$S.Ti1) / alpha(gv, zi, beta, f.02, pi))
# 
# h.star.1 = data.frame(Y.i = Y_support[row(ind)[ind]],
#                       T.i = X_support[col(ind)[ind]],
#                       D.i = 1,
#                       R.i = 1,
#                       g.Yi = sapply(Y.i, g_density, gv = gv),
#                       f1.Ti = sapply(T.i,
#                                      f.01_density,
#                                      mu = mu(eta(zi, beta))),
#                       valid.obs = 1 * (T.i <= Y.i + tau))
# h.star.1$prob = pi * h.star.1$valid.obs *
#   ((h.star.1$g.Yi * h.star.1$f1.Ti) / alpha(gv, zi, beta, f.02, pi))
# 
# h.star.2 = data.frame(Y.i = Y_support[row(ind)[ind]],
#                       T.i = X_support[col(ind)[ind]],
#                       D.i = 1,
#                       R.i = 2,
#                       g.Yi = sapply(Y.i, g_density, gv = gv),
#                       f2.Ti = sapply(T.i,
#                                      f.02_density,
#                                      f.02 = f.02),
#                       valid.obs = 1 * (T.i <= Y.i + tau))
# h.star.2$prob = (1 - pi) * h.star.2$valid.obs *
#   ((h.star.2$g.Yi * h.star.2$f2.Ti) / alpha(gv, zi, beta, f.02, pi))
# 
# 
# h.dens =
#   rbind(h.bar[,c("Y.i", "T.i", "D.i", "R.i", "prob")],
#         h.star.1[,c("Y.i", "T.i", "D.i", "R.i", "prob")],
#         h.star.2[,c("Y.i", "T.i", "D.i", "R.i", "prob")])
# 
# h.dens$U = cumsum(h.dens$prob)
# h.dens = h.dens[h.dens$prob > 0,]
# 
# sim.U = runif(1)
# smp <- findInterval(sim.U, h.dens$U) + 1
# 
# Ti <- h.dens$T.i[smp]
# Yi <- h.dens$Y.i[smp]
# Di <- h.dens$D.i[smp]
# Ri <- h.dens$R.i[smp]
# c(Yi, Ti, Di, Ri)
# 
# Ti = c()
# Yi = c()
# Di = c()
# Ri = c()
# 
# Ti <- append(Ti, h.dens$T.i[smp])
# Yi <- append(Yi, h.dens$Y.i[smp])
# Di <- append(Di, h.dens$D.i[smp])
# Ri <- append(Ri, h.dens$R.i[smp])

sim_sample.seed = function(n, omega, Delta, m, e,
                           beta, gv, z.cov, seed.start = 1){
  
  #global input to function
  tau = e - (m + Delta + 1)
  xi = min(omega, e - 1)
  X_support = c( (Delta + 1) : (xi) )
  Y_support = c( (Delta + 1) : (Delta + m) )
  
  Ti = c()
  Yi = c()
  Di = c()
  
  #begin function
  for(i in c(1:n)){
    
    zi = as.numeric(z.cov[i,])
    
    #build current version of h-dist w cens
    ind <- outer(Y_support, X_support, "<=")
    
    Y.i = Y_support[row(ind)[ind]]
    T.i = X_support[col(ind)[ind]]
    
    h.bar = data.frame(Y.i = Y.i,
                       T.i = T.i,
                       D.i = 0,
                       g.Yi = sapply(Y.i, g_density, gv = gv),
                       S.Ti1 = sapply((T.i + 1),
                                      S_function,
                                      mu = mu(eta(zi, beta))),
                       valid.obs = 1 * (T.i == Y.i + tau))
    h.bar$prob = h.bar$valid.obs *
      ((h.bar$g.Yi * h.bar$S.Ti1) / alpha(gv, zi, beta))
    
    h.star = data.frame(Y.i = Y.i,
                        T.i = T.i,
                        D.i = 1,
                        g.Yi = sapply(Y.i, g_density, gv = gv),
                        f.Ti = sapply(T.i,
                                       f_density,
                                       mu = mu(eta(zi, beta))),
                          valid.obs = 1 * (T.i <= Y.i + tau))
    h.star$prob = h.star$valid.obs *
      ((h.star$g.Yi * h.star$f.Ti) / alpha(gv, zi, beta))
    
    h.dens =
      rbind(h.bar[,c("Y.i", "T.i", "D.i", "prob")],
            h.star[,c("Y.i", "T.i", "D.i", "prob")])
    
    h.dens$U = cumsum(h.dens$prob)
    h.dens = h.dens[h.dens$prob > 0,]
    
    set.seed(i + n *(seed.start - 1))
    smp <- findInterval(runif(1), h.dens$U) + 1
    Ti <- append(Ti, h.dens$T.i[smp])
    Yi <- append(Yi, h.dens$Y.i[smp])
    Di <- append(Di, h.dens$D.i[smp])
    
  }
  
  xtd_dat = data.frame("Yi" = Yi,
                       "Ti" = Ti,
                       "Di" = Di)
  
  return(cbind(xtd_dat, z.cov))
  
}


# df = sim_sample.seed(n = 1000,
#                 omega,
#                 Delta,
#                 m,
#                 e,
#                 beta,
#                 gv,
#                 z.cov,
#                 seed.start = 1)


######################################################################################
######################################################################################
######################################################################################
# Data analysis scripts to reproduce application results in the manuscript:
# 
# "Discrete Time-to-Event Regression Analysis Under Left-Truncation
# with Applications to Consumer Finance"
#
# LAUTIER, POZDNYAKOV, YAN
# 2026
#
# Computer and R version details
# _                                
# platform       x86_64-w64-mingw32               
# arch           x86_64                           
# os             mingw32                          
# crt            ucrt                             
# system         x86_64, mingw32                  
# status                                          
# major          4                                
# minor          5.1                              
# year           2025                             
# month          06                               
# day            13                               
# svn rev        88306                            
# language       R                                
# version.string R version 4.5.1 (2025-06-13 ucrt)
# nickname       Great Square Root
#
######################################################################################
######################################################################################
#
# RUN-TIME WARNING: may exceed 10 hours for 1,000 reps depending on computing power
#
######################################################################################
######################################################################################
######################################################################################
######################################################################################
# INSTRUCTIONS
#
# supporting files:
# "./processed-data/aart-2017-37mo.csv
#
# './code/h-dist-cens-simulation-function-beta-parameters.R')
# './code/censoring.h.star.param.est.R')
# './code/censoring.hessian.std.error.R')
#
#The code must be run sequentially downwards.
#Results will appear in a newly created 'results' folder.
#
#
######################################################################################
######################################################################################
######################################################################################
######################################################################################

require('ggplot2')
require('extrafont')

######################################################################################
######################################################################################
######################################################################################
######################################################################################

#where results will be stored
dir.create('./results/')

################################################################################
rm(list=ls())

#set-up of problem

#inverse link
inv.link = function(p.est){
  
  return( -log( 1 / p.est - 1) )
}

#sample size
n = 1000

#trapezoid
omega = 8
Delta = 0
m = 5
e = 10

#parameters
beta.true = c(0.5, 0, 0, 0, 0)
gv.true = c(0.3, 0.25, 0.2, 0.15, 0.10)

#there is a seed in this sim data formula for reproducibility
source('./code/h-dist-cens-simulation-function-beta-parameters.R')
source('./code/censoring.h.star.param.est.R')
source('./code/censoring.hessian.std.error.R')


num.reps = 1000
rep.start = 1

Omega = rep(NA,1000)

start_time <- Sys.time()
for(r in c(rep.start:num.reps)){
  
  ##############################################################################
  # [1] generate random sample when H0 is true
  ##############################################################################
  
  #covariates
  set.seed(r)
  z.cov = data.frame(z0 = rep(1,n),
                     z1 = rnorm(n, 0, 0.1),
                     z2 = rnorm(n, 0, 0.1),
                     z3 = rnorm(n, 0, 0.1),
                     z4 = rnorm(n, 0, 0.1))
  
  sim.data = sim_sample.seed(n,
                             omega,
                             Delta,
                             m,
                             e,
                             beta.true,
                             gv.true,
                             z.cov,
                             seed.start = r)
  
  ##############################################################################
  # [2] estimate unrestricted param MLE
  ##############################################################################
  
  Z.df = cbind("z0" = rep(1,n), sim.data[,c("z1", "z2", "z3", "z4")])
  Z = as.matrix(Z.df)
  J = rep(1, n)
  
  #step 1: get initial values
  p.est = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
  g.0 = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p.est)
  beta.0 = inv.link(p.est)
  
  #step 2: use NR method to get beta estimates
  gv = g.0 
  B0 = c(beta.0, rep(0,length(beta.true)-1))
  
  B.hist = matrix(NA, nrow = length(B0), ncol = 200)
  B.hist[,1] = B0
  B.est = B0
  
  #Newton-Raphson method:
  for(i in 2:100){
    
    B.new = B.est - solve(l2B_function(B.est)) %*% lB_function(B.est)
    B.hist[,i] = B.new
    
    if(max(abs(B.hist[,i]-B.hist[,i-1])) < tol){B.est <- B.hist[,i-1]; break}
    else(B.est = B.new)
    
    #print(i)
    #print(B.hist[,1:i])
    
  }
  
  idx = min(which(is.na(B.hist[1,]))) - 1
  beta.est = B.hist[,idx]
  
  #step 3: update the gv estimates
  G0 = g.0[1:(length(g.0) - 1)]
  
  G.hist = matrix(NA, nrow = length(G0), ncol = 200)
  G.hist[,1] = G0
  G.est = G0
  
  #Newton-Raphson method:
  for(i in 2:100){
    
    G.new = G.est - solve(F2g_function(G.est)) %*% F_function(G.est)
    G.hist[,i] = G.new
    
    if(max(abs(G.hist[,i]-G.hist[,i-1])) < tol){G.est <- G.hist[,i-1]; break}
    else(G.est = G.new)
    
    #print(i)
    #print(G.hist[,1:i])
    
  }
  
  idx = min(which(is.na(G.hist[1,]))) - 1
  g.est = c(G.hist[,idx], 1 - sum(G.hist[,idx]))
  
  #step 4: iterate until convergence
  s4.g0 = g.est
  s4.beta.0 = beta.est
  
  cur.g.est = s4.g0
  cur.beta.est = s4.beta.0
  
  
  THETA.hist = matrix(NA,
                      nrow = length(c(cur.g.est, cur.beta.est)),
                      ncol = 200)
  THETA.hist[,1] = c(cur.g.est, cur.beta.est)
  
  step.g.est = cur.g.est
  step.B.est = cur.beta.est
  
  for(j in 2:100){
    
    #update parameters values from previous estimate
    gv = step.g.est
    B0 = step.B.est
    
    B.hist = matrix(NA, nrow = length(B0), ncol = 200)
    B.hist[,1] = B0
    B.est = B0
    
    #Newton-Raphson method:
    for(i in 2:100){
      
      B.new = B.est - solve(l2B_function(B.est)) %*% lB_function(B.est)
      B.hist[,i] = B.new
      
      if(max(abs(B.hist[,i]-B.hist[,i-1])) < tol){B.est <- B.hist[,i-1]; break}
      else(B.est = B.new)
      
      #print(i)
      #print(B.hist[,1:i])
      
    }
    
    idx = min(which(is.na(B.hist[1,]))) - 1
    new.beta.est = B.hist[,idx]
    
    #calc updated derivatives of log-like
    w_data = c()
    for(i in c(1:n)){
      w_data = append(w_data,
                      W_tau_function.star(u = sim.data$Ti[i],
                                          d = sim.data$Di[i],
                                          zi = as.numeric(Z.df[i,]),
                                          theta = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)))
    }
    
    W = diag(w_data)
    
    cur.deriv = c(sapply(c((Delta + 1):(Delta + m - 1)), dl.dgv,
                         theta.est = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)),
                  t(t(Z) %*% W %*% J))
    
    THETA.hist[,j] = c(step.g.est, new.beta.est)
    
    if(max(abs(cur.deriv)) < global.tol){
      step.beta.est <- THETA.hist[c((length(step.g.est)+1):(length(c(cur.g.est, cur.beta.est)))),j-1]; break}
    else(step.beta.est = new.beta.est)
    
    #update g parameters
    G0 = gv[1:(length(step.g.est)-1)]
    beta.est = step.beta.est
    
    G.hist = matrix(NA, nrow = length(G0), ncol = 200)
    G.hist[,1] = G0
    G.est = G0
    
    #Newton-Raphson method:
    for(i in 2:100){
      
      G.new = G.est - solve(F2g_function(G.est)) %*% F_function(G.est)
      G.hist[,i] = G.new
      
      if(max(abs(G.hist[,i]-G.hist[,i-1])) < tol){G.est <- G.hist[,i-1]; break}
      else(G.est = G.new)
      
      #print(i)
      #print(G.hist[,1:i])
      
    }
    
    idx = min(which(is.na(G.hist[1,]))) - 1
    new.g.est = c(G.hist[,idx], 1 - sum(G.hist[,idx]))
    
    w_data = c()
    for(i in c(1:n)){
      w_data = append(w_data,
                      W_tau_function.star(u = sim.data$Ti[i],
                                          d = sim.data$Di[i],
                                          zi = as.numeric(Z.df[i,]),
                                          theta = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)))
    }
    
    W = diag(w_data)
    
    cur.deriv = c(sapply(c((Delta + 1):(Delta + m - 1)), dl.dgv,
                         theta.est = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)),
                  t(t(Z) %*% W %*% J))
    
    THETA.hist[,j] = c(new.g.est, step.beta.est)
    
    if(max(abs(cur.deriv)) < global.tol){step.g.est <- THETA.hist[c(1:(length(step.g.est))),j-1]; break}
    else(step.g.est = new.g.est)
    
    if( sum(abs(THETA.hist[,j-1] - THETA.hist[j])) == 0 ){break}
    
    #print(j)
    #print( THETA.hist[,j] )
    #print( cur.deriv )
    
  }
  ################################################################STOPPING POINT
  idx = min(which(is.na(THETA.hist[1,]))) - 1
  
  theta.est = THETA.hist[,idx]
  
  ##############################################################################
  # [3] calculate unrestricted likelihood
  ##############################################################################
  
  L = c()
  data = sim.data
  
  D.0 = data[data$Di == 0,]
  D.1 = data[data$Di == 1,]
  
  for(i in c(1:(nrow(D.0)))){
    
    cTi = D.0[i, "Ti"]
    cYi = D.0[i, "Yi"]
    czi = c(1, 
            D.0[i, "z1"],
            D.0[i, "z2"],
            D.0[i, "z3"],
            D.0[i, "z4"])
    
    L = append(L, log( h_bar(cTi, cYi, czi,
                             gv = theta.est[1:5],
                             beta = theta.est[6:10]) ) )
    
  }
  
  for(i in c(1:(nrow(D.1)))){
    
    cTi = D.1[i, "Ti"]
    cYi = D.1[i, "Yi"]
    czi = c(1, 
            D.1[i, "z1"],
            D.1[i, "z2"],
            D.1[i, "z3"],
            D.1[i, "z4"])
    
    L = append(L, log( h_star(cTi, cYi, czi,
                              gv = theta.est[1:5],
                              beta = theta.est[6:10]) ) )
    
  }
  
  ell.1 = sum(L)
  
  ##############################################################################
  # [4] estimate restricted param MLE (AOAS)
  ##############################################################################
  
  p_hat = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
  G_hat = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)
  theta.est.0 = c(p_hat, G_hat)
  
  ##############################################################################
  # [5] calculate restricted likelihood
  ##############################################################################
  
  Li = c()
  for(i in c(1:n)){
    
    if(sim.data$Di[i] == 1){
      contrib = g_Y.aoas(sim.data$Yi[i], theta.est.0) * f_X.aoas(sim.data$Ti[i], theta.est.0)
      Li = append(Li, log(contrib))
    }
    
    if(sim.data$Di[i] == 0){
      contrib = g_Y.aoas(sim.data$Yi[i], theta.est.0) * S_X.aoas(sim.data$Ti[i] + 1, theta.est.0)
      Li = append(Li, log(contrib))
    }
  }
  
  n = nrow(sim.data)
  alp = alpha.aoas(theta.est.0)
  
  ell.0 = -n * log(alp) + sum(Li)
  
  ##############################################################################
  # [6] calculate chi-sq test statistic
  ##############################################################################
  
  Omega[r] = -2 * (ell.0 - ell.1)
  
  ##############################################################################
  # [7] save results + update user
  ##############################################################################
  
  write.csv(Omega, "./results/cens-sim-study-results-manu-chi-sq.csv")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  if( (r/100) %in% c(1:(num.reps/100)) ){
    print( paste("complete reps: ", r, " of ", num.reps,
                 " ### elapsed time: ", round(elapsed_time, 5), sep = "") )
  }
  
}

#require('ggplot2')
#require('extrafont')

ss.results = read.csv("./results/cens-sim-study-results-manu-chi-sq.csv")
ss.results = ss.results[,-1]

deg.free = 4 #per sim study set-up

#empirical type I error
sum(ss.results > qchisq(0.95, deg.free)) / length(ss.results)

#plot results
df = data.frame("sim_result" = ss.results)

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dchisq, args = list(df = deg.free)) +
  #xlab(TeX(("$\\sqrt{ n }(\\hat{p}_n - p_0)$"))) +
  xlab(
    expression( Omega [' '*tau*', '*italic(n)]
                ~ ' (dashed) versus ' ~
                  chi^2
                ~ ' distribution (solid)' )) +
  ylab("Density Height") +
  theme_bw() +
  theme(axis.title.x=element_text(size=10, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=10, family="Times New Roman"),
        axis.text.x=element_text(size=10, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman"),
        legend.text=element_text(size=10, family="Times New Roman"),
        legend.position = "bottom")

ggsave("./manuscript/chi-sq-sim-cens.pdf",height=4,width=6,device = cairo_pdf)
#set-up of problem

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

#inverse link
inv.link = function(p.est){
  
  return( -log( 1 / p.est - 1) )
}

#specific to set-up
best.p = function(p){
  
  
  v1 = sapply(X = X_support - 1,
              dbinom,
              size = (omega - (Delta + 1) + 1) - 1,
              prob = p)
  
  return( sum( (v1 - U_est)^2 ) )
  
}

#sample size
n = 1000

#trapezoid
omega = 12
Delta = 0
m = 8

#parameters
beta.true = c(0.5, 0.5, 1, -1.5, -0.5)
gv.true = c(0.3, 0.2, 0.13, 0.10, 0.09, 0.07, 0.06, 0.05)

#there is a seed in this sim data formula for reproducibility
source('./code/h-star-simulation-function-beta-parameters.R')
source('./code/h.star.param.est.R')
source('./code/hessian.std.error.R')

num.reps = 1000
res.mat = matrix(NA, nrow = num.reps, ncol = 2 * (length(c(beta.true, gv.true)) - 1))

rep.start = 1 #for breaking into parts r = 1355 stopped

start_time <- Sys.time()
for(r in c(rep.start:num.reps)){
  
  ##############################################################################
  #generate random sample
  ##############################################################################
  
  #covariates
  set.seed(r)
  x.cov = data.frame(x0 = rep(1,n),
                     x1 = rnorm(n, 0, 0.1),
                     x2 = rnorm(n, 0, 0.1),
                     x3 = rnorm(n, 0, 0.1),
                     x4 = rnorm(n, 0, 0.1))
  
  sim.data = sim_sample.seed(n, omega, Delta, m, beta.true, gv.true, x.cov, seed.start = r)
  
  ##############################################################################
  #parameter estimation
  ##############################################################################
  
  X.df = cbind("x0" = rep(1,n), sim.data[,c("x1", "x2", "x3", "x4")])
  X = as.matrix(X.df)
  J = rep(1, n)
  
  #step 1: get initial values
  rev_haz_est = sapply(c( (Delta + 1) : (Delta + m)), bnx)
  G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
  g.0 = G_est
  
  haz_est = sapply(c( (Delta + 1) : omega), lnx)
  U_est = sapply(c( (Delta + 1) : omega), f_est)
  
  p.est = optimize(best.p, c(0,1))$minimum
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
                      W_function.star(u = sim.data$Xi[i],
                                      xi = as.numeric(X.df[i,]),
                                      theta = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)))
    }
    
    W = diag(w_data)
    
    cur.deriv = c(sapply(c((Delta + 1):(Delta + m - 1)), dl.dgv,
                         theta.est = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)),
                  t(t(X) %*% W %*% J))
    
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
                      W_function.star(u = sim.data$Xi[i],
                                      xi = as.numeric(X.df[i,]),
                                      theta = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)))
    }
    
    W = diag(w_data)
    
    cur.deriv = c(sapply(c((Delta + 1):(Delta + m - 1)), dl.dgv,
                         theta.est = c(step.g.est[1:(length(step.g.est)-1)], new.beta.est)),
                  t(t(X) %*% W %*% J))
    
    THETA.hist[,j] = c(new.g.est, step.beta.est)
    
    if(max(abs(cur.deriv)) < global.tol){step.g.est <- THETA.hist[c(1:(length(step.g.est))),j-1]; break}
    else(step.g.est = new.g.est)
    
    if( sum(abs(THETA.hist[,j-1] - THETA.hist[j])) == 0 ){break}
    
    #print(j)
    #print( THETA.hist[,j] )
    #print( cur.deriv )
    
  }
  
  idx = min(which(is.na(THETA.hist[1,]))) - 1
  
  theta.est = THETA.hist[,idx]
  
  #end_time <- Sys.time()
  #elapsed_time <- end_time - start_time
  #print(elapsed_time)
  
  ##############################################################################
  #standard error estimation
  ##############################################################################
  
  B = matrix(NA, nrow = (length(Y_support) - 1), ncol = (length(Y_support) - 1))
  for(v in Y_support[1:(m-1)]){
    for(v.star in Y_support[1:(m-1)]){
      b = c()
      for(i in c(1:nrow(sim.data))){
        
        b = append(b,
                   B_function.star(v, v.star,
                                   theta = theta.est[-length(gv.true)],
                                   yi = sim.data$Yi[i],
                                   xi = as.numeric(X.df[i,])))
        
      }
      B[v,v.star] = sum(b)
    }
  }
  
  #check
  #rownames(B) = c("g1", "g2")
  #colnames(B) = c("g1", "g2")
  #B
  
  z_data = c()
  for(i in c(1:n)){
    z_data = append(z_data,
                    Z_function.star(u = sim.data$Xi[i],
                                    xi = as.numeric(X.df[i,]),
                                    theta = theta.est[-length(gv.true)]))
  }
  
  Z = diag(z_data)
  
  #check
  #t(X) %*% Z %*% X
  
  dl.dgv.db = matrix(NA, nrow = m - 1, ncol = length(beta.true))
  #colnames(dl.dgv.db) = c("b0", "b1", "b2")
  #rownames(dl.dgv.db) = c("g1", "g2")
  
  for(v in c( (Delta + 1):(Delta + m - 1)) ){
    
    a_data = c()
    for(i in c(1:n)){
      a_data = append(a_data,
                      A_function.star(v = v,
                                      xi = as.numeric(X.df[i,]),
                                      theta = theta.est[-length(gv.true)]))
    }
    
    A = diag(a_data)
    
    dl.dgv.db[(v - Delta), ] = t(J) %*% A %*% X
    
  }
  
  #formulaic hessian
  H = cbind(rbind(B, t(dl.dgv.db)),
            rbind(dl.dgv.db, t(X) %*% Z %*% X))
  
  colnames(H) <- NULL
  rownames(H) <- NULL
  
  res.mat[r, ] = c(theta.est[-length(gv.true)], diag(solve(-H)))
  
  write.csv(res.mat, "./results/sim-study-results-manu.csv")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  if( (r/100) %in% c(1:(num.reps/10)) ){
    print( paste("complete reps: ", r, " of ", num.reps,
                 " ### elapsed time: ", round(elapsed_time, 5), sep = "") )
  }
  
}

#summarize results

ss.results = read.csv("./results/sim-study-results-manu.csv")

ss.results = ss.results[,-1]

# g5 = c()
# for(i in c(1:1000)){
#   
#   g5 = append(g5, 1 - sum( ss.results[i,1:7] ) )
#   
# }
# g6 = ss.results$V5
# g7 = ss.results$V6
# 
# ss.results$V5 = g5
# ss.results$V6 = g6
# ss.results$V7 = g7

p.lab = c(paste("g",c(1:7),sep=""),paste("b",c(0:4),sep=""))
se.lab = paste("se.", p.lab, sep = "")

colnames(ss.results) = c(p.lab, se.lab)

beta.true = c(0.5, 0.5, 1, -1.5, -0.5)
gv.true = c(0.3, 0.2, 0.13, 0.10, 0.09, 0.07, 0.06, 0.05)
gv.true = gv.true[-length(gv.true)]

manu.table = data.frame("true" = c(gv.true, beta.true))

manu.table$eBias = abs( colMeans(ss.results)[1:length(p.lab)] - c(gv.true, beta.true) )

# var(ss.results$g1); mean(ss.results$se.g1)
# var(ss.results$g2); mean(ss.results$se.g2)
# var(ss.results$g3); mean(ss.results$se.g3)
# var(ss.results$g4); mean(ss.results$se.g4)
# var(ss.results$g5); mean(ss.results$se.g5)
# var(ss.results$g6); mean(ss.results$se.g6)
# var(ss.results$g7); mean(ss.results$se.g7)
# var(ss.results$b0); mean(ss.results$se.b0)
# var(ss.results$b1); mean(ss.results$se.b1)
# var(ss.results$b2); mean(ss.results$se.b2)
# var(ss.results$b3); mean(ss.results$se.b3)
# var(ss.results$b4); mean(ss.results$se.b4)

manu.table$eSE = c(sd(ss.results$g1),
                   sd(ss.results$g2),
                   sd(ss.results$g3),
                   sd(ss.results$g4),
                   sd(ss.results$g5),
                   sd(ss.results$g6),
                   sd(ss.results$g7),
                   sd(ss.results$b0),
                   sd(ss.results$b1),
                   sd(ss.results$b2),
                   sd(ss.results$b3),
                   sd(ss.results$b4))

manu.table$Thm21 = c(mean( sqrt(ss.results$se.g1)),
                     mean( sqrt(ss.results$se.g2)),
                     mean( sqrt(ss.results$se.g3)),
                     mean( sqrt(ss.results$se.g4)),
                     mean( sqrt(ss.results$se.g5)),
                     mean( sqrt(ss.results$se.g6)),
                     mean( sqrt(ss.results$se.g7)),
                     mean( sqrt(ss.results$se.b0)),
                     mean( sqrt(ss.results$se.b1)),
                     mean( sqrt(ss.results$se.b2)),
                     mean( sqrt(ss.results$se.b3)),
                     mean( sqrt(ss.results$se.b4)))

boxplot(ss.results[,1:length(p.lab)])

#coverage probabilities
cp = c()

ci.low = ss.results$g1 - qnorm(0.975) * sqrt(ss.results$se.g1)
ci.high = ss.results$g1 + qnorm(0.975) * sqrt(ss.results$se.g1)
cp = append(cp, sum( (gv.true[1] > ci.low) & (gv.true[1] < ci.high) ) / 1000)

ci.low = ss.results$g2 - qnorm(0.975) * sqrt(ss.results$se.g2)
ci.high = ss.results$g2 + qnorm(0.975) * sqrt(ss.results$se.g2)
cp = append(cp, sum( (gv.true[2] > ci.low) & (gv.true[2] < ci.high) ) / 1000)

ci.low = ss.results$g3 - qnorm(0.975) * sqrt(ss.results$se.g3)
ci.high = ss.results$g3 + qnorm(0.975) * sqrt(ss.results$se.g3)
cp = append(cp, sum( (gv.true[3] > ci.low) & (gv.true[3] < ci.high) ) / 1000)

ci.low = ss.results$g4 - qnorm(0.975) * sqrt(ss.results$se.g4)
ci.high = ss.results$g4 + qnorm(0.975) * sqrt(ss.results$se.g4)
cp = append(cp, sum( (gv.true[4] > ci.low) & (gv.true[4] < ci.high) ) / 1000)

ci.low = ss.results$g5 - qnorm(0.975) * sqrt(ss.results$se.g5)
ci.high = ss.results$g5 + qnorm(0.975) * sqrt(ss.results$se.g5)
cp = append(cp, sum( (gv.true[5] > ci.low) & (gv.true[5] < ci.high) ) / 1000)

ci.low = ss.results$g6 - qnorm(0.975) * sqrt(ss.results$se.g6)
ci.high = ss.results$g6 + qnorm(0.975) * sqrt(ss.results$se.g6)
cp = append(cp, sum( (gv.true[6] > ci.low) & (gv.true[6] < ci.high) ) / 1000)

ci.low = ss.results$g7 - qnorm(0.975) * sqrt(ss.results$se.g7)
ci.high = ss.results$g7 + qnorm(0.975) * sqrt(ss.results$se.g7)
cp = append(cp, sum( (gv.true[7] > ci.low) & (gv.true[7] < ci.high) ) / 1000)

ci.low = ss.results$b0 - qnorm(0.975) * sqrt(ss.results$se.b0)
ci.high = ss.results$b0 + qnorm(0.975) * sqrt(ss.results$se.b0)
cp = append(cp, sum( (beta.true[1] > ci.low) & (beta.true[1] < ci.high) ) / 1000)

ci.low = ss.results$b1 - qnorm(0.975) * sqrt(ss.results$se.b1)
ci.high = ss.results$b1 + qnorm(0.975) * sqrt(ss.results$se.b1)
cp = append(cp, sum( (beta.true[2] > ci.low) & (beta.true[2] < ci.high) ) / 1000)

ci.low = ss.results$b2 - qnorm(0.975) * sqrt(ss.results$se.b2)
ci.high = ss.results$b2 + qnorm(0.975) * sqrt(ss.results$se.b2)
cp = append(cp, sum( (beta.true[3] > ci.low) & (beta.true[3] < ci.high) ) / 1000)

ci.low = ss.results$b3 - qnorm(0.975) * sqrt(ss.results$se.b3)
ci.high = ss.results$b3 + qnorm(0.975) * sqrt(ss.results$se.b3)
cp = append(cp, sum( (beta.true[4] > ci.low) & (beta.true[4] < ci.high) ) / 1000)

ci.low = ss.results$b4 - qnorm(0.975) * sqrt(ss.results$se.b4)
ci.high = ss.results$b4 + qnorm(0.975) * sqrt(ss.results$se.b4)
cp = append(cp, sum( (beta.true[5] > ci.low) & (beta.true[5] < ci.high) ) / 1000)

manu.table$cp = cp





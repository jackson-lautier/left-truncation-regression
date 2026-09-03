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
# RUN-TIME WARNING: may exceed 2 hours depending on computing power
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

######################################################################################
# 60 MONTHS
######################################################################################

rm(list=ls())

start_time <- Sys.time()

reg.data = read.csv("./processed-data/aart-2017-60mo-robust.csv")
reg.data = reg.data[,-1]

#trapezoid
omega = max(reg.data$Ti)
Delta = min( min(reg.data$Yi), min(reg.data$Ti) ) - 1
m = max(reg.data$Yi) - Delta
tau = unique( (reg.data$Ti - reg.data$Yi)[reg.data$Di == 0] )
e = tau + (m + Delta + 1)

#inverse link
inv.link = function(p.est){
  
  return( -log( 1 / p.est - 1) )
}

#sample size
n = nrow(reg.data)

#helper functions
source('./code/h-dist-cens-simulation-function-beta-parameters.R')
source('./code/censoring.h.star.param.est.R')
source('./code/censoring.hessian.std.error.R')

################################################################################
# begin regression analysis
################################################################################

sim.data = reg.data


##############################################################################
#parameter estimation
##############################################################################

Z.df = cbind("z0" = rep(1,n), sim.data[,c(4:ncol(sim.data))])
Z = as.matrix(Z.df)
J = rep(1, n)

#step 1: get initial values
p.est = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
g.0 = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p.est)
beta.0 = inv.link(p.est)

#step 2: use NR method to get beta estimates
gv = g.0 
B0 = c(beta.0, rep(0, ncol(Z) - 1))

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

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

##############################################################################
#standard error estimation
##############################################################################

# B = matrix(NA, nrow = (length(Y_support) - 1), ncol = (length(Y_support) - 1))
# for(v in Y_support[1:(m-1)]){
#   for(v.star in Y_support[1:(m-1)]){
#     b = c()
#     for(i in c(1:nrow(sim.data))){
#       
#       b = append(b,
#                  B_function.star(v, v.star,
#                                  theta = theta.est[-length(gv.true)],
#                                  yi = sim.data$Yi[i],
#                                  zi = as.numeric(Z.df[i,])))
#       
#     }
#     B[v,v.star] = sum(b)
#   }
# }
B = matrix(0, nrow = length(Y_support) - 1, ncol = length(Y_support) - 1)
for(i in c(1:nrow(sim.data))){
  
  B = B + B_function.star(theta = theta.est[-length(Y_support)],
                          yi = sim.data$Yi[i],
                          zi = as.numeric(Z.df[i,]))
  
}

#check
#rownames(B) = c("g1", "g2")
#colnames(B) = c("g1", "g2")
#B

a_data = c()
for(i in c(1:n)){
  a_data = append(a_data,
                  A_tau_function.star(u = sim.data$Ti[i],
                                      d = sim.data$Di[i],
                                      zi = as.numeric(Z.df[i,]),
                                      theta = theta.est[-length(Y_support)]))
}

A = diag(a_data)

#check
#t(X) %*% Z %*% X

# dl.dgv.db = matrix(NA, nrow = m - 1, ncol = length(beta.true))
# #colnames(dl.dgv.db) = c("b0", "b1", "b2")
# #rownames(dl.dgv.db) = c("g1", "g2")
# 
# for(v in c( (Delta + 1):(Delta + m - 1)) ){
#   
#   d_data = c()
#   for(i in c(1:n)){
#     d_data = append(d_data,
#                     D_function.star(v = v,
#                                     zi = as.numeric(Z.df[i,]),
#                                     theta = theta.est[-length(gv.true)]))
#   }
#   
#   D = diag(d_data)
#   
#   dl.dgv.db[(v - Delta), ] = t(J) %*% D %*% Z
#   
# }
D = matrix(NA, nrow = n, ncol = length(Y_support) - 1)
for(i in c(1:n)){
  
  D[i,] = D_function.star(theta = theta.est[-length(Y_support)],
                          zi = as.numeric(Z.df[i,]))
  
}
dl.dgv.db = matrix(NA, nrow = m - 1, ncol = ncol(Z))
for(v in c( (Delta + 1):(Delta + m - 1)) ){
  
  D.star = diag(D[,v - Delta])
  
  dl.dgv.db[(v - Delta), ] = t(J) %*% D.star %*% Z
}

#formulaic hessian
H = cbind(rbind(B, t(dl.dgv.db)),
          rbind(dl.dgv.db, t(Z) %*% A %*% Z))

colnames(H) <- NULL
rownames(H) <- NULL

point.est = theta.est[-length(Y_support)]
st.errs = sqrt(diag(solve(-H)))

point.est + qnorm(0.975) * st.errs
point.est - qnorm(0.975) * st.errs

param = c(paste("g", c((Delta + 1):(m + Delta - 1)), sep = ""),
          colnames(Z.df))

results = data.frame("param" = param,
                     "point.est" = point.est,
                     "std.error" = st.errs)

write.csv(H, "./results/aart-hessian-60mo-robust.csv")
write.csv(results, "./results/aart-est-60mo-robust.csv")

end_time <- Sys.time()
elapsed_time <- end_time - start_time

################################################################################
# summarize results

hessian = read.csv("./results/aart-hessian-60mo-robust.csv")
hessian = hessian[,-1]
aart.est = read.csv("./results/aart-est-60mo-robust.csv")
aart.est = aart.est[,-1]

aart.est$wald.test = aart.est$point.est / aart.est$std.error
aart.est$p.value = 2 * pnorm( abs(aart.est$wald.test), lower.tail = FALSE)

aart.est$sig.code =
  ifelse( (aart.est$p.value <= 0.001), "***",
          ifelse( (aart.est$p.value <= 0.01) & (aart.est$p.value > 0.001), "**",
                  ifelse( (aart.est$p.value <= 0.05) & (aart.est$p.value > 0.01), "*",
                          ifelse( (aart.est$p.value <= 0.10) & (aart.est$p.value > 0.05), ".", ""))))

write.csv(aart.est, "./results/aart-est-60mo-robust.csv")

#summary
aart.est









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
# "./processed-data/aart-2017-36mo.csv
#
# './code/h-star-simulation-function-beta-parameters.R')
# './code/h.star.param.est.R')
# './code/hessian.std.error.R')
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
# 36 MONTHS
######################################################################################

rm(list=ls())

start_time <- Sys.time()

reg.data = read.csv("./processed-data/aart-2017-36mo.csv")
reg.data = reg.data[,-1]

#trapezoid
omega = max(reg.data$Xi)
Delta = min( min(reg.data$Yi), min(reg.data$Xi) ) - 1
m = max(reg.data$Yi) - Delta

X_support = c( (Delta + 1) : (omega) )
Y_support = c( (Delta + 1) : (Delta + m) )

#density function for lifetime; assume binomial
f_density = function(u, mu){
  
  size = (omega - (Delta + 1) + 1) - 1
  x = u - (Delta + 1)
  
  return( dbinom(x, size, mu) )
  
}

#link function; assume exponential
mu = function(eta){
  
  return( exp( eta ) / ( 1 + exp(eta) ) )
  
}

#covariates
eta = function(xi, beta){
  
  return( t(xi) %*% beta )
  
}

#inverse link
inv.link = function(p.est){
  
  return( -log( 1 / p.est - 1) )
}

#for initial value estimation
best.p = function(p){
  
  
  v1 = sapply(X = X_support - 1,
              dbinom,
              size = (omega - (Delta + 1) + 1) - 1,
              prob = p)
  
  return( sum( (v1 - U_est)^2 ) )
  
}

#sample size
n = nrow(reg.data)

#helper functions
source('./code/h-star-simulation-function-beta-parameters.R')
source('./code/h.star.param.est.R')
source('./code/hessian.std.error.R')

################################################################################
# begin regression analysis
################################################################################

sim.data = reg.data


##############################################################################
#parameter estimation
##############################################################################

Z.df = cbind("z0" = rep(1,n), sim.data[,c(3:ncol(sim.data))])
Z = as.matrix(Z.df)
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
                    W_function.star(u = sim.data$Xi[i],
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
                    W_function.star(u = sim.data$Xi[i],
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
#                                  xi = as.numeric(X.df[i,])))
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
                  A_function.star(u = sim.data$Xi[i],
                                  zi = as.numeric(Z.df[i,]),
                                  theta = theta.est[-length(Y_support)]))
}

A = diag(a_data)

#check
#t(X) %*% Z %*% X

#dl.dgv.db = matrix(NA, nrow = m - 1, ncol = length(beta.true))
#colnames(dl.dgv.db) = c("b0", "b1", "b2")
#rownames(dl.dgv.db) = c("g1", "g2")

# for(v in c( (Delta + 1):(Delta + m - 1)) ){
#   
#   a_data = c()
#   for(i in c(1:n)){
#     a_data = append(a_data,
#                     A_function.star(v = v,
#                                     xi = as.numeric(X.df[i,]),
#                                     theta = theta.est[-length(gv.true)]))
#   }
#   
#   A = diag(a_data)
#   
#   dl.dgv.db[(v - Delta), ] = t(J) %*% A %*% X
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

write.csv(H, "./results/aart-hessian-36mo.csv")
write.csv(results, "./results/aart-est-36mo.csv")

end_time <- Sys.time()
elapsed_time <- end_time - start_time

################################################################################
# summarize results

hessian = read.csv("./results/aart-hessian-36mo.csv")
hessian = hessian[,-1]
aart.est = read.csv("./results/aart-est-36mo.csv")
aart.est = aart.est[,-1]

aart.est$wald.test = aart.est$point.est / aart.est$std.error
aart.est$p.value = 2 * pnorm( abs(aart.est$wald.test), lower.tail = FALSE)

aart.est$sig.code =
  ifelse( (aart.est$p.value <= 0.001), "***",
        ifelse( (aart.est$p.value <= 0.01) & (aart.est$p.value > 0.001), "**",
                ifelse( (aart.est$p.value <= 0.05) & (aart.est$p.value > 0.01), "*",
                        ifelse( (aart.est$p.value <= 0.10) & (aart.est$p.value > 0.05), ".", ""))))

write.csv(aart.est, "./results/aart-est-36mo.csv")

#Table 2, 36mo
aart.est

#uniform dist hypothesis test
reg.data = read.csv("./processed-data/aart-2017-36mo.csv")
reg.data = reg.data[,-1]

#trapezoid
omega = max(reg.data$Xi)
Delta = min( min(reg.data$Yi), min(reg.data$Xi) ) - 1
m = max(reg.data$Yi) - Delta

obs_data = reg.data[,1:2]
colnames(obs_data) = c("X", "Y")

#G estimates & hypo test
Cnx <- function(x) {
  ind_x <- ifelse( (obs_data$Y) <= x & (x <= obs_data$X),1,0 )
  return( (1/nrow(obs_data)) * (sum(ind_x)) )
}
bnx <- function(y) {
  return( (1/nrow(obs_data)) * (1/Cnx(y)) * (sum(obs_data$Y == y)) )
}

gnx <- function(x) {
  return(sum(obs_data$Y == x) / nrow(obs_data))
}
cnuv <- function(u,v) {
  return( (sum(ifelse( (obs_data$Y <= min(u,v)) & (obs_data$X >= max(u,v)),1,0))) / nrow(obs_data))
}

y = sort(unique(obs_data$Y))

sapply(y, gnx)

cur_B_vec_var = vector()
for (i in c((Delta + 1):(m+Delta))) {
  cur_B_vec_var = append(cur_B_vec_var, (gnx(i) * cnuv(i-1,i)) / Cnx(i)^3 )
}

len.y = m + Delta - (Delta + 1) + 1

n = nrow(obs_data)
B_H0 = c(rep(1,len.y)) / c(1:len.y)
B_H0 = B_H0[2:len.y]

cur_B_vec = sapply(y, bnx)
cur_B_vec = cur_B_vec[2:len.y]

var_diag <- vector()
for (i in c(1:(len.y - 1))) {
  V = ((cur_B_vec[i])^2 * (1 - cur_B_vec[i])) / gnx(y[i+1])
  var_diag = append(var_diag,V)
}

SIG = diag(var_diag,(len.y - 1),(len.y - 1))

cur_Z_stat = 
  t(sqrt(n) * (cur_B_vec - B_H0) ) %*%
  solve(SIG) %*%
  (sqrt(n) * (cur_B_vec - B_H0) )

pchisq(cur_Z_stat,length(B_H0),lower.tail = FALSE)


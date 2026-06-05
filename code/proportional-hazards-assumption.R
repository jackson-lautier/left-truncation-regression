######################################################################################
######################################################################################
######################################################################################
# Data analysis scripts to reproduce application results in the manuscript:
# 
# "Discrete Time-to-Event Regression Analysis Under Left-Truncation"
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
# INSTRUCTIONS
#
# supporting files:
# "./processed-data/aart-2017-36mo.csv
#
# The code must be run sequentially downwards.
# Results will appear in a newly created 'results' folder.
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

rm(list=ls())

reg.data = read.csv("./processed-data/aart-2017-60mo.csv")
reg.data = reg.data[,-1]

#trapezoid
omega = max(reg.data$Ti)
Delta = min( min(reg.data$Yi), min(reg.data$Ti) ) - 1
m = max(reg.data$Yi) - Delta

X_support = c( (Delta + 1) : (omega) )
Y_support = c( (Delta + 1) : (Delta + m) )

ind.variable = reg.data$new.used

# these estimates follow from
# https://doi.org/10.1016/j.insmatheco.2023.02.003

#IME point estimates + CIs
obs_data <- data.frame(reg.data$Ti, reg.data$Yi, reg.data$Di)
#obs_data <- data.frame(mlease_term$Xc,mlease_term$Y,mlease_term$C)
names(obs_data)[names(obs_data) == 'reg.data.Ti'] <- 'Z'
names(obs_data)[names(obs_data) == 'reg.data.Yi'] <- 'Y'
names(obs_data)[names(obs_data) == 'reg.data.Di'] <- 'C'

obs_data = obs_data[ind.variable == 0,]
n = nrow(obs_data)

#get unique observations of Z
#note: Z = min(X_i, C_i)
z = sort(unique(obs_data$Z))

f_star <- function(x) {
  res = sum( (obs_data$C) * (obs_data$Z == x))
  return(res/n)
}

est_haz <- function(x) {
  num = sum( (obs_data$C) * (obs_data$Z == x))
  den = sum( (x >= obs_data$Y) * (x <= obs_data$Z))
  return(num/den)
}

est_C <- function(x) {
  ans = sum( (x >= obs_data$Y) * (x <= obs_data$Z))
  return(ans/n)
}

lam_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  lam_hat = append(lam_hat,est_haz(i))
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}

Var_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

est_dist = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1
est_dist$S_hat[nrow(est_dist)] = 0

plot_df.0 = data.frame("age" = X_support[4:58],
                       "ind" = as.character(rep("Used",length(X_support[4:58]))),
                       "haz.rate" = est_dist$lam_hat[4:58],
                       "ci.low" = est_dist$CI_lower[4:58],
                       "ci.high" = est_dist$CI_upper[4:58])



#IME point estimates + CIs
obs_data <- data.frame(reg.data$Ti, reg.data$Yi, reg.data$Di)
#obs_data <- data.frame(mlease_term$Xc,mlease_term$Y,mlease_term$C)
names(obs_data)[names(obs_data) == 'reg.data.Ti'] <- 'Z'
names(obs_data)[names(obs_data) == 'reg.data.Yi'] <- 'Y'
names(obs_data)[names(obs_data) == 'reg.data.Di'] <- 'C'

obs_data = obs_data[ind.variable == 1,]
n = nrow(obs_data)

#get unique observations of Z
#note: Z = min(X_i, C_i)
z = sort(unique(obs_data$Z))

f_star <- function(x) {
  res = sum( (obs_data$C) * (obs_data$Z == x))
  return(res/n)
}

est_haz <- function(x) {
  num = sum( (obs_data$C) * (obs_data$Z == x))
  den = sum( (x >= obs_data$Y) * (x <= obs_data$Z))
  return(num/den)
}

est_C <- function(x) {
  ans = sum( (x >= obs_data$Y) * (x <= obs_data$Z))
  return(ans/n)
}

lam_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  lam_hat = append(lam_hat,est_haz(i))
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}

Var_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

est_dist = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1
est_dist$S_hat[nrow(est_dist)] = 0

#prepare plot data
plot_df.1 = data.frame("age" = X_support[4:58],
                       "ind" = as.character(rep("New",length(X_support[4:58]))),
                       "haz.rate" = est_dist$lam_hat[4:58],
                       "ci.low" = est_dist$CI_lower[4:58],
                       "ci.high" = est_dist$CI_upper[4:58])

plot_df = rbind(plot_df.1, plot_df.0)


ggplot() +
  geom_line(data=plot_df, aes(x = age,
                              y = haz.rate,
                              group = ind,
                              colour = ind)) +
  geom_point(data=plot_df, aes(x = age,
                               y = haz.rate,
                               group = ind,
                               colour = ind,
                               shape = ind)) +
  geom_ribbon(data=plot_df, aes(x=age,
                                ymin=ci.low,
                                ymax=ci.high,
                                group = ind,
                                fill = ind), alpha = 0.4) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme_bw() +
  theme(axis.title.x=element_text(size=10, family="Times New Roman"),
        axis.title.y=element_text(size=10,family="Times New Roman"),
        strip.text=element_text(size=10,family="Times New Roman"),
        axis.text=element_text(size=10,family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        #legend.title=element_text(size=10, family="Times New Roman"),
        legend.title=element_blank(),
        legend.position = c(0.1, 0.8))


ggsave("./results/aart-prop-haz.pdf",height=4,width=6,device = cairo_pdf)
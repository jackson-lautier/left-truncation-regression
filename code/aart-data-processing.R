######################################################################################
######################################################################################
######################################################################################
# Data processing scripts to produce the data used in the manuscript:
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

######################################################################################
######################################################################################
######################################################################################
######################################################################################
# INSTRUCTIONS
#
# supporting files:
# ".\raw-data\aart173_compiledr.csv"
#
# ".\code\default_time.R"
#
#The code must be run sequentially downwards.
#As the new, cleaned files are prepared, they will be saved in a new
#folder 'processed-data' in the wd.
#For data analysis, proceed directly to 'data_analysis.R'.
#
#
################################################################################
################################################################################
################################################################################
################################################################################
require('lubridate')


################################################################################
################################################################################
################################################################################
################################################################################

#where processed data will be stored
dir.create('./processed-data/')

################################################################################
# DATA SUMMARIES
################################################################################

rm(list=ls())
path = "./raw-data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

#count
length(unique(aart$assetNumber))

#origination dates
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)

#credit scores
summary(as.numeric(aart$obligorCreditScore))

#loan terms
table(aart$originalLoanTerm)

#original loan balance
summary(as.numeric(aart$reportingPeriodBeginningLoanBalanceAmount))

#manufacturers
length(unique(aart$vehicleManufacturerName))

#manufacturers
length(unique(aart$obligorGeographicLocation))

################################################################################
# DATA SUMMARIES
################################################################################

source("./code/default_time.R")

l.terms = sort( unique(aart$originalLoanTerm) )

#l.terms = c(25, 26, 37, 38, 43, 44, 49, 50, 55, 56, 61, 62, 64, 65, 67, 68, 73,
#            74, 76, 77)

results = matrix(NA, nrow = length(l.terms), ncol = 5)
colnames(results) = c("loan.term", "count", "cens", "repay", "def")
rownames(results) = l.terms

for(l in l.terms){
  
  loan_term_c = l
  len_obs_window = 43 #num. mnths in obs. window
  
  path = "./raw-data/"
  aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))
  
  aart <- aart[aart$originalLoanTerm %in% loan_term_c,]
  aart$originalLoanTerm = l
  
  #calculate remaining payments
  aart_trust_start_date = "06-01-2017"
  date <- paste(aart$originationDate,"-01",sep="")
  date <- as.Date(date, "%m/%Y-%d")
  age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
  aart$initialLoanAge = age
  
  aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge
  
  delta = l - max(aart$remainingTermtoMaturityNumber) 
  M = l - min(aart$remainingTermtoMaturityNumber) - delta
  T_start = M + delta + aart$remainingTermtoMaturityNumber - l 
  Y = M + delta - T_start + 1
  
  ######################################################################
  ######################################################################
  ######################################################################
  # algorithm to find loan outcomes (def, repay, cens)
  ######################################################################
  ######################################################################
  ######################################################################
  X = vector()
  C = vector()
  D = vector()
  R = vector()
  
  for (j in c(1:nrow(aart))) {
    c_bond = default_time(aart[j,])
    X = append(X, c_bond[1])
    C = append(C, c_bond[2])
    R = append(R, c_bond[3])
    D = append(D, c_bond[4])
  }
  
  #shift back to the original timeline
  Xc = M + delta + X - T_start + 1
  
  ######################################################################
  ######################################################################
  aart = cbind(aart,Y,X,Xc,C,D,R)
  
  results[as.character(l), "loan.term"] = l
  results[as.character(l), "count"] = nrow(aart)
  results[as.character(l), "cens"] = sum(aart$C)
  results[as.character(l), "repay"] = sum(aart$R)
  results[as.character(l), "def"] = sum(aart$D)
  
  print(results)
  
}

manu.table = as.data.frame(results)
manu.table$cens.per = manu.table$cens / manu.table$count
manu.table$repay.per = manu.table$repay / manu.table$count
manu.table$def.per = manu.table$def / manu.table$count

sum(results[,"def"]) / sum(results[,"count"])

write.csv(manu.table, "./processed-data/loan-term-table.csv",
          row.names = FALSE)

################################################################################
# 36-MO LOANS
################################################################################

rm(list=ls())

source("./code/default_time.R")

loan_term_c = c(37, 38) 
len_obs_window = 43 #num. mnths in obs. window

path = "./raw-data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

aart <- aart[aart$originalLoanTerm %in% loan_term_c,]
aart$exactLoanTerm = aart$originalLoanTerm
aart$originalLoanTerm = 36

#calculate remaining payments
aart_trust_start_date = "06-01-2017"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = 36 - max(aart$remainingTermtoMaturityNumber) 
M = 36 - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - 36 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #5 loans

aart = aart[!(aart$assetNumber %in% bad_data),]

#check all censored
sum(aart$C); nrow(aart)
aart = aart[aart$C == 0, ]
#check all repayments
sum(aart$R); nrow(aart)
aart = aart[aart$R == 1, ]

table(aart$Xc)
#remove loans with terms greater than 41 months (possible extensions)
aart = aart[aart$Xc <= 41, ]
table(aart$Xc)

#check all delta + 1 <= Y <= M + delta
table(aart$Y)

aart = aart[!(aart$obligorCreditScoreType == "None"),]
aart = aart[!(aart$obligorIncomeVerificationLevelCode == 3),]
#final count: 2,756

obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
table(obs_data$aart.C) #no right-censored observations
#obs_data = obs_data[,1:2]

names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Xi'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Yi'
#names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

omega = max(obs_data$Xi)

#categorical variables check
table(aart$assetSubjectDemandIndicator) #all false
table(aart$coObligorIndicator)
table(aart$obligorCreditScoreType)
table(aart$obligorEmployementVerificationCode) #all NA
table(aart$obligorIncomeVerificationLevelCode)
table(aart$obligorGeographicLocation)
table(aart$originatorName) #all Ally Bank
table(aart$originalInterestRateTypeCode) #all type 1
table(aart$primaryLoanServicerName) #all Ally Bank
table(aart$subvented)
table(aart$vehicleModelName)
table(aart$vehicleNewUsedCode)
table(aart$vehicleManufacturerName)
table(aart$vehicleModelYear)
table(aart$vehicleTypeCode)
table(aart$vehicleValueSourceCode)

#continuous variables to include
summary(as.numeric(aart$obligorCreditScore))
summary(aart$originalInterestRatePercentage)
summary(aart$paymentToIncomePercentage)
summary(log(aart$vehicleValueAmount))
#summary(aart$vehicleModelYear) #treat as continuous

s.df = data.frame("credit.score" = as.numeric(aart$obligorCreditScore),
                  "interest.rate" = aart$originalInterestRatePercentage,
                  "pti" = aart$paymentToIncomePercentage,
                  "veh.value" = log(aart$vehicleValueAmount))
#"veh.year" = aart$vehicleModelYear)

s.df = scale(s.df)
s.df = as.data.frame(s.df)

#categorical variables to include
co.sign = ifelse(aart$coObligorIndicator == "True", 1, 0)
new.used = ifelse(aart$vehicleNewUsedCode == 1, 1, 0) #1 is new

#subvention
table(aart$subvented)

subvent.rate = ifelse(aart$subvented == 1, 1, 0)
subvent.cash = ifelse(aart$subvented == 2, 1, 0)
# 18 - Subvented (Item 3(c)(14))
# If the value in the subvented field is equal to "1," the related receivable is a "subvented receivable." Some receivables are orginated under incentive programs sponsored by vehicle manufacturers, for which the financing rates are below the standard rates at which Ally Bank otherwise offers financing under retail contracts. Those receivables are referred to as subvented receivables.
# Because the rates on the subvented receivables are lower than would otherwise be offered by Ally Bank, Ally Bank purchases those subvented receivables from the dealers and the applicable manufacturer pays the present value of the difference between the customer's subvented rate and Ally Bank's standard rate to Ally Financial on behalf of the dealer selling the related vehicle and Ally Finanical forwards such amount to Ally Bank.
# Subvention is not taken into account by Ally Bank when determining the credit scoring.
# The prospectus does not report cash subvention as a subvented receivable; however, [this filed will show cash subvention, indicated as a "2," if the obligor received cash in connection with the purchase of the related vehicle.

#credit score type
table(aart$obligorCreditScoreType)
#aart = aart[!(aart$obligorCreditScoreType == "None"),]

cred.score.type = ifelse(aart$obligorCreditScoreType == "Commercial Bureau", 1, 0)
#26 - Obligor credit score type (Item 3(e)(1))
#The data in the obligor credit score type field will be "Not Available" or "None" if the
#related obligor is a business without an individual co-obligor and there is no commercial
#bureau score available.  In the case of a business obligor without a co-obligor, a 
#commercial bureau score will be included in this field if it is available.  If the related
#obligor is a business with a co-obligor that has a consumer bureau score, that co-obligor's
#consumer bureau score will be reported.  The obligor credit score type will differ from the
#prospectus because the prospectus classifies the obligor on the receivable as a business
#even if there is an individual co-obligor with a consumer bureau score.


#income verification
table(aart$obligorIncomeVerificationLevelCode)
#aart = aart[!(aart$obligorIncomeVerificationLevelCode == 3),]

inc.verif = ifelse(aart$obligorIncomeVerificationLevelCode == 2, 1, 0)

#28 - Obligor income verification level (Item 3(e)(3))
#In some cases, income will be "not stated, not verified," which is generally as a result
#of the credit application not including the relevant information when it was submitted.
#With respect to a receivable with a business obligor that also has an individual co-obligor,
#the verification level attributable to the individual co-obligor is reported for this field.
#If a receivable with a business obligor does not also have an individual co-obligor, the
#financial information of the related business may be reviewed and will be reported in this field.

cor(inc.verif, cred.score.type)

#exactly correlated; drop cred score type


#vehicle type
table(aart$vehicleTypeCode)

#23 - Vehicle type (Item 3(d)(5))
#The methodology used to populate the vehicle type field conforms with the methodology used
#in creating the tables set forth in "The Receivables Pool" in the prospectus; however, for
#purposes of the prospectus disclosure, sport utility vehicles will be classified as trucks.

head(aart[,c("vehicleTypeCode", "vehicleManufacturerName", "vehicleModelName")])
#1 = sedan
#2 = pick-up truck
#3 = SUV

veh.pick.up = ifelse(aart$vehicleTypeCode == 2, 1, 0)
veh.suv = ifelse(aart$vehicleTypeCode == 3, 1, 0)

#vehicle source
table(aart$vehicleValueSourceCode)

#25 - Source of vehicle value (Item 3(d)(7))
#If the source of vehicle value field is "invoice price," the vehicle value identified
#in field 24 above will be the sum of the manufacturer invoice price and the invoice price
#for any dealer installed accessories. With respect to used vehicles, the vehicle value is
#generally either based on the NADA guide, which will be marked "other," or the Kelley Blue Book.

head(aart[,c("vehicleTypeCode", "vehicleManufacturerName", "vehicleValueSourceCode", "vehicleNewUsedCode")])

used.veh.value.source.1 = ifelse(aart$vehicleValueSourceCode == 3, 1, 0)
used.veh.value.source.2 = ifelse(aart$vehicleValueSourceCode == 98, 1, 0)


#create regression data
reg.data = cbind(obs_data, s.df,
                 co.sign,
                 new.used,
                 subvent.rate,
                 subvent.cash,
                 #cred.score.type,
                 inc.verif,
                 veh.pick.up,
                 veh.suv)
#used.veh.value.source.1,
#used.veh.value.source.2)

reg.data = reg.data[,-3] #no right-censoring

cor_df <- round(cor(reg.data[,3:13]), 2)
library(reshape2)
melted_cor <- melt(cor_df)
library(ggplot2)
ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), size = 5) +
  scale_fill_gradient2(low = "blue", high = "red", limit = c(-1, 1), name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())

#drop income verification indicator
reg.data = reg.data[,-which(colnames(reg.data) == "inc.verif")]

cor_df <- round(cor(reg.data[,3:12]), 2)
#library(reshape2)
melted_cor <- melt(cor_df)
#library(ggplot2)
ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), size = 5) +
  scale_fill_gradient2(low = "blue", high = "red", limit = c(-1, 1), name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())

write.csv(reg.data, './processed-data/aart-2017-36mo.csv')

reg.data$exact.loan.term = aart$exactLoanTerm

write.csv(reg.data, './processed-data/aart-2017-36mo-robust.csv')


################################################################################
# 60 MONTH LOANS
################################################################################

rm(list=ls())

source("./code/default_time.R")

loan_term_c = c(61, 62) 
len_obs_window = 43 #num. mnths in obs. window

path = "./raw-data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

aart <- aart[aart$originalLoanTerm %in% loan_term_c,]
aart$exactLoanTerm = aart$originalLoanTerm
aart$originalLoanTerm = 60

#calculate remaining payments
aart_trust_start_date = "06-01-2017"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = 60 - max(aart$remainingTermtoMaturityNumber) 
M = 60 - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - 60 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #57 loans

aart = aart[!(aart$assetNumber %in% bad_data),]

#check all censored
#sum(aart$C); nrow(aart)
#aart = aart[aart$C == 0, ]
#censoring OK per revision

#remove defaults ~ 3.22%
sum(aart$D); nrow(aart)
aart = aart[aart$D == 0, ]

table(aart$Xc)
#remove loans with terms greater than 68 months (possible extensions)
aart = aart[aart$Xc <= 68, ]
table(aart$Xc)

#check all delta + 1 <= Y <= M + delta
table(aart$Y)
aart = aart[aart$Y <= 64, ]
table(aart$Y)

aart = aart[!(aart$obligorCreditScoreType == "None"),]
aart = aart[!(aart$obligorIncomeVerificationLevelCode == 3),]
#final count: 17,016

#percentage censored
sum(aart$C) / nrow(aart)

obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
obs_data$D = 1 - obs_data$aart.C
#table(obs_data$aart.C) #no right-censored observations
#obs_data = obs_data[,1:2]

names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Ti'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Yi'
names(obs_data)[names(obs_data) == 'D'] <- 'Di'
obs_data = obs_data[,-3]
n = nrow(obs_data)

#omega = max(obs_data$Xi)

#categorical variables check
table(aart$assetSubjectDemandIndicator) #all false
table(aart$coObligorIndicator)
table(aart$obligorCreditScoreType)
table(aart$obligorEmployementVerificationCode) #all NA
table(aart$obligorIncomeVerificationLevelCode)
table(aart$obligorGeographicLocation)
table(aart$originatorName) #all Ally Bank
table(aart$originalInterestRateTypeCode) #all type 1
table(aart$primaryLoanServicerName) #all Ally Bank
table(aart$subvented)
table(aart$vehicleModelName)
table(aart$vehicleNewUsedCode)
table(aart$vehicleManufacturerName)
table(aart$vehicleModelYear)
table(aart$vehicleTypeCode)
table(aart$vehicleValueSourceCode)

#continuous variables to include
summary(as.numeric(aart$obligorCreditScore))
summary(aart$originalInterestRatePercentage)
summary(aart$paymentToIncomePercentage)
summary(log(aart$vehicleValueAmount))
#summary(aart$vehicleModelYear) #treat as continuous

s.df = data.frame("credit.score" = as.numeric(aart$obligorCreditScore),
                  "interest.rate" = aart$originalInterestRatePercentage,
                  "pti" = aart$paymentToIncomePercentage,
                  "veh.value" = log(aart$vehicleValueAmount))
#"veh.year" = aart$vehicleModelYear)

s.df = scale(s.df)

#save for example calcs in supplement
scale.attr = as.data.frame(attr(s.df, "scaled:scale"))
scale.attr = cbind(scale.attr, attr(s.df, "scaled:center"))
write.csv(scale.attr, "./processed-data/aart-60mo-scale-attr.csv")

s.df = as.data.frame(s.df)

#categorical variables to include
co.sign = ifelse(aart$coObligorIndicator == "True", 1, 0)
new.used = ifelse(aart$vehicleNewUsedCode == 1, 1, 0) #1 is new

#subvention
table(aart$subvented)

subvent.rate = ifelse(aart$subvented == 1, 1, 0)
subvent.cash = ifelse(aart$subvented == 2, 1, 0)
# 18 - Subvented (Item 3(c)(14))
# If the value in the subvented field is equal to "1," the related receivable is a "subvented receivable." Some receivables are orginated under incentive programs sponsored by vehicle manufacturers, for which the financing rates are below the standard rates at which Ally Bank otherwise offers financing under retail contracts. Those receivables are referred to as subvented receivables.
# Because the rates on the subvented receivables are lower than would otherwise be offered by Ally Bank, Ally Bank purchases those subvented receivables from the dealers and the applicable manufacturer pays the present value of the difference between the customer's subvented rate and Ally Bank's standard rate to Ally Financial on behalf of the dealer selling the related vehicle and Ally Finanical forwards such amount to Ally Bank.
# Subvention is not taken into account by Ally Bank when determining the credit scoring.
# The prospectus does not report cash subvention as a subvented receivable; however, [this filed will show cash subvention, indicated as a "2," if the obligor received cash in connection with the purchase of the related vehicle.

#credit score type
table(aart$obligorCreditScoreType)
#aart = aart[!(aart$obligorCreditScoreType == "None"),]

cred.score.type = ifelse(aart$obligorCreditScoreType == "Commercial Bureau", 1, 0)
#26 - Obligor credit score type (Item 3(e)(1))
#The data in the obligor credit score type field will be "Not Available" or "None" if the
#related obligor is a business without an individual co-obligor and there is no commercial
#bureau score available.  In the case of a business obligor without a co-obligor, a 
#commercial bureau score will be included in this field if it is available.  If the related
#obligor is a business with a co-obligor that has a consumer bureau score, that co-obligor's
#consumer bureau score will be reported.  The obligor credit score type will differ from the
#prospectus because the prospectus classifies the obligor on the receivable as a business
#even if there is an individual co-obligor with a consumer bureau score.


#income verification
table(aart$obligorIncomeVerificationLevelCode)
#aart = aart[!(aart$obligorIncomeVerificationLevelCode == 3),]

inc.verif = ifelse(aart$obligorIncomeVerificationLevelCode == 2, 1, 0)

#28 - Obligor income verification level (Item 3(e)(3))
#In some cases, income will be "not stated, not verified," which is generally as a result
#of the credit application not including the relevant information when it was submitted.
#With respect to a receivable with a business obligor that also has an individual co-obligor,
#the verification level attributable to the individual co-obligor is reported for this field.
#If a receivable with a business obligor does not also have an individual co-obligor, the
#financial information of the related business may be reviewed and will be reported in this field.

cor(inc.verif, cred.score.type)

#almost exactly correlated; drop cred score type


#vehicle type
table(aart$vehicleTypeCode)

#23 - Vehicle type (Item 3(d)(5))
#The methodology used to populate the vehicle type field conforms with the methodology used
#in creating the tables set forth in "The Receivables Pool" in the prospectus; however, for
#purposes of the prospectus disclosure, sport utility vehicles will be classified as trucks.

head(aart[,c("vehicleTypeCode", "vehicleManufacturerName", "vehicleModelName")])
#1 = sedan
#2 = pick-up truck
#3 = SUV

veh.pick.up = ifelse(aart$vehicleTypeCode == 2, 1, 0)
veh.suv = ifelse(aart$vehicleTypeCode == 3, 1, 0)

#vehicle source
table(aart$vehicleValueSourceCode)

#25 - Source of vehicle value (Item 3(d)(7))
#If the source of vehicle value field is "invoice price," the vehicle value identified
#in field 24 above will be the sum of the manufacturer invoice price and the invoice price
#for any dealer installed accessories. With respect to used vehicles, the vehicle value is
#generally either based on the NADA guide, which will be marked "other," or the Kelley Blue Book.

head(aart[,c("vehicleTypeCode", "vehicleManufacturerName", "vehicleValueSourceCode", "vehicleNewUsedCode")])

used.veh.value.source.1 = ifelse(aart$vehicleValueSourceCode == 3, 1, 0)
used.veh.value.source.2 = ifelse(aart$vehicleValueSourceCode == 98, 1, 0)


#create regression data
reg.data = cbind(obs_data, s.df,
                 co.sign,
                 new.used,
                 subvent.rate,
                 subvent.cash,
                 #cred.score.type,
                 inc.verif,
                 veh.pick.up,
                 veh.suv)
#used.veh.value.source.1,
#used.veh.value.source.2)

#reg.data = reg.data[,-3] #no right-censoring

cor_df <- round(cor(reg.data[,4:14]), 2)
library(reshape2)
melted_cor <- melt(cor_df)
library(ggplot2)
ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), size = 5) +
  scale_fill_gradient2(low = "blue", high = "red", limit = c(-1, 1), name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())

#drop income verification indicator
reg.data = reg.data[,-which(colnames(reg.data) == "inc.verif")]

cor_df <- round(cor(reg.data[,4:13]), 2)
#library(reshape2)
melted_cor <- melt(cor_df)
#library(ggplot2)
ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), size = 5) +
  scale_fill_gradient2(low = "blue", high = "red", limit = c(-1, 1), name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())

write.csv(reg.data, './processed-data/aart-2017-60mo.csv')

reg.data$exact.loan.term = aart$exactLoanTerm

write.csv(reg.data, './processed-data/aart-2017-60mo-robust.csv')











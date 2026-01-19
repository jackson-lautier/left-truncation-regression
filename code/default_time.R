default_time = function(bond_data) {
  
  #get current bond loan performance
  bal_vec <- c(bond_data$reportingPeriodBeginningLoanBalanceAmount)
  pmt_vec <- c(bond_data$totalActualAmountPaid)
  prc_vec <- c(bond_data$actualPrincipalCollectedAmount)
  for (i in c(1:(len_obs_window))) {
    c_bal = as.numeric(bond_data[1,paste("BAL",i,sep="")])
    c_pmt = as.numeric(bond_data[1,paste("PMT",i,sep="")])
    c_prc = as.numeric(bond_data[1,paste("PRC",i,sep="")])
    bal_vec = append(bal_vec,c_bal)
    pmt_vec = append(pmt_vec,c_pmt)
    prc_vec = append(prc_vec,c_prc)
  }
  bal_vec = as.vector(na.omit(bal_vec))
  pmt_vec = as.vector(na.omit(pmt_vec))
  prc_vec = as.vector(na.omit(prc_vec))
  
  #get the time zero loan balance to check for defaults
  init_bal = 
    min(
      bond_data$reportingPeriodActualEndBalanceAmount,
      bal_vec[1]
    )
  
  paid_princ = sum(prc_vec) + 10
  
  #repayment check: total principal paid
  if( paid_princ >= init_bal){
    D = 0 #does not default
    R = 1 #repaid in full
    C = 0 #not censored
    X = min(which(bal_vec == 0),
            length(bal_vec),
            len_obs_window) #time-of-repayment
  }
  
  #default or censor check: 3 consec. zero payments
  if ( paid_princ < init_bal ) {
    run = rle(pmt_vec)
    idx = which((run[[2]]==0 )&(run[[1]] >= 3))
    
    #censored if statement
    if(length(idx) == 0){
      D = 0 #does not default
      R = 0 #does not complete repayment
      C = 1 #censored
      X = len_obs_window
    }
    
    #default if statement (not all missed payments)
    if( (length(idx) > 0) & (idx[1] > 1)) {
      D = 1 #defaulted
      R = 0 #does not repay
      C = 0 #not censored
      X = sum(run[[1]][1:(idx[1]-1)]) + 1 #first occurence of zero pmt
    }
    
    if( (length(idx) > 0) & (idx[1] == 1)) {
      D = 1 #defaulted
      R = 0 #does not repay
      C = 0 #not censored
      X = 1 #all missed payments
    }
  }
  
  #bad data check
  #if( (length(bal_vec) < len_obs_window) & (C = 1) ){
  #  D = 0
  #  R = 1
  #  C = 0
  #  X = which(bal_vec == 0)
  #}
  
  return(c(X,C,R,D))
}
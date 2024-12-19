triplm = function(Y, Xk, mytype='ls', ridged=0) {
  n = length(Y)
  p = ncol(Xk)
  
  # simple linear regression
  if(is.null(p)) {
    p = 1
  }
  
  # prepend the intercept term
  X = cbind(1, Xk)
  
  # construct betahat estimate according to regression type
  if(mytype == 'ls') {
    betahat = solve(t(X) %*% X) %*% t(X) %*% Y
  }
  else if(mytype == 'ridge') {
    betahat = solve((t(X) %*% X) + (ridged^2 * diag(p + 1))) %*% t(X) %*% Y
  }
  else if(mytype == 'lasso') {
    library(glmnet)
    lamideal = cv.glmnet(Xk, Y, alpha=1)$lambda.min
    M = glmnet(Xk, Y, alpha = 1, lambda = lamideal)
    betahat = coef(M)
  }
  else {
    stop("Invalid 'mytype' argument. Choose 'ls', 'ridge', or 'lasso'.")
  }
  
  # y_hat and residuals
  yhat = X %*% betahat
  resid = Y - yhat
  
  # sum of squares
  SSE = sum(resid^2)
  MSE = SSE / (n - p - 1)
  RMSE = sqrt(SSE / n)
  SST = var(Y) * (n - 1)
  MST = SST / (n - 1)
  SSM = SST - SSE
  MSM = SSM / (p)
  
  # F statistic
  Fstat = MSM / MSE
  pval = pf(Fstat, p, n-p-1, lower.tail=F)
  
  # r-squared
  r2 = 1 - SSE/SST
  r2adj = 1 - MSE/MST
  
  # covariance matrix and standard errors
  sig2 = MSE
  sig_bhat = sig2 * solve(t(X) %*% X)   # covariance matrix
  se_bhat = sqrt(diag(sig_bhat))

  # hat matrix and leverages
  H = X %*% solve(t(X) %*% X) %*% t(X)
  lev = diag(H)
  
  # standardized residuals
  sresid = (resid - 0) / (sqrt(MSE) * sqrt(1 - lev))
  
  # correlation matrix
  corr = cor(cbind(Y, Xk))
  
  # standardized coefficients
  betahat_stand = (betahat - 0) / se_bhat
  
  # standardized data coefficients
  Z = cbind(1, scale(Xk))
  zetahat = solve(t(Z) %*% Z) %*% t(Z) %*% Y
  
  # pseudo determinant of design matrix X
  determ = sqrt(det(t(X) %*% X))
  
  # VIF and MCI
  myvifs = c()
  mymcis = c()
  if(p == 1){
    myvifs = append(myvifs, 1)
    mymcis = append(mymcis, 1)
  }
  else{
    for(k in 1:p){
      Xp = X[, -c(k+1)]
      Yp = X[, k+1]
      
      my_betahat = solve(t(Xp) %*% Xp) %*% t(Xp) %*% Yp
      my_yhat = Xp %*% my_betahat
      my_resid = Yp - my_yhat
      
      my_SSE = sum(my_resid^2)
      my_SST = var(Yp) * (length(Yp) - 1)
      my_r2 = 1 - my_SSE/my_SST
      
      my_vif = 1 / (1 - my_r2)
      myvifs = append(myvifs, my_vif)
      my_mci = sqrt(my_vif)
      mymcis = append(mymcis, my_mci)
    }
  }
  
  # partial F statistic (only single variable removed)
  myFs = c()
  if(p == 1){
    myFs = append(myFs, "N/A")
  }
  else{
    for(k in 1:p){
      Xp = X[, -c(k+1)]

      my_betahat = solve(t(Xp) %*% Xp) %*% t(Xp) %*% Y
      my_yhat = Xp %*% my_betahat
      my_resid = Y - my_yhat
      
      my_SSE = sum(my_resid^2)
      my_Fstat = ((n - p - 1) / (p - (p - 1))) * ((my_SSE - SSE) / SSE)
      myFs = append(myFs, my_Fstat)
    }
  }
  
  # added variable plot slopes
  myAVs = c()
  if(p == 1){
    myAVs = append(myAVs, "N/A")
  }
  else{
    for(k in 1:p){
      X_rem = X[, -c(k+1)]
      X_incl = Xk[, c(k)]
      
      y_beta_hat = solve(t(X_rem) %*% X_rem) %*% t(X_rem) %*% Y
      k_beta_hat = solve(t(X_rem) %*% X_rem) %*% t(X_rem) %*% X_incl
      
      y_hat = X_rem %*% y_beta_hat
      k_hat = X_rem %*% k_beta_hat
      
      yk_H = X_rem %*% solve(t(X_rem) %*% X_rem) %*% t(X_rem)
      yk_lev = diag(yk_H)
      
      y_resid = Y - y_hat
      y_SSE = sum(y_resid^2)
      y_MSE = y_SSE / (n - p - 2)
      y_sresid = (y_resid - 0) / (sqrt(y_MSE) * sqrt(1 - yk_lev))
      
      k_resid = X_incl - k_hat
      k_SSE = sum(k_resid^2)
      k_MSE = k_SSE / (n - p - 2)
      k_sresid = (k_resid - 0) / (sqrt(k_MSE) * sqrt(1 - yk_lev))
      
      incl_vs_y = cor(y_sresid, k_sresid) * (sd(y_sresid) / sd(k_sresid))

      myAVs = append(myAVs, incl_vs_y)
    }
  }
  
  # objects to return
  results = list('sample_size'            = n,
                 'num_predictors'         = p,
                 'beta_hat'               = betahat,
                 'SE_beta_hat'            = se_bhat,
                 'standardized_betahat'   = betahat_stand,
                 'Y_hat'                  = yhat,
                 'residuals'              = resid,
                 'SSE'                    = SSE,
                 'MSE'                    = MSE,
                 'RMSE'                   = RMSE,
                 'SSM'                    = SSM,
                 'MSM'                    = MSM,
                 'SST'                    = SST,
                 'MST'                    = MST,
                 'F_stat'                 = Fstat,
                 'p_val'                  = pval,
                 'r2'                     = r2,
                 'r2_adj'                 = r2adj,
                 'hat_matrix'             = H,
                 'leverages'              = lev,
                 'standardized_resids'    = sresid,
                 'correlation_matrix'     = corr,
                 'zeta_hat'               = zetahat,
                 'pseud_determinant'      = determ,
                 'VIF'                    = myvifs,
                 'MCI'                    = mymcis,
                 'partial_F_stat'         = myFs,
                 'added_variable_slopes'  = myAVs
                 )
  
  return(results)
}
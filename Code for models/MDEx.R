#data format: 
#data <- data.frame(
#  p = c(0.625/n,0.25,0.5, 0.75,1-0.625/n),
#  q = (c(a,q1,m,q3,b))
#  )


#MDEp
#method can be one of ("Nelder-Mead", "BFG", "CG", "SANN","Brent")
normal<-function(data,n,startm,startsd,method="Nelder-Mead"){
  
  N<-n
  start_params <- c(mu = startm, sigma = startsd)
  
  
  objective_function_w <- function(params, q, p, N) {
    mu <- params[1]
    sigma <- params[2]
    
    # Ensure sigma is positive
    if (sigma <= 0) return(Inf)
    
    # Calculate predicted probabilities
    q_pred <- qnorm(p, mean = mu, sd = sigma)
    
    x<-qnorm(p,mu,sigma)
    # Calculate density at the quantile
    f <- dnorm(x, mean = mu, sd = sigma)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }

    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q- q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals
    
  }
  
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = method,
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
} #normal
weibull<-function(data,n,starta,startb,method="Nelder-Mead"){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    q_pred <- qweibull(p, shape = bsb, scale = bsa)
    
    x<-qweibull(p,shape = bsb, scale = bsa)
    # Calculate density at the quantile
    f <- dweibull(x, shape = bsb, scale = bsa)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q- q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = method)
  result_w
} #weibull
lognormal<-function(data,n,startm,startsd,method="Nelder-Mead"){
  
  N<-n
  start_params <- c(mu = startm, sigma = startsd)
  
  
  objective_function_w <- function(params, q, p, N) {
    mu <- params[1]
    sigma <- params[2]
    
    
    # Calculate predicted probabilities
    q_pred <- qlnorm(p, meanlog = mu, sdlog = sigma)
    
    x<-qlnorm(p,meanlog=mu,sdlog=sigma)
    # Calculate density at the quantile
    f <- dlnorm(x, meanlog = mu, sdlog = sigma)
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = method,
    hessian = FALSE
  )
  result_w
}#lognormal
beta<-function(data,n,starta,startb,method="Nelder-Mead"){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    q_pred <- qbeta(p, shape1 = bsa, shape2 = bsb)
    
    x<-qbeta(p,shape1 = bsa,shape2 = bsb)
    # Calculate density at the quantile
    f <- dbeta(x, shape1 = bsa, shape2 = bsb)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals 
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = method,  
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
}#beta
egamma<-function(data,n,starta,startb,method="Nelder-Mead"){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    
    # Calculate predicted probabilities
    q_pred <- qgamma(p, shape = bsa, scale = bsb)
    
    x<-qgamma(p,shape = bsa, scale = bsb)
    # Calculate density at the quantile
    f <- dgamma(x, shape = bsa, scale = bsb)
    
    # Calculate variance and weights
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/( f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights%*% residuals
    weighted_residuals
    
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = method,  # Method that allows for bounds on parameters
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
}#gamma

#L-BFGS-B Method
#lowsd:lower bound of SD; upsd:upper bound of SD
#THe mean will set to the range of [q1,q3]
normal<-function(data,n,startm,startsd,lowsd=1e-3,upsd=50){
  
  N<-n
  start_params <- c(mu = startm, sigma = startsd)
  
  
  objective_function_w <- function(params, q, p, N) {
    mu <- params[1]
    sigma <- params[2]
    
    # Ensure sigma is positive
    if (sigma <= 0) return(Inf)
    
    # Calculate predicted probabilities
    q_pred <- qnorm(p, mean = mu, sd = sigma)
    
    x<-qnorm(p,mu,sigma)
    # Calculate density at the quantile
    f <- dnorm(x, mean = mu, sd = sigma)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
  
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights%*% residuals
    weighted_residuals
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "L-BFGS-B",
    lower = c(data$q[2],lowsd),
    upper = c(data$q[4],upsd),
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
} #normal

#lowa, lowb: lower bound of shape parameter and scale paramter;
#upa, upb:upper bound of shape parameter and scale paramter
weibull<-function(data,n,starta,startb,lowa=1e-3,upa=100,lowb=1e-3,upb=100){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    q_pred <- qweibull(p, shape = bsb, scale = bsa)
    
    x<-qweibull(p,shape = bsb, scale = bsa)
    # Calculate density at the quantile
    f <- dweibull(x, shape = bsb, scale = bsa)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q- q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "L-BFGS-B",
    lower = c(lowa, lowb),
    upper = c(upa, upb))
  
  result_w
} #weibull

#lowsd:lower bound of sigma,upds:upper bound of sigma
#location paramter will be in the range of [log(q2),log(q3)]
lognormal<-function(data,n,startm,startsd,lowsd=1e-3,upsd=10){
  
  N<-n
  start_params <- c(mu = startm, sigma = startsd)
  
  
  objective_function_w <- function(params, q, p, N) {
    mu <- params[1]
    sigma <- params[2]
    
    
    # Calculate predicted probabilities
    # Calculate predicted probabilities
    # Calculate predicted probabilities
    q_pred <- qlnorm(p, meanlog = mu, sdlog = sigma)
    
    x<-qlnorm(p,meanlog=mu,sdlog=sigma)
    # Calculate density at the quantile
    f <- dlnorm(x, meanlog = mu, sdlog = sigma)
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "L-BFGS-B",
    lower = c(log(data$q[2]), lowsd),
    upper = c(log(data$q[4]), upsd),
    hessian = FALSE
  )
  result_w
}#lognormal

#lowa, lowb: lower bound of shape parameters a and b, respectively;
#upa, upb:upper bound of shape parameter and scale paramter
beta<-function(data,n,starta,startb,lowa=1e-3,upa=40,lowb=1e-3,upb=40){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    q_pred <- qbeta(p, shape1 = bsa, shape2 = bsb)
    
    x<-qbeta(p,shape1 = bsa,shape2 = bsb)
    # Calculate density at the quantile
    f <- dbeta(x, shape1 = bsa, shape2 = bsb)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals 
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "L-BFGS-B",  
    lower = c(lowa,lowb),
    upper = c(upb,upb), 
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
}#beta

#lowa, lowb: lower bound of shape parameter and scale paramter, respectively;
#upa, upb:upper bound of shape parameter and scale paramter
egamma<-function(data,n,starta,startb,lowa=1e-3,upa=100,lowb=1e-3,upb=100){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    
    # Calculate predicted probabilities
    q_pred <- qgamma(p, shape = bsa, scale = bsb)
    
    x<-qgamma(p,shape = bsa, scale = bsb)
    # Calculate density at the quantile
    f <- dgamma(x, shape = bsa, scale = bsb)
    
    # Calculate variance and weights
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/( f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- qr.solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- q - q_pred
    weighted_residuals <- t(residuals) %*% weights%*% residuals
    weighted_residuals
    
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "L-BFGS-B",  # Method that allows for bounds on parameters
    lower = c(lowa,lowb),
    upper = c(upb,upb), 
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
}#gamma


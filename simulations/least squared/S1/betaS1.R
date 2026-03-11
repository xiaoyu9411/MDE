
library(stats) 
library(devtools)
devtools::install_github("stmcg/estmeansd",lib="/usr3/graduate/rainie")
library(metaBLUE,lib.loc ="/usr3/graduate/rainie" )
library(estmeansd,lib.loc ="/usr3/graduate/rainie" )
library(pracma)




get.scenario <- function(min.val, q1.val, med.val, q3.val, max.val) {
  if (!(missing(min.val) | missing(q1.val) | missing(med.val) | missing(q3.val)
        | missing(max.val))) {
    if (!any(is.na(c(min.val, q1.val, med.val, q3.val, max.val)))) {
      return("S3")
    }
  }
  if (!(missing(q1.val) | missing(med.val) | missing(q3.val))) {
    if (!any(is.na(c(q1.val, med.val, q3.val)))) {
      return("S2")
    }
  }
  if (!(missing(min.val) | missing(med.val) | missing(max.val))) {
    if (!any(is.na(c(min.val, med.val, max.val)))) {
      return("S1")
    }
  }
  stop("Summary measures not in appropriate form. See documentation for
       appropriate forms.")
}
set.qe.fit.control <- function(quants, n, scenario){
  con <- list()
  if (scenario == "S1" | scenario == "S2") {
    con$norm.mu.bounds <- c(quants[1], quants[3])
    med.val <- quants[2]
    if (min(quants) > 0) {
      con$lnorm.mu.start <- log(med.val)
      con$lnorm.mu.bounds <- c(log(quants[1]), log(quants[3]))
    }
  }
  if (scenario == "S3") {
    con$norm.mu.bounds <- c(quants[2], quants[4])
    med.val <- quants[3]
    if (min(quants) > 0) {
      con$lnorm.mu.start <- log(med.val)
      con$lnorm.mu.bounds <- c(log(quants[2]), log(quants[4]))
    }
  }
  
  con$norm.sigma.bounds <- c(1e-3, 50)
  con$lnorm.sigma.bounds = c(1e-3, 10)
  
  mean.hat <- metaBLUE::Luo.mean(quants, n, scenario)$muhat
  sd.hat <- metaBLUE::Wan.std(quants, n, scenario)$sigmahat
  
  con$norm.mu.start <- mean.hat
  con$norm.sigma.start <- sd.hat
  
  con$lnorm.mu.start <- log(mean.hat / sqrt(1 + (sd.hat/mean.hat)^2))
  con$lnorm.sigma.start <- sqrt(log( 1 + (sd.hat/mean.hat)^2))
  
  con$gamma.shape.start <- mean.hat^2/sd.hat^2
  con$gamma.rate.start <- mean.hat/sd.hat^2
  
  con$beta.shape1.start <- mean.hat *
    (((mean.hat * (1 - mean.hat)) / (sd.hat^2)) - 1)
  con$beta.shape2.start <- con$beta.shape1.start * (1 - mean.hat) / mean.hat
  
  start.val <- (mean.hat/sd.hat)^1.086
  obj.fun <- function(k, mean.hat, sd.hat) {
    temp1 <- gamma((k + 2)/k)
    temp2 <- gamma((k + 1)/k)
    temp3 <- sqrt((temp1 / (temp2^2)) - 1)
    ((sd.hat / mean.hat) - temp3)^2
  }
  no.fit <- function(e) {
    return(list(par = NA))
  }
  get.weibull.start <- tryCatch({
    stats::nlminb(start = start.val, objective = obj.fun, lower = 1e-3,
                  mean.hat = mean.hat, sd.hat = sd.hat)
  },
  error = no.fit,
  warning = no.fit
  )
  if (!is.na(get.weibull.start$par)){
    con$weibull.shape.start <- get.weibull.start$par
  } else {
    con$weibull.shape.start <- start.val
  }
  con$weibull.scale.start <- mean.hat / gamma((con$weibull.shape.start + 1)
                                              / con$weibull.shape.start)
  
  con$gamma.shape.bounds <- c(1e-3, 1e2)
  con$gamma.rate.bounds <- c(1e-3, 1e2)
  con$weibull.shape.bounds <- c(1e-3, 1e2)
  con$weibull.scale.bounds <- c(1e-3, 1e2)
  con$beta.shape1.bounds <- c(10^(-3), 40)
  con$beta.shape2.bounds <- c(10^(-3), 40)
  return(con)
}
get.mean.sd <- function(x, family) {
  if (!(family %in% c("normal", "log-normal", "gamma", "weibull", "beta"))) {
    stop("family must be either normal, log-normal, gamma, Weibull, or beta.")
  }
  if (family == "normal") {
    par <- unname(x$norm.par)
    est.mean <- par[1]
    est.sd <- par[2]
  }
  else if (family == "log-normal") {
    par <- unname(x$lnorm.par)
    est.mean <- exp(par[1] + par[2]^2 / 2)
    est.sd <- sqrt((exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2))
  }
  else if (family == "gamma") {
    par <- unname(x$gamma.par)
    est.mean <- par[1] / par[2]
    est.sd <- sqrt(par[1] / (par[2]^2))
  }
  else if (family == "weibull") {
    par <- unname(x$weibull.par)
    est.mean <- par[2]* gamma(1 + 1/par[1])
    est.sd <- sqrt(par[2]^2 * (gamma(1 + 2 / par[1]) -
                                 (gamma(1 + 1 / par[1]))^2))
  }
  else if (family == "beta") {
    par <- unname(x$beta.par)
    est.mean <- par[1]/(par[1]+par[2])
    est.sd <- sqrt(par[1] * par[2] / ((par[1] + par[2])^2 *
                                        (par[1] + par[2] + 1)))
  }
  return(list(est.mean = est.mean, est.sd = est.sd))
}
get.num.input <- function(min.val, q1.val, med.val, q3.val, max.val, n){
  res <- list()
  if (!missing(min.val)){
    res$min.val <- min.val
  }
  if (!missing(q1.val)){
    res$q1.val <- q1.val
  }
  if (!missing(med.val)){
    res$med.val <- med.val
  }
  if (!missing(q3.val)){
    res$q3.val <- q3.val
  }
  if (!missing(max.val)){
    res$max.val <- max.val
  }
  if (!missing(n)){
    res$n <- n
  }
  return(res)
}

#weighted QE
wqe.fit<-function(min.val, q1.val, med.val, q3.val, max.val, n)
{
  scenario <- get.scenario(min.val, q1.val, med.val, q3.val, max.val)
  
  if (scenario == "S1") {
    probs <- c(0.625/ n, 0.5, 1 - 0.625 / n)
    quants <- c(min.val, med.val, max.val)
  } else if (scenario == "S2") {
    probs <- c(0.25, 0.5, 0.75)
    quants <- c(q1.val, med.val, q3.val)
  } else if (scenario == "S3") {
    probs <- c(0.625/ n, 0.25, 0.5, 0.75, 1 - 0.625/ n)
    quants <- c(min.val, q1.val, med.val, q3.val, max.val)}
  
  if (min(quants == 0)) {
    quants[quants == 0] <- 10^(-2)
  }
  
  con <-  set.qe.fit.control(quants, n, scenario)
  no.fit <- function(e) {
    return(list(par = NA, value = NA))
  }
  
  fit.norm <- tryCatch({
    N<-n
    start_params <- c(mu = con$norm.mu.start,sigma = con$norm.sigma.start)
    p<-probs
    q<-quants 
    
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
          if(i==j){
            cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
            cov[j,i]<-cov[i,j]}
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$norm.mu.bound[1], con$norm.sigma.bound[1]),
      upper = c(con$norm.mu.bound[2], con$norm.sigma.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.n<-qnorm(probs, mean = as.numeric(fit.norm$par[1]), sd = as.numeric(fit.norm$par[2]))
  residual.n<-(quants -pre.n)^2
  value.n<-sum(residual.n)
  
  fit.lnorm <- tryCatch({
    N<-n
    start_params <- c(mu = con$lnorm.mu.start,sigma = con$lnorm.sigma.start)
    p<-probs
    q<-quants 
    
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
          if(i==j){
            cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
            cov[j,i]<-cov[i,j]}
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$lnorm.mu.bound[1], con$lnorm.sigma.bound[1]),
      upper = c(con$lnorm.mu.bound[2], con$lnorm.sigma.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.l<-qlnorm(probs, meanlog = as.numeric(fit.lnorm$par[1]), sdlog = as.numeric(fit.lnorm$par[2]))
  residual.l<-(quants -pre.l)^2
  value.l<-sum(residual.l)
  
  
  fit.weibull <- tryCatch({
    N<-n
    start_params <- c(bsa = con$weibull.shape.start, bsb = con$weibull.scale.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      q_pred <- qweibull(p, shape = bsa, scale = bsb)
      
      x<-qweibull(p,shape = bsa, scale = bsb)
      # Calculate density at the quantile
      f <- dweibull(x, shape = bsa, scale = bsb)
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights
      for (i in 1:length(p)){
        for (j in i:length(p)){
          if(i==j){
            cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
            cov[j,i]<-cov[i,j]}
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$weibull.shape.bound[1],con$weibull.scale.bound[1]),
      upper = c(con$weibull.shape.bound[2],con$weibull.scale.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.w<-qweibull(probs,shape = as.numeric(fit.weibull$par[1]), scale = as.numeric(fit.weibull$par[2]))
  residual.w<-(quants -pre.w)^2
  value.w<-sum(residual.w)
  
  fit.beta <- tryCatch({
    N<-n
    start_params <- c(con$beta.shape1.start, con$beta.shape2.start)
    p<-probs
    q<-quants 
    
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
          if(i==j){
            cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
            cov[j,i]<-cov[i,j]}
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$beta.shape1.bound[1],con$beta.shape2.bound[1]),
      upper = c(con$beta.shape1.bound[2],con$beta.shape2.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.b<-qbeta(probs,shape1 = as.numeric(fit.beta$par[1]), shape2 = as.numeric(fit.beta$par[2]))
  residual.b<-(quants -pre.b)^2
  value.b<-sum(residual.b)
  
  
  fit.gamma <- tryCatch({
    N<-n
    start_params <- c(con$gamma.shape.start, con$gamma.rate.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      q_pred <- qgamma(p, shape = bsa, rate = bsb)
      
      x<-qgamma(p,shape = bsa, rate = bsb)
      # Calculate density at the quantile
      f <- dgamma(x, shape = bsa, rate = bsb)
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights
      for (i in 1:length(p)){
        for (j in i:length(p)){
          if(i==j){
            cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
            cov[j,i]<-cov[i,j]}
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$gamma.shape.bound[1], con$gamma.rate.bound[1]),
      upper = c(con$gamma.shape.bound[2], con$gamma.rate.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.g<-qgamma(probs,shape = as.numeric(fit.gamma$par[1]), rate = as.numeric(fit.gamma$par[2]))
  residual.g<-(quants -pre.g)^2
  value.g<-sum(residual.g)
  
  values <- c(value.n, value.l, value.g,
              value.w, value.b)
  
  names(values) <- c("normal", "log-normal", "gamma", "weibull", "beta")
  norm.par <- fit.norm$par
  lnorm.par <- fit.lnorm$par
  gamma.par <- fit.gamma$par
  weibull.par <- fit.weibull$par
  beta.par <- fit.beta$par
  num.input <- get.num.input(min.val, q1.val, med.val, q3.val, max.val, n)
  output <- list(norm.par = norm.par, lnorm.par = lnorm.par,
                 gamma.par = gamma.par, weibull.par = weibull.par,
                 beta.par = beta.par, values = values,
                 num.input = num.input, scenario = scenario)
  class(output) <- "wqe.fit"
  return(output)
  
}
wqe.mean.sd <- function(min.val, q1.val, med.val, q3.val, max.val, n) {
  args <- as.list(environment())
  x <- wqe.fit(min.val = min.val, q1.val = q1.val, med.val = med.val,
               q3.val = q3.val, max.val = max.val, n = n)
  selected.dist <- names(which.min(x$values))
  ests <- get.mean.sd(x, selected.dist)
  output <- list(est.mean = ests$est.mean, est.sd = ests$est.sd,
                 selected.dist = selected.dist, values = x$values,args=args,
                 scenario = x$scenario, fitted.dists = x)
  for (dist.name in names(x$values)){
    ests.all <- get.mean.sd(x, dist.name)
    output[paste0(dist.name, '.est.mean')] <- ests.all$est.mean
    output[paste0(dist.name, '.est.sd')] <- ests.all$est.sd
  }
  class(output) <- "wqe.mean.sd"
  return(output)
}

#MDEx

MDEx.fit<-function(min.val, q1.val, med.val, q3.val, max.val, n)
{
  scenario <- get.scenario(min.val, q1.val, med.val, q3.val, max.val)
  
  if (scenario == "S1") {
    probs <- c(0.625 / n, 0.5, 1 - 0.625 / n)
    quants <- c(min.val, med.val, max.val)
  } else if (scenario == "S2") {
    probs <- c(0.25, 0.5, 0.75)
    quants <- c(q1.val, med.val, q3.val)
  } else if (scenario == "S3") {
    probs <- c(0.625 / n, 0.25, 0.5, 0.75, 1 - 0.625/ n)
    quants <- c(min.val, q1.val, med.val, q3.val, max.val)}
  
  if (min(quants == 0)) {
    quants[quants == 0] <- 10^(-2)
  }
  
  con <-  set.qe.fit.control(quants, n, scenario)
  no.fit <- function(e) {
    return(list(par = NA, value = NA))
  }
  
  
  fit.norm <- tryCatch({
    N<-n
    start_params <- c(mu = con$norm.mu.start,sigma = con$norm.sigma.start)
    p<-probs
    q<-quants 
    
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
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$norm.mu.bound[1], con$norm.sigma.bound[1]),
      upper = c(con$norm.mu.bound[2], con$norm.sigma.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.n<-qnorm(probs, mean = as.numeric(fit.norm$par[1]), sd = as.numeric(fit.norm$par[2]))
  residual.n<-(quants -pre.n)^2
  value.n<-sum(residual.n)
  
  fit.lnorm <- tryCatch({
    N<-n
    start_params <- c(mu = con$lnorm.mu.start,sigma = con$lnorm.sigma.start)
    p<-probs
    q<-quants 
    
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
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$lnorm.mu.bound[1], con$lnorm.sigma.bound[1]),
      upper = c(con$lnorm.mu.bound[2], con$lnorm.sigma.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.l<-qlnorm(probs, meanlog = as.numeric(fit.lnorm$par[1]), sdlog = as.numeric(fit.lnorm$par[2]))
  residual.l<-(quants -pre.l)^2
  value.l<-sum(residual.l)
  
  
  fit.weibull <- tryCatch({
    N<-n
    start_params <- c(bsa = con$weibull.shape.start, bsb = con$weibull.scale.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      q_pred <- qweibull(p, shape = bsa, scale = bsb)
      
      x<-qweibull(p,shape = bsa, scale = bsb)
      # Calculate density at the quantile
      f <- dweibull(x, shape = bsa, scale = bsb)
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
          cov[j,i]<-cov[i,j]
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$weibull.shape.bound[1],con$weibull.scale.bound[1]),
      upper = c(con$weibull.shape.bound[2],con$weibull.scale.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.w<-qweibull(probs,shape = as.numeric(fit.weibull$par[1]), scale = as.numeric(fit.weibull$par[2]))
  residual.w<-(quants -pre.w)^2
  value.w<-sum(residual.w)
  
  fit.beta <- tryCatch({
    N<-n
    start_params <- c(con$beta.shape1.start, con$beta.shape2.start)
    p<-probs
    q<-quants 
    
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
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$beta.shape1.bound[1],con$beta.shape2.bound[1]),
      upper = c(con$beta.shape1.bound[2],con$beta.shape2.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.b<-qbeta(probs,shape1 = as.numeric(fit.beta$par[1]), shape2 = as.numeric(fit.beta$par[2]))
  residual.b<-(quants -pre.b)^2
  value.b<-sum(residual.b)
  
  
  fit.gamma <- tryCatch({
    N<-n
    start_params <- c(con$gamma.shape.start, con$gamma.rate.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      q_pred <- qgamma(p, shape = bsa, rate = bsb)
      
      x<-qgamma(p,shape = bsa, rate = bsb)
      # Calculate density at the quantile
      f <- dgamma(x, shape = bsa, rate = bsb)
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p[i]*(1-p[j])/(f[i]*f[j])
          cov[j,i]<-cov[i,j]
        }
      }
      weights <- qr.solve(cov)
      # Calculate weighted sum of squared residuals
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights%*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$gamma.shape.bound[1], con$gamma.rate.bound[1]),
      upper = c(con$gamma.shape.bound[2], con$gamma.rate.bound[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.g<-qgamma(probs,shape = as.numeric(fit.gamma$par[1]), rate = as.numeric(fit.gamma$par[2]))
  residual.g<-(quants -pre.g)^2
  value.g<-sum(residual.g)
  
  values <- c(value.n, value.l, value.g,
              value.w, value.b)
  
  names(values) <- c("normal", "log-normal", "gamma", "weibull", "beta")
  norm.par <- fit.norm$par
  lnorm.par <- fit.lnorm$par
  gamma.par <- fit.gamma$par
  weibull.par <- fit.weibull$par
  beta.par <- fit.beta$par
  num.input <- get.num.input(min.val, q1.val, med.val, q3.val, max.val, n)
  output <- list(norm.par = norm.par, lnorm.par = lnorm.par,
                 gamma.par = gamma.par, weibull.par = weibull.par,
                 beta.par = beta.par, values = values,
                 num.input = num.input, scenario = scenario)
  class(output) <- "MDE.fit"
  return(output)
  
}
MDEx.mean.sd <- function(min.val, q1.val, med.val, q3.val, max.val, n) {
  args <- as.list(environment())
  x <- MDEx.fit(min.val = min.val, q1.val = q1.val, med.val = med.val,
                q3.val = q3.val, max.val = max.val, n = n)
  selected.dist <- names(which.min(x$values))
  ests <- get.mean.sd(x, selected.dist)
  output <- list(est.mean = ests$est.mean, est.sd = ests$est.sd,
                 selected.dist = selected.dist, values = x$values,args=args,
                 scenario = x$scenario, fitted.dists = x)
  for (dist.name in names(x$values)){
    ests.all <- get.mean.sd(x, dist.name)
    output[paste0(dist.name, '.est.mean')] <- ests.all$est.mean
    output[paste0(dist.name, '.est.sd')] <- ests.all$est.sd
  }
  class(output) <- "MDEx.mean.sd"
  return(output)
}




set.seed(1456)


#lognormal distribution
sa<-4
sb<-2
truemean<-8*sa/(sb+sa)
truesd<-8*sqrt((sa*sb)/((sa+sb)^2*(sa+sb+1)))
nnum<-seq(from=30,to=1000,by=10)
nsim<-10000
it<-0
rbmluo<-c()
rbmqe<-c()
rbmwqe<-c()
rbmmdex<-c()

rbsdluo<-c()
rbsdqe<-c()
rbsdwqe<-c()
rbsdmdex<-c()

msemluo<-c()
msemqe<-c()
msemwqe<-c()
msemmdex<-c()
msemsample<-c()

msesdluo<-c()
msesdqe<-c()
msesdwqe<-c()
msesdmdex<-c()
msesdsample<-c()

pqe<-c()
pwqe<-c()
pmdex<-c()

disqe<-c()
diswqe<-c()
dismdex<-c()


it<-0
xbar<-c()
sd<-c()

summarydat<-matrix(nrow=length(nnum),ncol=25)




for (i in nnum ){
  
  mluo<-c()
  sdshi<-c()
  mqe<-c()
  sdqe<-c()
  xbar<-c()
  sd<-c()
  mwqe<-c()
  sdwqe<-c()
  mmdex<-c()
  sdmdex<-c()
  disqeall<-c()
  diswqeall<-c()
  dismdexall<-c()
  
  for (j in 1:nsim){
    n<-i
    datset<-8*rbeta(n,shape1=sa,shape2=sb)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])/8
    q1<-as.numeric(s[2])/8
    m<-as.numeric(s[3])/8
    q3<-as.numeric(s[5])/8
    b<-as.numeric(s[6])/8
    
    #luo's method
    mluo[j]<-8*(((4/(4+n^0.75))*((a+b)/2))+((n^0.75)/(4+n^0.75))*m)
    #Shi's method
    theta1<-2*qnorm((n-0.375)/(n+0.25))
    sdshi[j]<-8*((b-a)/theta1)
    
    #qe
    ss<-qe.mean.sd(min.val = a,med.val = m,max.val = b,n=n)
    mqe[j]<-8*ss$est.mean
    sdqe[j]<-8*ss$est.sd
    disqeall[j]<-ss$selected.dist
    #wqe
    wqe<-wqe.mean.sd(min.val = a,med.val = m,max.val = b,n=n)
    mwqe[j]<-8*wqe$est.mean
    sdwqe[j]<-8*wqe$est.sd
    diswqeall[j]<-wqe$selected.dist
    
    #mdex
    MDEx<-MDEx.mean.sd(min.val = a,med.val = m,max.val = b,n=n)
    mmdex[j]<-8*MDEx$est.mean
    sdmdex[j]<-8*MDEx$est.sd
    dismdexall[j]<-MDEx$selected.dist
  }
  
  it<-it+1
  
  rbmluo[it]<-mean((mluo-truemean)/truemean)
  rbmqe[it]<-mean((mqe[(!is.na(mqe))&mqe<100]-truemean)/truemean)
  rbmwqe[it]<-mean((mwqe[(!is.na(mwqe))&mwqe<100]-truemean)/truemean)
  rbmmdex[it]<-mean((mmdex[(!is.na(mmdex))&mmdex<100]-truemean)/truemean)
  
  rbsdluo[it]<-mean((sdshi-truesd)/truesd)
  rbsdqe[it]<-mean((sdqe[(!is.na(sdqe))&sdqe<100]-truesd)/truesd)
  rbsdwqe[it]<-mean((sdwqe[(!is.na(sdwqe))&sdwqe<100]-truesd)/truesd)
  rbsdmdex[it]<-mean((sdmdex[(!is.na(sdmdex))&sdmdex<100]-truesd)/truesd)
  
  msemluo[it]<-mean((mluo-truemean)^2)
  msemqe[it]<-mean((mqe[(!is.na(mqe))&mqe<100]-truemean)^2)
  msemwqe[it]<-mean((mwqe[(!is.na(mwqe))&mwqe<100]-truemean)^2)
  msemmdex[it]<-mean((mmdex[(!is.na(mmdex))&mmdex<100]-truemean)^2)
  msemsample[it]<-mean((xbar-truemean)^2)
  
  msesdluo[it]<-mean((sdshi-truesd)^2)
  msesdqe[it]<-mean((sdqe[(!is.na(sdqe))&sdqe<100]-truesd)^2)
  msesdwqe[it]<-mean((sdwqe[(!is.na(sdwqe))&sdwqe<100]-truesd)^2)
  msesdmdex[it]<-mean((sdmdex[(!is.na(sdmdex))&sdmdex<100]-truesd)^2)
  msesdsample[it]<-mean((sd-truesd)^2)
  
  pqe[it]<-mean(is.na(mqe))
  pwqe[it]<-mean(is.na(mwqe)|mwqe>100)
  pmdex[it]<-mean(is.na(mmdex)|mmdex>100)
  
  disqe[it]<-mean(disqeall=="beta")
  diswqe[it]<-mean(diswqeall=="beta")
  dismdex[it]<-mean(dismdexall=="beta")
  
  
  summarydat[it,]<-c(i,rbmluo[it],rbsdluo[it],rbmqe[it],rbsdqe[it],
                     rbmwqe[it],rbsdwqe[it],rbmmdex[it],rbsdmdex[it],
                     msemluo[it],msesdluo[it],msemqe[it],msesdqe[it],
                     msemwqe[it],msesdwqe[it],msemmdex[it],msesdmdex[it],
                     msemsample[it],msesdsample[it],
                     pqe[it],pwqe[it],pmdex[it],
                     disqe[it],diswqe[it],dismdex[it])
  write.table(summarydat[it,],file="resultbetaS1.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}

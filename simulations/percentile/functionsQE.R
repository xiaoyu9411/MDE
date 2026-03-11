library(stats) 
library(metaBLUE )
library(estmeansd)
library(pracma)

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
set.qe.fit.control <- function(quants, n, scenario, twosample_default){
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
  
  if (twosample_default){
    con$norm.mu.start = med.val
    con$norm.sigma.start = 1
    con$lnorm.sigma.start = 1
    con$gamma.shape.start = 1
    con$gamma.rate.start = 1
    con$gamma.shape.bounds = c(1e-3, 40)
    con$gamma.rate.bounds = c(1e-3, 40)
    con$weibull.shape.start = 1
    con$weibull.scale.start = 1
    con$weibull.shape.bounds = c(1e-3, 50)
    con$weibull.scale.bounds = c(1e-3, 50)
  } else {
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
  }
  return(con)
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
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep="", collapse=" ")
}
qe.fit <- function(min.val, q1.val, med.val, q3.val, max.val, n,p,
                   two.sample.default = FALSE, qe.fit.control = list()) {
  
  scenario <- get.scenario(min.val, q1.val, med.val, q3.val, max.val)
  check_errors(min.val = min.val, q1.val = q1.val, med.val = med.val,
               q3.val = q3.val, max.val = max.val, n = n, scenario = scenario)
  
  if (scenario == "S1") {
    probs <- c(p / n, 0.5, 1 - p / n)
    quants <- c(min.val, med.val, max.val)
  } else if (scenario == "S2") {
    probs <- c(0.25, 0.5, 0.75)
    quants <- c(q1.val, med.val, q3.val)
  } else if (scenario == "S3") {
    probs <- c(p / n, 0.25, 0.5, 0.75, 1 - p / n)
    quants <- c(min.val, q1.val, med.val, q3.val, max.val)
  }
  
  if (min(quants == 0)) {
    quants[quants == 0] <- 10^(-2)
  }
  
  con <-  set.qe.fit.control(quants, n, scenario, two.sample.default)
  con[names(qe.fit.control)] <- qe.fit.control
  
  S.theta.norm <- function(theta) {
    summand <- sum((stats::qnorm(p = probs, mean = theta[1],
                                 sd = theta[2]) - quants)^2)
  }
  S.theta.lnorm <- function(theta) {
    summand <- sum((stats::qlnorm(p = probs, meanlog = theta[1],
                                  sdlog = theta[2]) - quants)^2)
  }
  S.theta.gamma <- function(theta) {
    summand <- sum((stats::qgamma(p = probs, shape = theta[1],
                                  rate = theta[2]) - quants)^2)
  }
  S.theta.weibull <- function(theta) {
    summand <- sum((stats::qweibull(p = probs, shape = theta[1],
                                    scale = theta[2]) - quants)^2)
  }
  S.theta.beta <- function(theta) {
    summand <- sum((stats::qbeta(p = probs, shape1 = theta[1],
                                 shape2 = theta[2]) - quants)^2)
  }
  
  no.fit <- function(e) {
    return(list(par = NA, value = NA))
  }
  
  fit.norm <- tryCatch({
    stats::optim(par = c(con$norm.mu.start, con$norm.sigma.start),
                 S.theta.norm, method = "L-BFGS-B",
                 lower = c(con$norm.mu.bounds[1], con$norm.sigma.bounds[1]),
                 upper = c(con$norm.mu.bounds[2], con$norm.sigma.bounds[2]))
  },
  error = no.fit
  )
  
  if (min(quants) < 0) {
    message("Only fit the normal distribution because of negative quantiles.")
    fit.lnorm <- fit.gamma <- fit.weibull <- fit.beta <- no.fit()
  } else {
    fit.lnorm <- tryCatch({
      stats::optim(par = c(con$lnorm.mu.start, con$lnorm.sigma.start),
                   S.theta.lnorm, method = "L-BFGS-B",
                   lower = c(con$lnorm.mu.bounds[1], con$lnorm.sigma.bounds[1]),
                   upper = c(con$lnorm.mu.bounds[2], con$lnorm.sigma.bounds[2]))
    },
    error = no.fit
    )
    fit.gamma <- tryCatch({
      stats::optim(par = c(con$gamma.shape.start, con$gamma.rate.start),
                   S.theta.gamma, method = "L-BFGS-B",
                   lower = c(con$gamma.shape.bounds[1], con$gamma.rate.bounds[1]),
                   upper = c(con$gamma.shape.bounds[2], con$gamma.rate.bounds[2]))
    },
    error = no.fit
    )
    fit.weibull <- tryCatch({
      stats::optim(par = c(con$weibull.shape.start, con$weibull.scale.start),
                   S.theta.weibull, method = "L-BFGS-B",
                   lower = c(con$weibull.shape.bounds[1],
                             con$weibull.scale.bounds[1]),
                   upper = c(con$weibull.shape.bounds[2],
                             con$weibull.scale.bounds[2]))
    },
    error = no.fit
    )
    if (two.sample.default){
      fit.beta <- no.fit(1)
    } else {
      fit.beta <- tryCatch({
        stats::optim(par = c(con$beta.shape1.start, con$beta.shape2.start),
                     S.theta.beta, method = "L-BFGS-B",
                     lower = c(con$beta.shape1.bounds[1],
                               con$beta.shape2.bounds[1]),
                     upper = c(con$beta.shape1.bounds[2],
                               con$beta.shape2.bounds[2]))
      },
      error = no.fit,
      warning = no.fit
      )
    }
  }
  
  values <- c(fit.norm$value, fit.lnorm$value, fit.gamma$value,
              fit.weibull$value, fit.beta$value)
  names(values) <- c("normal", "log-normal", "gamma", "weibull", "beta")
  norm.par <- fit.norm$par
  lnorm.par <- fit.lnorm$par
  gamma.par <- fit.gamma$par
  weibull.par <- fit.weibull$par
  beta.par <- fit.beta$par
  if (length(norm.par) != 1) {
    names(norm.par) <- c("mu", "sigma")
  }
  if (length(lnorm.par) != 1) {
    names(lnorm.par) <- c("mu", "sigma")
  }
  if (length(gamma.par) != 1) {
    names(gamma.par) <- c("shape", "rate")
  }
  if (length(weibull.par) != 1) {
    names(weibull.par) <- c("shape", "scale")
  }
  if (length(beta.par) != 1) {
    names(beta.par) <- c("shape1", "shape2")
  }
  
  num.input <- get.num.input(min.val, q1.val, med.val, q3.val, max.val, n)
  output <- list(norm.par = norm.par, lnorm.par = lnorm.par,
                 gamma.par = gamma.par, weibull.par = weibull.par,
                 beta.par = beta.par, values = values,
                 num.input = num.input, scenario = scenario)
  class(output) <- "qe.fit"
  return(output)
}
check_errors <- function(min.val, q1.val, med.val, q3.val, max.val, n, scenario){
  if (missing(n)) {
    stop("Need to specify n")
  }
  if (is.na(n) | n < 3 | n > 1e6){
    stop("Value of n must be between 3 and 1,000,000")
  }
  if (scenario == 'S1'){
    quants <- c(min.val, med.val, max.val)
    if (!is.numeric(quants)){
      stop(paste('The quantiles are not of type numeric. The values of min.val, med.val, and max.val are:', min.val, med.val, max.val))
    }
    if (is.unsorted(quants)){
      stop(paste('The quantiles are not in increasing order. The values of min.val, med.val, and max.val are:', min.val, med.val, max.val))
    }
    if (length(unique(c(quants))) < 3){
      warning(paste('Some of the quantiles are equal to each other. This can result in unexpected behaviour of the method. The values of min.val, med.val, and max.val are:', min.val, med.val, max.val))
    }
  } else if (scenario == 'S2'){
    quants <- c(q1.val, med.val, q3.val)
    if (!is.numeric(quants)){
      stop(paste('The quantiles are not of type numeric. The values of q1.val, med.val, and q3.val are:', q1.val, med.val, q3.val))
    }
    if (is.unsorted(quants)){
      stop(paste('The quantiles are not in increasing order. The values of q1.val, med.val, and q3.val are:', q1.val, med.val, q3.val))
    }
    if (length(unique(quants)) < 3){
      warning(paste('Some of the quantiles are equal to each other. This can result in unexpected behaviour of the method. The values of q1.val, med.val, and q3.val are:', q1.val, med.val, q3.val))
    }
  } else if (scenario == 'S3'){
    quants <- c(min.val, q1.val, med.val, q3.val, max.val)
    if (!is.numeric(quants)){
      stop(paste('The quantiles are not of type numeric. The values of min.val, q1.val, med.val, q3.val, and max.val are:', min.val, q1.val, med.val, q3.val, max.val))
    }
    if (is.unsorted(quants)){
      stop(paste('The quantiles are not in increasing order. The values of min.val, q1.val, med.val, q3.val, and max.val are:', min.val, q1.val, med.val, q3.val, max.val))
    }
    if (length(unique(quants)) < 5){
      warning(paste('Some of the quantiles are equal to each other. This can result in unexpected behaviour of the method. The values of min.val, q1.val, med.val, q3.val, and max.val are:', min.val, q1.val, med.val, q3.val, max.val))
    }
  }
}
qe.mean.sd2 <- function(min.val, q1.val, med.val, q3.val, max.val, n,p,
                        qe.fit.control = list()) {
  args <- as.list(environment())
  x <- qe.fit(min.val = min.val, q1.val = q1.val, med.val = med.val,
              q3.val = q3.val, max.val = max.val, n = n,p,
              qe.fit.control = qe.fit.control)
  selected.dist <- names(which.min(x$values))
  ests <- get.mean.sd(x, selected.dist)
  output <- list(est.mean = ests$est.mean, est.sd = ests$est.sd,
                 selected.dist = selected.dist, values = x$values,
                 args = args, scenario = x$scenario, fitted.dists = x)
  for (dist.name in names(x$values)){
    ests.all <- get.mean.sd(x, dist.name)
    output[paste0(dist.name, '.est.mean')] <- ests.all$est.mean
    output[paste0(dist.name, '.est.sd')] <- ests.all$est.sd
  }
  class(output) <- "qe.mean.sd"
  return(output)
}
summary.qe.mean.sd <- function(object, digits = 5, ...) {
  res.mat <- matrix(nrow = length(object$values), ncol = 3)
  rownames(res.mat) <- names(object$values[order(object$values)])
  colnames(res.mat) <- c("Mean", "SD", "SS")
  for (dist.name in row.names(res.mat)){
    res.mat[dist.name, "Mean"] <- object[[paste0(dist.name, ".est.mean")]]
    res.mat[dist.name, "SD"] <- object[[paste0(dist.name, ".est.sd")]]
    res.mat[dist.name, "SS"] <- unname(object$values[dist.name])
  }
  res.mat <- round(res.mat, digits = digits)
  return(res.mat)
}


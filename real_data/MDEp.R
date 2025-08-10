
library(meta)
library(estmeansd)
library(pracma)


normal<-function(data,n,startm,startsd,lowsd=1e-3,upsd=50){
  
  N<-n
  start_params <- c(mu = startm, sigma = startsd)
  
  
  objective_function_w <- function(params, q, p, N) {
    mu <- params[1]
    sigma <- params[2]
    
    # Ensure sigma is positive
    if (sigma <= 0) return(Inf)
    
    # Calculate predicted probabilities
    p_pred <- pnorm(q, mean = mu, sd = sigma)
    
    x<-qnorm(p,mu,sigma)
    # Calculate density at the quantile
    FF <- pnorm(x, mean = mu, sd = sigma)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-FF[i]*(1-FF[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- p - p_pred
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
    lower = c(data$q[2],lowsd),
    upper = c(data$q[4],upsd),
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
} #normal
weibull<-function(data,n,starta,startb,lowa=1e-3,upa=100,lowb=1e-3,upb=100){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    p_pred <- pweibull(q, shape = bsb, scale = bsa)
    
    x<-qweibull(p,shape = bsb, scale = bsa)
    # Calculate density at the quantile
    FF <- pweibull(x, shape = bsb, scale = bsa)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-FF[i]*(1-FF[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- p - p_pred
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
lognormal<-function(data,n,startm,startsd,lowsd=1e-3,upsd=10){
  
  N<-n
  start_params <- c(mu = startm, sigma = startsd)
  
  
  objective_function_w <- function(params, q, p, N) {
    mu <- params[1]
    sigma <- params[2]
    
    
    # Calculate predicted probabilities
    # Calculate predicted probabilities
    # Calculate predicted probabilities
    p_pred <- plnorm(q, meanlog = mu, sdlog = sigma)
    
    x<-qlnorm(p,meanlog=mu,sdlog=sigma)
    # Calculate density at the quantile
    FF <- plnorm(x, meanlog = mu, sdlog = sigma)
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-FF[i]*(1-FF[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- p - p_pred
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
beta<-function(data,n,starta,startb){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    p_pred <- pbeta(q, shape1 = bsa, shape2 = bsb)
    
    x<-qbeta(p,shape1 = bsa,shape2 = bsb)
    # Calculate density at the quantile
    FF <- pbeta(x, shape1 = bsa, shape2 = bsb)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-FF[i]*(1-FF[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- p - p_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
    weighted_residuals 
  }
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "Nelder-Mead",  
    control = list(trace = 0),# Trace for debugging
    hessian = FALSE
  )
  result_w
}#beta
egamma<-function(data,n,starta,startb,lowa=0,upa=50,lowb=0,upb=1/50){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    
    # Calculate predicted probabilities
    p_pred <- pgamma(q, shape = bsa, scale = bsb)
    
    x<-qgamma(p,shape = bsa, scale = bsb)
    # Calculate density at the quantile
    FF <- pgamma(x, shape = bsa, scale = bsb)
    
    # Calculate variance and weights
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-FF[i]*(1-FF[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- qr.solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- p - p_pred
    weighted_residuals <- t(residuals) %*% weights %*% residuals
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
    hessian = TRUE
  )
  result_w
}#gamma

ds<-read.csv("realdata.csv") #Data for PHQ-9 scores

colnames(ds)<-c("Study","a","q1","m","q3","b","n")


dstrue<-read.csv("truevalues.csv")#Include true mean and 95% CI
dstrue$SE<-(dstrue$U-dstrue$L)/(2*1.96)

model<-c()
modelqe<-c()

for (i in 1:nrow(ds)){

  value<-rep(9999,5)
  dat<-ds[i,2:7]
  a<-as.numeric(dat[1])/27
  q1<-as.numeric(dat[2])/27
  m<-as.numeric(dat[3])/27
  q3<-as.numeric(dat[4])/27
  b<-as.numeric(dat[5])/27
  n<-as.numeric(dat[6])
  data <- data.frame(
    p = c(0.625/n, 0.25,0.5, 0.75,1-0.625/n),
    q = (c(a,q1,m,q3,b))
  )
  
 
  ##model
  
  mm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
  theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
  theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
  var<-((b-a)/theta1+(q3-q1)/theta2)^2
  startmn<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
  startsdn<-(b-a)/theta1+(q3-q1)/theta2
  startmln<-log(mm/sqrt(1+var/mm^2))
  startsdln<-sqrt(log(1+var/mm^2))
  startbw<-(sqrt(var)/mm)^(-1.086)
  startaw<-mm/gamma(1+1/startbw)
  startab<-mm*(mm*(1-mm)/var-1)
  startbb<-(1-mm)*(mm*(1-mm)/var-1)
  startbg<-var/mm
  startag<-mm/startbg
  
  
  er1<-tryCatch(normal(data,n,startmn,startsdn),error=function(e) 1)
  if (typeof(er1)=="list"){
    fit1<-normal(data,n,startmn,startsdn)
    if(fit1$convergence==0){
      pre<-pnorm(data$q, mean = as.numeric(fit1$par[1]), sd = as.numeric(fit1$par[2]))
      residual<-(data$p-pre)^2
      value[1]<-sum(residual)}}
  
  
  er2<-tryCatch(weibull(data,n,startaw,startbw),error=function(e) 1)
  if (typeof(er2)=="list"){
    fit2<-weibull(data,n,startaw,startbw)
    if (fit2$convergence==0){
      pre<-pweibull(data$q, scale = as.numeric(fit2$par[1]), shape = as.numeric(fit2$par[2]))
      residual<-(data$p-pre)^2
      value[2]<-sum(residual)}}
  
  er3<-tryCatch(lognormal(data,n,startmln,startsdln),error=function(e) 1)
  if (typeof(er3)=="list"){
    fit3<-lognormal(data,n,startmln,startsdln)
    pre<- plnorm(data$q, meanlog = as.numeric(fit3$par[1]), sdlog = as.numeric(fit3$par[2]))
    residual<-(data$p-pre)^2
    value[3]<-sum(residual)}
  
  er4<-tryCatch(beta(data,n,startab,startbb),error=function(e) 1)
  if (typeof(er4)=="list"){
    fit4<-beta(data,n,startab,startbb)
    pre<- pbeta(data$q, shape1 = as.numeric(fit4$par[1]), shape2 = as.numeric(fit4$par[2]))
    residual<-(data$p-pre)^2
    value[4]<-sum(residual)
    }
  
  er5<-tryCatch(egamma(data,n,startag,startbg),error=function(e) 1)
  if(typeof(er5)=="list"){
    fit5<-egamma(data,n,startag,startbg)
    if (fit5$convergence==0){
      pre<- pgamma(data$q, shape = as.numeric(fit5$par[1]), scale = as.numeric(fit5$par[2]))
      residual<-(data$p-pre)^2
      value[5]<-sum(residual)
    }}
  
  if (which(value==min(value))==1){fit<-fit1
  model[i]="normal"}
  if (which(value==min(value))==2){fit<-fit2
  model[i]="weibull"}
  if (which(value==min(value))==3){fit<-fit3
  model[i]="lognormal"}
  if (which(value==min(value))==4){fit<-fit4
  model[i]="beta"}
  if (which(value==min(value))==5){fit<-fit5
  model[i]="gamma"}
  
  
  
  if (model[i]=="normal"){
    fit<-normal(data,n,startmn,startsdn)
    ds$mw[i]<-fit$par[1]
    ds$sdw[i]<-fit$par[2]
  }
  
  
  
  if(model[i]=="weibull"){
    fit<-weibull(data,n,startaw,startbw)
    ds$mw[i]<-fit$par[1]*gamma(1+1/fit$par[2])
    ds$sdw[i]<-as.numeric(sqrt(fit$par[1]^2*(gamma(1+2/fit$par[2])-(gamma(1+1/fit$par[2]))^2)))
  }
  
  if(model[i]=="gamma"){
    fit<-egamma(data,n,startag,startbg)
    ds$mw[i]<-as.numeric(fit$par[1]*fit$par[2])
    ds$sdw[i]<-sqrt(as.numeric(fit$par[1]*fit$par[2]^2))
  }
  
  if(model[i]=="beta"){
    fit<-beta(data,n,startab,startbb)
    ds$mw[i]<-as.numeric(fit$par[1]/(fit$par[2]+fit$par[1]))
    ds$sdw[i]<-as.numeric(sqrt((fit$par[1]*fit$par[2])/((fit$par[1]+fit$par[2])^2*(fit$par[1]+fit$par[2]+1))))
  }
  
  if(model[i]=="lognormal"){
    fit<-lognormal(data,n,startmln,startsdln)
    ds$mw[i]<-exp(fit$par[1]+fit$par[2]^2/2)
    ds$sdw[i]<-sqrt((exp(fit$par[2]^2)-1)*(exp(2*fit$par[1]+fit$par[2]^2)))
  }
  ds$mw[i]<-ds$mw[i]*27
  ds$sdw[i]<-ds$sdw[i]*27
               
  a1<-as.numeric(dat[1])/27
  q11<-as.numeric(dat[2])/27
  m1<-as.numeric(dat[3])/27
  q31<-as.numeric(dat[4])/27
  b1<-as.numeric(dat[5])/27
         
               
  ds$luo[i]<- (((2.2/(2.2+n^0.75))*((a1+b1)/2))+((0.7-0.72/(n^0.55))*((q11+q31)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m1)*27
  
  theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
  theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
  ds$shi[i] <-((b1-a1)/theta1+(q31-q11)/theta2)*27   
  
  ds$mbc[i]<-(as.numeric(bc.mean.sd(min.val = a1,q1.val = q11,med.val = m1,q3.val = q31,max.val = b1,n=n)[1]))*27
  ds$sdbc[i]<-(as.numeric(bc.mean.sd(min.val = a1,q1.val = q11,med.val = m1,q3.val = q31,max.val = b1,n=n)[2]))*27
  
  ds$mqe[i]<-(as.numeric(qe.mean.sd(min.val = a1,q1.val = q11,med.val = m1,q3.val = q31,max.val = b1,n=n)[1]))*27
  ds$sdqe[i]<-(as.numeric(qe.mean.sd(min.val = a1,q1.val = q11,med.val = m1,q3.val = q31,max.val = b1,n=n)[2]))*27
  ss<-summary(qe.mean.sd(min.val = a1,q1.val = q11,med.val = m1,q3.val = q31,max.val = b1,n=n))
  modelqe[i]<-rownames(ss)[1]
  }


#meta analysis
#our model
mymeta1<-metagen(TE=mw,seTE=sdw/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta1)

mymeta2<-metagen(TE=luo,seTE=shi/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta2)

mymeta3<-metagen(TE=mqe,seTE=sdqe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta3)

mymeta4<-metagen(TE=Mean,seTE=SE,studlab=Study,data=dstrue,method.tau="REML")
summary(mymeta5)


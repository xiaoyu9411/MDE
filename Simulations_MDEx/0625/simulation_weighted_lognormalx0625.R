library(stats) 
library(devtools)
library(pracma)
#weighted model function

weightednls<-function(data,n,startm,startsd){
  
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
    weights <- qr.solve(cov)
    
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
    lower = c(log(data$q[2]), 1e-3),
    upper = c(log(data$q[4]), 10),
    hessian = FALSE
  )
  result_w
}


set.seed(1456)
#lognormal distribution
meant<-1.5
sdt<-0.5
truemean<-exp(meant+sdt^2/2)
truesd<-sqrt((exp(sdt^2)-1)*(exp(2*meant+sdt^2)))
nnum<-seq(from=30,to=1000,by=10)
nsim<-10000

it<-0
rbmeanw<-c()
rbsdw<-c()


msemeanw<-c()
msesdw<-c()
msesamplemean<-c()
msesamplesd<-c()


it<-0
xbar<-c()
sd<-c()
pw<-c()
summarydat<-matrix(nrow=length(nnum),ncol=8)


for (i in nnum ){
  mw<-c()
  sdw<-c()
  xbar<-c()
  for (j in 1:nsim){
    n<-i
    datset<-rlnorm(n,meanlog=meant,sdlog=sdt)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])
    q025<-as.numeric(quantile(datset,p=0.025))
    q1<-as.numeric(s[2])
    m<-as.numeric(s[3])
    q3<-as.numeric(s[5])
    q975<-as.numeric(quantile(datset,p=0.975))
    b<-as.numeric(s[6])
    
  mm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    var<-((b-a)/theta1+(q3-q1)/theta2)^2
    startm<-log(mm/sqrt(1+var/mm^2))
    startsd<-sqrt(log(1+var/mm^2))
    #startm<-meant
    #startm<-sdt
    
    data <- data.frame(
      p = c(0.625/n, 0.25, 0.5, 0.75,1-0.625/n),
      q = c(a,q1,m,q3,b )
    )
    #weighted model
    er<-tryCatch(weightednls(data,n,startm,startsd),error=function(e) 1)
    if (typeof(er)=="list"){
    fit <- weightednls(data,n,startm,startsd)
    if (fit$convergence==0){
    mw[j]<-exp(fit$par[1]+fit$par[2]^2/2)
    sdw[j]<-sqrt((exp(fit$par[2]^2)-1)*(exp(2*fit$par[1]+fit$par[2]^2)))
    }
    }
    
  }
  
  it<-it+1
  
  pw[it]<-mean(is.na(mw)|mw>100)
  
  rbmeanw[it]<-mean((mw[(!is.na(mw))&mw<100]-truemean)/truemean)
  rbsdw[it]<-mean((sdw[(!is.na(sdw))&sdw<100]-truesd)/truesd)
  msemeanw[it]<-sqrt(mean((mw[(!is.na(mw))&mw<100]-truemean)^2))
  msesdw[it]<-sqrt(mean((sdw[(!is.na(sdw))&sdw<100]-truesd)^2))
  msesamplemean[it]<-sqrt(mean((xbar-truemean)^2))
  msesamplesd[it]<-sqrt(mean((sd-truesd)^2))
  
  
  summarydat[it,]<-c(i,rbmeanw[it],rbsdw[it],msemeanw[it],msesdw[it],msesamplemean[it],msesamplesd[it],pw[it])
    write.table(summarydat[it,],file="resultlognormalx0625.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}

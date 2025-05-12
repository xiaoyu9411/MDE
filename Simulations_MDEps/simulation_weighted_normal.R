library(stats) 
library(devtools)
devtools::install_github("stmcg/estmeansd",lib="/usr3/graduate/rainie")
library(metaBLUE,lib.loc ="/usr3/graduate/rainie" )
library(estmeansd,lib.loc ="/usr3/graduate/rainie" )
library(pracma)
#weighted model function

weightednls<-function(data,n,startm,startsd){
  
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
    f <- dnorm(x, mean = mu, sd = sigma)
    
    cov<-matrix(NA,nrow=length(p),ncol=length(p))
    # Calculate variance and weights
    for (i in 1:length(p)){
      for (j in i:length(p)){
        cov[i,j]<-p[i]*(1-p[j])/(N * f[i]*f[j])
        cov[j,i]<-cov[i,j]
      }
    }
    
    
    #variance <- p * (1 - p) / (N * f^2)
    weights <- solve(cov)
    
    # Calculate weighted sum of squared residuals
    residuals <- p - p_pred
    weighted_residuals <- t(residuals) %*% sqrtm(weights)$B %*% residuals
    weighted_residuals

  }
  
  
  result_w <- optim(
    par = start_params,
    fn = objective_function_w,
    q = data$q,
    p = data$p,
    N = N,
    method = "L-BFGS-B",
    lower = c(data$q[2],1e-3),
    upper = c(data$q[4],50),
    control = list(trace = 0),# Trace for debugging
    hessian = TRUE
  )
  result_w
}


set.seed(135)


truemean<-5
truesd<-1
nnum<-seq(from=30,to=1000,by=10)
nsim<-10000
it<-0
rbmeanw<-c()
rbsdw<-c()
rbmeanluo<-c()
rbsdshi<-c()

rbmeanqe<-c()
rbsdqe<-c()

msemeanluo<-c()
msesdshi<-c()
msemeanw<-c()
msesdw<-c()
msesamplemean<-c()
msesamplesd<-c()

msemeanqe<-c()
msesdqe<-c()

it<-0
xbar<-c()
sd<-c()

pw<-c()
pqe<-c()
summarydat<-matrix(nrow=length(nnum),ncol=17)



for (i in nnum ){
  mw<-c()
  sdw<-c()
  xbar<-c()
  sd<-c()
  mluo<-c()
  sdshi<-c()

  mqe<-c()
  sdqe<-c()
  for (j in 1:nsim){
    n<-i
    datset<-rnorm(n,mean=truemean,sd=truesd)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])
    q1<-as.numeric(s[2])
    m<-as.numeric(s[3])
    q3<-as.numeric(s[5])
    b<-as.numeric(s[6])
    
    ##model
    mm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    var<-((b-a)/theta1+(q3-q1)/theta2)^2
    startm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    startsd<-(b-a)/theta1+(q3-q1)/theta2
    
    data <- data.frame(
      p = c(0.625/n,0.25, 0.5, 0.75, 1-0.625/n),
      q = c(a,q1,m,q3,b )
    )
    #weighted model
    fit <- weightednls(data,n,startm,startsd)
    if (fit$convergence==0){

  
    mw[j]<-fit$par[1]
    sdw[j]<-fit$par[2]
    }
    #luo's method
    mluo[j]<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    #Shi's method
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    sdshi[j]<-(b-a)/theta1+(q3-q1)/theta2
    

    ss<-summary(qe.mean.sd(min.val = a,q1.val = q1,med.val = m,q3.val = q3,max.val = b,n=n))
    smodel<-ss[rownames(ss)=="normal",]
    mqe[j]<-as.numeric(smodel[1])
    sdqe[j]<-as.numeric(smodel[2])
    
    
    
  }
  it<-it+1

  pw[it]<-mean(is.na(mw)|mw>100)
  pqe[it]<-mean(is.na(mqe))
  rbmeanw[it]<-mean((mw[(!is.na(mw))&mw<100]-truemean)/truemean)
  rbsdw[it]<-mean((sdw[(!is.na(sdw))&sdw<100]-truesd)/truesd)
  rbmeanluo[it]<-mean((mluo-truemean)/truemean)
  rbsdshi[it]<- mean((sdshi-truesd)/truesd)

  rbmeanqe[it]<-mean((mqe-truemean)/truemean)
  rbsdqe[it]<- mean((sdqe-truesd)/truesd)
  
  
  msemeanluo[it]<-sqrt(mean((mluo-truemean)^2))
  msesdshi[it]<-sqrt(mean((sdshi-truesd)^2))
  msemeanw[it]<-sqrt(mean((mw[(!is.na(mw))&mw<100]-truemean)^2))
  msesdw[it]<-sqrt(mean((sdw[(!is.na(sdw))&sdw<100]-truesd)^2))
  msesamplemean[it]<-sqrt(mean((xbar-truemean)^2))
  msesamplesd[it]<-sqrt(mean((sd-truesd)^2))
  

  msemeanqe[it]<-sqrt(mean((mqe-truemean)^2))
  msesdqe[it]<-sqrt(mean((sdqe-truesd)^2))
    summarydat[it,]<-c(i,rbmeanw[it],rbsdw[it],rbmeanluo[it],rbsdshi[it],rbmeanqe[it],rbsdqe[it],msemeanw[it],msesdw[it],msemeanluo[it],msesdshi[it],msemeanqe[it],msesdqe[it],msesamplemean[it],msesamplesd[it],pw[it],pqe[it])
  write.table(summarydat[it,],file="resultnormalw_T.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}

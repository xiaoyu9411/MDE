library(stats) 
library(devtools)
devtools::install_github("stmcg/estmeansd",lib="/usr3/graduate/rainie")
library(metaBLUE,lib.loc ="/usr3/graduate/rainie" )
library(estmeansd,lib.loc ="/usr3/graduate/rainie" )
library(pracma)
#weighted model function

weightednls<-function(data,n,starta,startb){
  
  N<-n
  start_params <- c(bsa = starta, bsb = startb)
  
  
  objective_function_w <- function(params, q, p, N) {
    bsa <- params[1]
    bsb <- params[2]
    
    
    # Calculate predicted probabilities
    p_pred <- pbeta(q, shape1 = bsa, shape2 = bsb)
    
    x<-qbeta(p,shape1 = bsa,shape2 = bsb)
    # Calculate density at the quantile
    f <- dbeta(x, shape1 = bsa, shape2 = bsb)
    
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
    method = "Nelder-Mead",  
    control = list(trace = 0),# Trace for debugging
    hessian = TRUE
  )
  result_w
}


set.seed(135)

#beta distribution
sa<-4
sb<-2
truemean<-8*sa/(sb+sa)
truesd<-8*sqrt((sa*sb)/((sa+sb)^2*(sa+sb+1)))
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
    datset<-8*rbeta(n,shape1=sa,shape2=sb)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])/8
    q1<-as.numeric(s[2])/8
    m<-as.numeric(s[3])/8
    q3<-as.numeric(s[5])/8
    b<-as.numeric(s[6])/8
    
    ##model
   

    mm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    var<-((b-a)/theta1+(q3-q1)/theta2)^2
    starta<-mm*(mm*(1-mm)/var-1)
    startb<-(1-mm)*(mm*(1-mm)/var-1)
    
    data <- data.frame(
      p = c(0.625/n, 0.25, 0.5, 0.75, 1-0.625/n),
      q = c(a,q1,m,q3,b )
    )
    #weighted model
    er<-tryCatch(weightednls(data,n,starta,startb),error=function(e) 1)
    if (typeof(er)=="list"){
    fit <- weightednls(data,n,starta,startb)
    if (fit$convergence==0){
  
    mw[j]<-8*as.numeric(fit$par[1]/(fit$par[2]+fit$par[1]))
    sdw[j]<-8*as.numeric(sqrt((fit$par[1]*fit$par[2])/((fit$par[1]+fit$par[2])^2*(fit$par[1]+fit$par[2]+1))))
    #luo's method
    }
    }
    mluo[j]<-8*(((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m)
    #Shi's method
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    sdshi[j]<-8*((b-a)/theta1+(q3-q1)/theta2)
    
    
    
    ss<-summary(qe.mean.sd(min.val = a,q1.val = q1,med.val = m,q3.val = q3,max.val = b,n=n))
    smodel<-ss[rownames(ss)=="beta",]
    mqe[j]<-8*as.numeric(smodel[1])
    sdqe[j]<-8*as.numeric(smodel[2])
    
  }
    it<-it+1

    pw[it]<-mean(is.na(mw)|mw>100)
    pqe[it]<-mean(is.na(mqe))
    rbmeanw[it]<-mean((mw[(!is.na(mw))&mw<100]-truemean)/truemean)
    rbsdw[it]<-mean((sdw[(!is.na(sdw))&sdw<100]-truesd)/truesd)
    rbmeanluo[it]<-mean((mluo-truemean)/truemean)
    rbsdshi[it]<- mean((sdshi-truesd)/truesd)

    rbmeanqe[it]<-mean((mqe[(!is.na(mqe))]-truemean)/truemean)
    rbsdqe[it]<- mean((sdqe[(!is.na(sdqe))]-truesd)/truesd)
    
    
    msemeanluo[it]<-sqrt(mean((mluo-truemean)^2))
    msesdshi[it]<-sqrt(mean((sdshi-truesd)^2))
    msemeanw[it]<-sqrt(mean((mw[(!is.na(mw))&mw<100]-truemean)^2))
    msesdw[it]<-sqrt(mean((sdw[(!is.na(sdw))&sdw<100]-truesd)^2))
    msesamplemean[it]<-sqrt(mean((xbar-truemean)^2))
    msesamplesd[it]<-sqrt(mean((sd-truesd)^2))
    

    msemeanqe[it]<-sqrt(mean((mqe[(!is.na(mqe))]-truemean)^2))
    msesdqe[it]<-sqrt(mean((sdqe[(!is.na(sdqe))]-truesd)^2))
    
    summarydat[it,]<-c(i,rbmeanw[it],rbsdw[it],rbmeanluo[it],rbsdshi[it],rbmeanqe[it],rbsdqe[it],msemeanw[it],msesdw[it],msemeanluo[it],msesdshi[it],msemeanqe[it],msesdqe[it],msesamplemean[it],msesamplesd[it],pw[it],pqe[it])
  write.table(summarydat[it,],file="resultbetacom_T.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}

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
    datset<-rlnorm(n,meanlog=meant,sdlog=sdt)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])
    q1<-as.numeric(s[2])
    m<-as.numeric(s[3])
    q3<-as.numeric(s[5])
    b<-as.numeric(s[6])
    
    mm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    var<-((b-a)/theta1+(q3-q1)/theta2)^2
    startm<-log(mm/sqrt(1+var/mm^2))
    startsd<-sqrt(log(1+var/mm^2))
    
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
    w1<-2.2/(2.2+n^0.75)
    w2<-0.7-0.72*n^(-0.55)
    mu3<-w1*(log(a)+log(b))/2+w2*(log(q1)+log(q3))/2+(1-w1-w2)*log(m)
    w3<-1/(1+0.07*n^0.6)
    epsi<-2*qnorm((n-0.375)/(n+0.25))
    ita<-2*qnorm((0.75*n-0.125)/(n+0.25))
    sig<-(w3*((log(b)-log(a))/epsi)+(1-w3)*((log(q3)-log(q1))/ita))^2*(1+0.28/(log(n))^2)^-1
    mluo[j]<-exp(mu3+sig/2)*(1+(0.405/n)*sig+(0.315/n)*sig^2)^-1
    #Shi's method
    sdshi[j]<-sqrt((exp(2*mu3+2*sig)*(1+(1.62/n*sig)+(5.04/n)*sig^2)^-1)-(exp(2*mu3+sig)*(1+(1.62/n*sig)+(1.26/n)*sig^2)^-1))
 

    ss<-summary(qe.mean.sd(min.val = a,q1.val = q1,med.val = m,q3.val = q3,max.val = b,n=n))
    smodel<-ss[rownames(ss)=="log-normal",]
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
    write.table(summarydat[it,],file="resultlognormalp.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}

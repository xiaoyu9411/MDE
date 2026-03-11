source("functions.R")
set.seed(1456)


#lognormal distribution
meant<-1.5
sdt<-0.5
truemean<-exp(meant+sdt^2/2)
truesd<-sqrt((exp(sdt^2)-1)*(exp(2*meant+sdt^2)))
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
    datset<-rlnorm(n,meanlog=meant,sdlog=sdt)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])

    m<-as.numeric(s[3])

    b<-as.numeric(s[6])
    
    #Shi's method
    w1<-4/(4+n^0.75)
    mu3<-w1*(log(a)+log(b))/2+(1-w1)*log(m)
    epsi<-2*qnorm((n-0.375)/(n+0.25))
    sig<-(((log(b)-log(a))/epsi)^2)*(1.01+0.25/((log(n))^2))^-1
    mluo[j]<-exp(mu3+sig/2)*(1+(0.565/n)*sig+(0.37/n)*sig^2)^-1
    sdshi[j]<-sqrt((exp(2*mu3+2*sig)*(1+(2.26/n*sig)+(5.92/n)*sig^2)^-1)-(exp(2*mu3+sig)*(1+(2.26/n*sig)+(1.48/n)*sig^2)^-1))
    
    
    #qe
    ss<-qe.mean.sd(min.val = a,med.val = m,max.val = b,n=n)
    mqe[j]<-ss$est.mean
    sdqe[j]<-ss$est.sd
    disqeall[j]<-ss$selected.dist
    #wqe
    wqe<-wqe.mean.sd(min.val = a,med.val = m,max.val = b,n=n)
    mwqe[j]<-wqe$est.mean
    sdwqe[j]<-wqe$est.sd
    diswqeall[j]<-wqe$selected.dist
    
    #mdex
    MDEx<-MDEx.mean.sd(min.val = a,med.val = m,max.val = b,n=n)
    mmdex[j]<-MDEx$est.mean
    sdmdex[j]<-MDEx$est.sd
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
  
  disqe[it]<-mean(disqeall=="log-normal")
  diswqe[it]<-mean(diswqeall=="log-normal")
  dismdex[it]<-mean(dismdexall=="log-normal")
  
  
  summarydat[it,]<-c(i,rbmluo[it],rbsdluo[it],rbmqe[it],rbsdqe[it],
                     rbmwqe[it],rbsdwqe[it],rbmmdex[it],rbsdmdex[it],
                     msemluo[it],msesdluo[it],msemqe[it],msesdqe[it],
                     msemwqe[it],msesdwqe[it],msemmdex[it],msesdmdex[it],
                     msemsample[it],msesdsample[it],
                     pqe[it],pwqe[it],pmdex[it],
                     disqe[it],diswqe[it],dismdex[it])
  write.table(summarydat[it,],file="resultlognormalS1.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}

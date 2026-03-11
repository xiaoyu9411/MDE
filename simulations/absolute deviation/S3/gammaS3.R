source("functions.R")
set.seed(1456)


#gamma distribution
sa<-2.5 #shape
sb<-0.5 #rate
truemean<-sa/sb
truesd<-sqrt(sa*(1/sb)^2)
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
    datset<-rgamma(n,shape=sa,rate=sb)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])
    q1<-as.numeric(s[2])
    m<-as.numeric(s[3])
    q3<-as.numeric(s[5])
    b<-as.numeric(s[6])
    
    #luo's method
    mluo[j]<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    #Shi's method
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    sdshi[j]<-(b-a)/theta1+(q3-q1)/theta2
    
    #qe
    ss<-qe.mean.sd(min.val = a,q1.val = q1,med.val = m,q3.val = q3,max.val = b,n=n)
    mqe[j]<-ss$est.mean
    sdqe[j]<-ss$est.sd
    disqeall[j]<-ss$selected.dist
    #wqe
    wqe<-wqe.mean.sd(min.val = a,q1.val = q1,med.val = m,q3.val = q3,max.val = b,n=n)
    mwqe[j]<-wqe$est.mean
    sdwqe[j]<-wqe$est.sd
    diswqeall[j]<-wqe$selected.dist
    
    #mdex
    MDEx<-MDEx.mean.sd(min.val = a,q1.val = q1,med.val = m,q3.val = q3,max.val = b,n=n)
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
  
  disqe[it]<-mean(disqeall=="gamma")
  diswqe[it]<-mean(diswqeall=="gamma")
  dismdex[it]<-mean(dismdexall=="gamma")
  
  
  summarydat[it,]<-c(i,rbmluo[it],rbsdluo[it],rbmqe[it],rbsdqe[it],
                     rbmwqe[it],rbsdwqe[it],rbmmdex[it],rbsdmdex[it],
                     msemluo[it],msesdluo[it],msemqe[it],msesdqe[it],
                     msemwqe[it],msesdwqe[it],msemmdex[it],msesdmdex[it],
                     msemsample[it],msesdsample[it],
                     pqe[it],pwqe[it],pmdex[it],
                     disqe[it],diswqe[it],dismdex[it])
  write.table(summarydat[it,],file="resultgammaS3.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}


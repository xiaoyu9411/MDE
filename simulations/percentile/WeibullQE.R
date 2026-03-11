source("functionsQE.R")

set.seed(1456)

#Weibull distribution
sa<-6 #scale
sb<-2#shape
truemean<-sa*gamma(1+1/sb)
truesd<-sqrt(sa^2*(gamma(1+2/sb)-(gamma(1+1/sb))^2))
nnum<-seq(from=30,to=1000,by=10)
nsim<-10000
it<-0
rbmean1<-c()
rbsd1<-c()
rbmean05<-c()
rbsd05<-c()
rbmean0625<-c()
rbsd0625<-c()
rbmeanluo<-c()
rbsdshi<-c()

msemean1<-c()
msesd1<-c()
msemean05<-c()
msesd05<-c()
msemean0625<-c()
msesd0625<-c()
msemeanluo<-c()
msesdshi<-c()
msesamplemean<-c()
msesamplesd<-c()

p1<-c()
p05<-c()
p0625<-c()


it<-0
xbar<-c()
sd<-c()
summarydat<-matrix(nrow=length(nnum),ncol=22)



for (i in nnum ){
  m1<-c()
  sd1<-c()
  xbar<-c()
  sd<-c()
  m05<-c()
  sd05<-c()
  m0625<-c()
  sd0625<-c()
  mluo<-c()
  sdshi<-c()
  for (j in 1:nsim){
    n<-i
    datset<-rweibull(n,shape=sb,scale=sa)
    xbar[j]<-mean(datset)
    sd[j]<-sd(datset)
    s<-summary(datset)
    a<-as.numeric(s[1])
    q1<-as.numeric(s[2])
    m<-as.numeric(s[3])
    q3<-as.numeric(s[5])
    b<-as.numeric(s[6])
    

    mluo[j]<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
    #Shi's method
    theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
    theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
    sdshi[j]<-(b-a)/theta1+(q3-q1)/theta2
    
    
    ss1<-summary.qe.mean.sd(qe.mean.sd2(a,q1,m,q3,b,n=n,p=1))
    smodel1<-ss1[rownames(ss1)=="weibull",]
    m1[j]<-as.numeric(smodel1[1])
    sd1[j]<-as.numeric(smodel1[2])
    
    ss05<-summary.qe.mean.sd(qe.mean.sd2(a,q1,m,q3,b,n=n,p=0.5))
    smodel05<-ss05[rownames(ss05)=="weibull",]
    m05[j]<-as.numeric(smodel05[1])
    sd05[j]<-as.numeric(smodel05[2])
    
    
    ss0625<-summary.qe.mean.sd(qe.mean.sd2(a,q1,m,q3,b,n=n,p=0.625))
    smodel0625<-ss0625[rownames(ss0625)=="weibull",]
    m0625[j]<-as.numeric(smodel0625[1])
    sd0625[j]<-as.numeric(smodel0625[2])
  }
  it<-it+1
  p1[it]<-mean(is.na(m1))
  p05[it]<-mean(is.na(m05))
  p0625[it]<-mean(is.na(m0625))
  rbmeanluo[it]<-mean((mluo-truemean)/truemean)
  rbsdshi[it]<- mean((sdshi-truesd)/truesd)
  rbmean1[it]<-mean((m1[(!is.na(m1))]-truemean)/truemean)
  rbsd1[it]<- mean((sd1[(!is.na(sd1))]-truesd)/truesd)
  rbmean05[it]<-mean((m05[(!is.na(m05))]-truemean)/truemean)
  rbsd05[it]<- mean((sd05[(!is.na(sd05))]-truesd)/truesd)
  rbmean0625[it]<-mean((m0625[(!is.na(m0625))]-truemean)/truemean)
  rbsd0625[it]<- mean((sd0625[(!is.na(sd0625))]-truesd)/truesd)
  

  msemean1[it]<-sqrt(mean((m1[(!is.na(m1))]-truemean)^2))
  msesd1[it]<-sqrt(mean((sd1[(!is.na(sd1))]-truesd)^2))
  msemean05[it]<-sqrt(mean((m05[(!is.na(m05))]-truemean)^2))
  msesd05[it]<-sqrt(mean((sd05[(!is.na(sd05))]-truesd)^2))
  msemean0625[it]<-sqrt(mean((m0625[(!is.na(m0625))]-truemean)^2))
  msesd0625[it]<-sqrt(mean((sd0625[(!is.na(sd0625))]-truesd)^2))
  msemeanluo[it]<-sqrt(mean((mluo-truemean)^2))
  msesdshi[it]<-sqrt(mean((sdshi-truesd)^2))
  msesamplemean[it]<-sqrt(mean((xbar-truemean)^2))
  msesamplesd[it]<-sqrt(mean((sd-truesd)^2))
  
  summarydat[it,]<-c(i,rbmeanluo[it],rbsdshi[it],rbmean1[it],rbsd1[it],rbmean05[it],rbsd05[it],rbmean0625[it],rbsd0625[it],msemeanluo[it],msesdshi[it],msemean1[it],msesd1[it],msemean05[it],msesd05[it],msemean0625[it],msesd0625[it],msesamplemean[it],msesamplesd[it],p1[it],p05[it],p0625[it])
  write.table(summarydat[it,],file="resultweibullQE.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}





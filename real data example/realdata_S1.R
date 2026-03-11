
source("functions.R")

ds<-read.csv("realdata.csv") #Data for PHQ-9 scores

colnames(ds)<-c("Study","a","q1","m","q3","b","n")


dstrue<-read.csv("truevalues.csv")#Include true mean and 95% CI
dstrue$SE<-(dstrue$U-dstrue$L)/(2*1.96)

model<-c()
modelqe<-c()
modelwqe<-c()
modelmde<-c()

for (i in 1:nrow(ds)){
  dat<-ds[i,2:7]
  a<-as.numeric(dat[1])
  q1<-as.numeric(dat[2])
  m<-as.numeric(dat[3])
  q3<-as.numeric(dat[4])
  b<-as.numeric(dat[5])
  n<-as.numeric(dat[6])
  
  
  a1<-as.numeric(dat[1])/27
  q11<-as.numeric(dat[2])/27
  m1<-as.numeric(dat[3])/27
  q31<-as.numeric(dat[4])/27
  b1<-as.numeric(dat[5])/27
  
  
  ds$luo[i]<- (((4/(4+n^0.75))*((a1+b1)/2))+((n^0.75)/(4+n^0.75))*m1)*27
  
  theta1<-2*qnorm((n-0.375)/(n+0.25))
  ds$shi[i] <-((b1-a1)/theta1)*27
  
  ds$mqe[i]<-(as.numeric(qe.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)[1]))*27
  ds$sdqe[i]<-(as.numeric(qe.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)[2]))*27
  ss<-summary(qe.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n))
  modelqe[i]<-rownames(ss)[1]
  
  ds$mwqe[i]<-(as.numeric(wqe.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)[1]))*27
  ds$sdwqe[i]<-(as.numeric(wqe.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)[2]))*27
  swqe<-wqe.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)
  modelwqe[i]<-swqe$selected.dist
  
  
  ds$mmde[i]<-(as.numeric(MDEx.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)[1]))*27
  ds$sdmde[i]<-(as.numeric(MDEx.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)[2]))*27
  smde<-MDEx.mean.sd(min.val = a1,med.val = m1,max.val = b1,n=n)
  modelmde[i]<-smde$selected.dist
}

settings.meta(common = FALSE, details = FALSE,
              digits = 2, digits.tau = 2,
              print.tau2 = FALSE, print.tau2.ci = FALSE, print.H = FALSE,
              print.tau.ci = FALSE, print.I2.ci = FALSE,
              print.Q = FALSE)

mymeta1<-metagen(TE=luo,seTE=shi/sqrt(n),studlab=Study,data=ds,method.tau="REML")
print(mymeta1) #table 2-Luo/Shi(Wan) S1 

mymeta2<-metagen(TE=mqe,seTE=sdqe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
print(mymeta2)#table 2-QE S1 

mymeta3<-metagen(TE=mwqe,seTE=sdwqe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta3)#table 2-wQE S1

mymeta4<-metagen(TE=mmde,seTE=sdmde/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta4) #table 2-MDE S1

mymeta5<-metagen(TE=Mean,seTE=SE,studlab=Study,data=dstrue,method.tau="REML")
summary(mymeta5) #table 2-True sample mean/SD

#Starting value calculation
#User can define their own starting value, or use the following method
#The method calculate the mean and variance using Luo et al. and Shi et al.'s methods
mm<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
theta1<-(2+0.14*(n^0.6))*qnorm((n-0.375)/(n+0.25))
theta2<-(2+(2/(0.07*n^0.6)))*qnorm((0.75*n-0.125)/(n+0.25))
var<-((b-a)/theta1+(q3-q1)/theta2)^2

#starting value for normal distribution 
startmn<-((2.2/(2.2+n^0.75))*((a+b)/2))+((0.7-0.72/(n^0.55))*((q1+q3)/2))+(0.3+0.72/(n^0.55)-(2.2/(2.2+n^0.75)))*m
startsdn<-(b-a)/theta1+(q3-q1)/theta2

#starting value for log-normal distribution
startmln<-log(mm/sqrt(1+var/mm^2))
startsdln<-sqrt(log(1+var/mm^2))

#starting value for Weibull distribution
startbw<-(sqrt(var)/mm)^(-1.086)
startaw<-mm/gamma(1+1/startbw)

#starting value for beta distribution
startab<-mm*(mm*(1-mm)/var-1)
startbb<-(1-mm)*(mm*(1-mm)/var-1)

#starting value for gamma distribution
startbg<-var/mm
startag<-mm/startbg

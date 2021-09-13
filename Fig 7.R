setwd("")
source("functions.R")

iteration=5000
samplesize=c(5,10,14,20,30,40) # number of fish to be sampled
SI=c(1,2,3,4,  7, 10, 14, 21, 28) # sampling interval in days
TT=c(1) # treatment threshold : average number of adult female sea lice per salmon
rho=0.005 # external pressure: No. of floating sea lice per salmon per day
alpha=0.02 # attachment rate: probability of copepodids that attach salmon per day
TE=0.95 # treatment efficacy : proportion of detached planktonic sea lice due to treatment
tfinal=553 # period of a production cycle in days
fish=10000 # number of fish in a pen
kappa=2.19 # dispersion parameter in the negative binomial distribution
pars  <- c(eta=592, epsilon=0.0476, PSU=26.4, Temp=8, alpha=alpha) 
yini  <- c(timetrack=0,ssAB=0,P=1,I=0,C=0,A=0,Ps=1,Is=0,Cs=0,As=0)
times <- seq(0, tfinal, by = 1)

samplingfun <- function(A,samplesize){mean(sample(rnbinom(fish,size=kappa,mu=A),samplesize))}

OP=matrix(NA,length(times)*4*length(TT)*length(samplesize)*length(SI)*iteration); dim(OP)=c(length(times),4,length(TT),length(samplesize),length(SI),iteration)
for(i in 1:iteration){for(w in 1:length(SI)){for(s in 1:length(samplesize)){for(l in 1:length(TT)){
  eventfun <- function(t, y, parms){
    with (as.list(y),{
      timetrack <- timetrack+1
      P <- P
      I <- I
      C <- ifelse(A>=TT[l],C*(1-TE),C)
      A <- ifelse(A>=TT[l],A*(1-TE),A)
      
      Ps <- Ps
      Is <- Is
      Cs <- ifelse(ssAB>=TT[l],Cs*(1-TE),Cs)
      As <- ifelse(ssAB>=TT[l],As*(1-TE),As)
      
      ssAB<-ifelse(timetrack%%SI[w]==0,samplingfun(As,samplesize[s]),0)
      
      return(c(timetrack,ssAB,P,I,C,A,Ps,Is,Cs,As))})}
  
  sealice <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      dtimetrack<-0
      dssAB<-0
      
      dP<-eta*epsilon*V(Temp,PSU)*A-gammaP(Temp)*P-mup(PSU)*P
      dI<-gammaP(Temp)*P-alpha*I-mui(PSU)*I
      dC<-alpha*I-gammaC(Temp)*C-muc(PSU)*C
      dA<-gammaC(Temp)*C/2-mua(PSU)*A
      
      dPs<-eta*epsilon*V(Temp,PSU)*As-gammaP(Temp)*Ps-mup(PSU)*Ps
      dIs<-gammaP(Temp)*Ps-alpha*Is-mui(PSU)*Is
      dCs<-alpha*Is-gammaC(Temp)*Cs-muc(PSU)*Cs
      dAs<-gammaC(Temp)*Cs/2-mua(PSU)*As
      return(list(c(dtimetrack,dssAB,dP,dI,dC,dA,dPs,dIs,dCs,dAs)))})}
  
  out  <- ode(yini, times, sealice, pars, events=list(func=eventfun,time=times))
  OP[,,l,s,w,i]<-out[,c(2,3,7,11)] # time, sampling, realtrue, sampletrue
}}};print(i)}

dim(OP); samplesize; SI

kk=matrix(NA,(tfinal-1)*length(TT)*length(samplesize)*length(SI)*iteration); dim(kk)=c(tfinal-1,length(TT),length(samplesize),length(SI),iteration)
Diff=matrix(NA,length(TT)*length(samplesize)*length(SI)*iteration*7); dim(Diff)=c(length(TT),length(samplesize),length(SI),iteration,7)
for(i in 1:iteration){for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(w in 1:length(SI)){for(t in 1:(tfinal-1)){
  kk[t,l,s,w,i]=ifelse(OP[t,4,l,s,w,i]>OP[t+1,4,l,s,w,i],OP[t,4,l,s,w,i]-TT[l],NA)}
  Diff[l,s,w,i,]=na.omit(kk[,l,s,w,i])[1:7]
}}}}

Summ=matrix(NA,length(TT)*length(SI)*3); dim(Summ)=c(length(TT),length(SI),3)
for(w in 1:length(SI)){for(l in 1:length(TT)){Summ[l,w,]=c(mean(Diff[l,w,w,,],na.rm=T),median(Diff[l,w,w,,],na.rm=T),var(as.numeric(Diff[l,w,w,,]),na.rm=T))}}
Summ=data.frame(Summ[1,,]); colnames(Summ)=c("Mean","Median","Var"); Summ

mDI=melt(Diff); dim(mDI)
colnames(mDI)=c("TT","SS","SI","Iter","Order","Diff"); head(mDI)
mDI[,1]=as.factor(TT[mDI[,1]]); mDI[,2]=as.factor(samplesize[mDI[,2]]); mDI[,3]=as.factor(SI[mDI[,3]])

ct=at=matrix(NA,length(TT)*length(samplesize)*length(SI)); dim(ct)=dim(at)=c(length(TT),length(samplesize),length(SI))
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(w in 1:length(SI)){ct[l,s,w]=median(Diff[l,s,w,,],na.rm = T)
at[l,s,w]=var(as.numeric(Diff[l,s,w,,]),na.rm=T)}}}

filled.contour(cex=0.4,x=samplesize,y=SI,z=ct[1,,],ylim=c(1,20),zlim=c(-0.51,0.4),color=terrain.colors,xlab="Sample Size",ylab="Sampling Interval (days)",key.title=title(main=""),main="",cex.lab=1.5)

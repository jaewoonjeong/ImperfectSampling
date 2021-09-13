setwd("")
source("functions.R")

iteration=5000
samplesize=c(10,30,50,100) # number of fish to be sampled
TT=c(0.2,0.5,1,3) # treatment threshold : average number of adult female sea lice per salmon
rho=0.005 # external pressure: No. of floating sea lice per salmon per day
alpha=0.02 # attachment rate: probability of copepodids that attach salmon per day
TE=0.95 # treatment efficacy : proportion of detached planktonic sea lice due to treatment
SI=7 # sampling interval in days
tfinal=553 # period of a production cycle in days
fish=10000 # number of fish in a pen
kappa=2.19 # dispersion parameter in the negative binomial distribution
pars  <- c(eta=592, epsilon=0.0476, PSU=26.4, Temp=8, alpha=alpha) 
yini  <- c(timetrack=0,ssAB=0,P=1,I=0,C=0,A=0,Ps=1,Is=0,Cs=0,As=0)
times <- seq(0, tfinal, by = 1)

samplingfun <- function(A,samplesize){mean(sample(rnbinom(fish,size=kappa,mu=A),samplesize))}

OP=matrix(NA,length(times)*4*length(TT)*length(samplesize)*iteration); dim(OP)=c(length(times),4,length(TT),length(samplesize),iteration)
for(i in 1:iteration){for(s in 1:length(samplesize)){for(l in 1:length(TT)){
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
      
      ssAB<-ifelse(timetrack%%SI==0,samplingfun(As,samplesize[s]),0)
      
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
  OP[,,l,s,i]<-out[,c(2,3,7,11)] # time, sampling, realtrue, sampletrue
}};print(i)}

Error=matrix(NA,(tfinal/7+1)*length(TT)*length(samplesize)*iteration); dim(Error)=c(tfinal/7+1,length(TT),length(samplesize),iteration)
for(i in 1:iteration){for(t in 1:(tfinal/7+1)){for(l in 1:length(TT)){for(s in 1:length(samplesize)){DF=data.frame(OP[seq(0,tfinal,SI)+1,2,l,s,i],c(0,OP[seq(0,tfinal,SI),4,l,s,i])); colnames(DF)=c("Monitored","True")
Error[t,l,s,i]=ifelse(DF[t,1]>=TT[l] & DF[t,2]<TT[l],"E",ifelse(DF[t,1]<TT[l] & DF[t,2]>=TT[l],"M",ifelse(DF[t,1]>=TT[l] & DF[t,2]>=TT[l], "C",ifelse(DF[t,1]>=TT[l] & DF[t,2]>=TT[l] & DF[t-1,1]<TT[l],"L",NA))))}}};print(i)}
for(t in 2:(tfinal/7+1)){for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){if(is.na(Error[t,l,s,i])==FALSE& is.na(Error[t-1,l,s,i])==FALSE & Error[t,l,s,i]=="C" & Error[t-1,l,s,i]=="M"){Error[t,l,s,i]="L"}}}}}

ErrorCount=matrix(NA,4*length(TT)*length(samplesize)); dim(ErrorCount)=c(length(TT),length(samplesize),4)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){
  ErrorCount[l,s,1]=ifelse(names(table(Error[,l,s,]))[1]=="C",table(Error[,l,s,])[1],0)
  ErrorCount[l,s,2]=ifelse(names(table(Error[,l,s,]))[1]=="E",table(Error[,l,s,])[1],ifelse(names(table(Error[,l,s,]))[2]=="E",table(Error[,l,s,])[2],0))
  ErrorCount[l,s,3]=ifelse(names(table(Error[,l,s,]))[1]=="L",table(Error[,l,s,])[1],ifelse(names(table(Error[,l,s,]))[2]=="L",table(Error[,l,s,])[2],ifelse(names(table(Error[,l,s,]))[3]=="L",table(Error[,l,s,])[3],0)))
  ErrorCount[l,s,4]=ifelse(names(table(Error[,l,s,]))[1]=="M",table(Error[,l,s,])[1],ifelse(names(table(Error[,l,s,]))[2]=="M",table(Error[,l,s,])[2],ifelse(names(table(Error[,l,s,]))[3]=="M",table(Error[,l,s,])[3],ifelse(names(table(Error[,l,s,]))[4]=="M",table(Error[,l,s,])[4],0))))
};print(l)}

Region=RegionOver=matrix(NA,length(TT)*length(samplesize)*iteration); dim(Region)=dim(RegionOver)=c(length(TT),length(samplesize),iteration)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){
  Region[l,s,i]=auc(x=times,y=OP[,4,l,s,i])
  RegionOver[l,s,i]=auc(x=times,y=ifelse(OP[,4,l,s,i]>TT[l],OP[,4,l,s,i],TT[l]))-auc(x=times,y=rep(TT[l],length(times)))
}}}

# difference in estimation of abundance
Chai=AbsChai=SamChai=SamAbsChai=matrix(NA,length(TT)*length(samplesize)*iteration)
dim(Chai)=dim(AbsChai)=dim(SamChai)=dim(SamAbsChai)=c(length(TT),length(samplesize),iteration)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){
  Chai[l,s,i]=mean(OP[seq(SI,tfinal,SI)+1,3,l,s,i]-OP[seq(SI,tfinal,SI)+1,4,l,s,i]) # realtrue - samptrue
  AbsChai[l,s,i]=mean(abs(OP[seq(SI,tfinal,SI)+1,3,l,s,i]-OP[seq(SI,tfinal,SI)+1,4,l,s,i])) # abs(realtrue - samptrue)
  SamChai[l,s,i]=mean(OP[seq(SI,tfinal,SI)+1,4,l,s,i]-OP[seq(SI,tfinal,SI)+1,2,l,s,i]) # samptrue - sampling
  SamAbsChai[l,s,i]=mean(abs(OP[seq(SI,tfinal,SI)+1,4,l,s,i]-OP[seq(SI,tfinal,SI)+1,2,l,s,i])) # abs(samptrue - Sampling)
}}}

AUCofTruesamp=matrix(NA,length(TT)*length(samplesize)*iteration*2); dim(AUCofTruesamp)=c(length(TT),length(samplesize),iteration,2)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){AUCofTruesamp[l,s,i,1]=auc(OP[,1,l,s,i],OP[,4,l,s,i])
AUCofTruesamp[l,s,i,2]=auc(OP[,1,l,s,i],ifelse(OP[,4,l,s,i]<TT[l],TT[l],OP[,4,l,s,i]))-auc(OP[,1,l,s,i],rep(TT[l],length(OP[,1,l,s,i])))}}}

# Difference at TX
DIFF=matrix(NA,80*length(TT)*length(samplesize)*iteration); dim(DIFF)=c(80,length(TT),length(samplesize),iteration)
for(i in 1:iteration){for(t in 1:80){for(l in 1:length(TT)){for(s in 1:length(samplesize)){DIFF[t,l,s,i]=ifelse(is.na(Error[t,l,s,i]) | Error[t,l,s,i]=="M",NA,c(OP[c(seq(0,tfinal,SI)+1)[t],2,l,s,i]-OP[c(seq(0,tfinal,SI)+1)[t],4,l,s,i]))}}};print(i)}

###################### Figure 1
par(mfcol=c(2,1),mar=c(0,4,0.5,1),oma=c(5,1.5,1,0))

for(i in 1){for(l in c(3)){for(s in c(2)){ #
  comparison=data.frame(OP[seq(SI,tfinal,SI)+1,c(2),l,s,i],OP[seq(SI,tfinal,SI),c(4),l,s,i])
  ErrorColor=ifelse(comparison[,1]>=TT[l] & comparison[,2]>=TT[l],6,ifelse(comparison[,1]>=TT[l] & comparison[,2]<TT[l],5,ifelse(comparison[,1]<TT[l] & comparison[,2]>=TT[l],4,ifelse(comparison[,1]>=TT[l] & comparison[,2]>=TT[l] ,7,3))))
  # no tx: 3, missing; 4, early; 5, correct; 6, late: 7
  for(t in 2:dim(comparison)[1]){if(ErrorColor[t]==6 & ErrorColor[t-1]==4){ErrorColor[t]=7}}  
  plot(x=times,y=OP[,4,l,s,i],type="l",main="",xlab="",ylab="Abundance",cex.lab=1.5,xlim=c(0,tfinal),ylim=c(0,TT[l]*2),lty=1,lwd=2,col=1,axes=F,las=1)
  axis(2,las=1,cex.axis=1.5); box(cex.axis=1.5)
  xx=c(times,rev(times)); yy=c(rep(0,length(times)),rev(OP[,4,l,s,i]))
  polygon(xx,yy,col="gray",border=NA,lty=1)
  xx=c(times,rev(times)); yy=c(rep(TT[l],length(times)),rev(ifelse(OP[,4,l,s,i]>TT[l],OP[,4,l,s,i],TT[l])))
  polygon(xx, yy, col = "black", border = "black", lty=1); abline(h=TT[l],lty=1,col=7)
  mtext("(A)",side=2,line=3,at=2,adj=1,cex=2,las=1)
  plot(x=seq(SI,tfinal,SI)+1,y=OP[seq(SI,tfinal,SI)+1,2,l,s,i],type="o",cex.axis=1.5,cex.lab=1.5,
       col=ifelse(ErrorColor==3,"gray",ifelse(ErrorColor==4,"gray",ifelse(ErrorColor==5,4,ifelse(ErrorColor==7,2,3)))),
       pch=19,cex=ifelse(ErrorColor==3,1,ifelse(ErrorColor==4,1,1.5)),main="",xlab="",ylab="Abundance",xlim=c(0,tfinal),las=1,ylim=c(0,TT[l]*2))
  points(x=which(ErrorColor==4)*SI+1,y=OP[which(ErrorColor==4)*SI+1,2,l,s,i],pch=4,col="black")
  abline(h=TT[l],lty=1,col=7) ; mtext(side=1,"Time in days",line=3,cex=1.5) 
  legend("topleft",title="Treatment",legend=c("Early","Timely","Late"),col=c(4,3,2),pch=19,cex=1.5);  par(new=F)
  mtext("(B)",side=2,line=3,at=2,adj=1,cex=2,las=1)
}}}

###################### Figure 3
diff=melt(DIFF); colnames(diff)=c("Time","Threshold","SampleSize","Iteration","Difference")
diff[,3]=as.factor(samplesize[diff[,3]])
diff[,2]=as.factor(ifelse(diff[,2]==1,"0.2",ifelse(diff[,2]==2,"0.5",ifelse(diff[,2]==3,"1","3"))))

p0.2<-ggplot(subset(diff,Threshold==0.2&SampleSize!=50),aes(x=Difference,y=SampleSize,group=SampleSize))+geom_density_ridges(alpha=0.7)+xlim(-1.2,1.2)+xlab("(Monitored Abundance) - (True Abundance)")+ylab("")+geom_vline(xintercept = 0,linetype="dashed",color="grey40",size=1)
p1<-ggplot(subset(diff,Threshold==1&SampleSize!=50),aes(x=Difference,y=SampleSize))+geom_density_ridges()+xlim(-1.2,1.2)+xlab("(Monitored Abundance) - (True Abundance)")+ylab("")+geom_vline(xintercept = 0,linetype="dashed",color="grey40",size=1)
p3<-ggplot(subset(diff,Threshold==3&SampleSize!=50),aes(x=Difference,y=SampleSize))+geom_density_ridges()+xlim(-1.2,1.2)+xlab("(Monitored Abundance) - (True Abundance)")+ylab("")+geom_vline(xintercept = 0,linetype="dashed",color="grey40",size=1)

###################### Figure 5
ee=matrix(NA,length(TT)*length(samplesize)*iteration*4); dim(ee)=c(length(TT),length(samplesize),iteration,4) 
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){ee[l,s,i,]=c(length(which(na.omit(Error[,l,s,i])=="C")),length(which(na.omit(Error[,l,s,i])=="E")),length(which(na.omit(Error[,l,s,i])=="L")),length(which(na.omit(Error[,l,s,i])=="M")))}}}

EEc=melt(ee[,,,1]); EEe=melt(ee[,,,2]); EEl=melt(ee[,,,3]); EEm=melt(ee[,,,4])
EEE=cbind(EEc,EEe[,4],EEl[,4],EEm[,4]) # data frame of error
colnames(EEE)=c("Threshold","SampleSize","Iteration","Timely","Early","Late","Missed")
EEE[,1]=as.factor(rep(TT,length(samplesize)*iteration))
EEE[,2]=as.factor(rep(rep(samplesize,each=length(TT)),iteration))
# TotalData of frequency of tx and AUC
TotalData=cbind(EEE,cbind(melt(AUCofTruesamp[,,,1]),melt(AUCofTruesamp[,,,2])[,4])[,c(4,5)])
colnames(TotalData)=c("Threshold","SampleSize","Iteration","Timely","Early","Late","Missed","AUC","AUCoverTH")
head(TotalData)

library(sqldf)
No=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==0 AND Missed==0')

E1=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==1 AND Late==0 AND Missed==0')
E2=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==2 AND Late==0 AND Missed==0')
E3=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==3 AND Late==0 AND Missed==0')
E4=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==4 AND Late==0 AND Missed==0')
E5=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==5 AND Late==0 AND Missed==0')

L1M1=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==1 AND Missed==1')
L2M2=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==2 AND Missed==2')
L3M3=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==3 AND Missed==3')
L1M2=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==1 AND Missed==2')
L2M3=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==2 AND Missed==3')
L3M4=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==0 AND Late==3 AND Missed==4')

B1=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==1 AND Late==1')
B2=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==2 AND Late==2')
B3=sqldf('SELECT Threshold, SampleSize, Iteration, Timely, Early, Late, AUC, AUCoverTH FROM TotalData WHERE Threshold==1 AND SampleSize==30 AND Early==3 AND Late==3')

par(mfcol=c(2,1),oma=c(1,1,1,1),mar=c(0.3,3,0,1))
boxplot(list(No[,7]/No[1,7],E1[,7]/No[1,7],E2[,7]/No[1,7],E3[,7]/No[1,7],B1[,7]/No[1,7],B2[,7]/No[1,7],B3[,7]/No[1,7],L1M1[,7]/No[1,7],L2M2[,7]/No[1,7],L3M3[,7]/No[1,7],L1M2[,7]/No[1,7],L2M3[,7]/No[1,7],L3M4[,7]/No[1,7]),main="",xlab="",ylab="",las=1,axes=F)
abline(h=1,lty=2); axis(2,las=1); box(); mtext("Relative IB",side=2,line = 3); abline(v=c(1.5,4.5,7.5,10.5),lty=3,col="grey")
boxplot(list(No[,8]/No[1,8],E1[,8]/No[1,8],E2[,8]/No[1,8],E3[,8]/No[1,8],B1[,8]/No[1,8],B2[,8]/No[1,8],B3[,8]/No[1,8],L1M1[,8]/No[1,8],L2M2[,8]/No[1,8],L3M3[,8]/No[1,8],L1M2[,8]/No[1,8],L2M3[,8]/No[1,8],L3M4[,8]/No[1,8]),main="",xlab="",ylab="",las=1,axes=F,ylim=c(0,4.5))
abline(h=1,lty=2); axis(2,las=1); box(); mtext("Relative IB Over Threshold",side=2,line = 3); abline(v=c(1.5,4.5,7.5,10.5),lty=3,col="grey")

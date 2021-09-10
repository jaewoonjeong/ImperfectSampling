setwd("")
source("functions.R")

iteration=1000
samplesize=seq(10,100,10) # number of fish to be sampled
TT=c(0.2,0.5,1,1.5,2,2.5,3) # treatment threshold : average number of adult female sea lice per salmon
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

Error=matrix(NA,(tfinal/7+1)*length(TT)*length(samplesize)*iteration); 
dim(Error)=c(tfinal/7+1,length(TT),length(samplesize),iteration)
for(i in 1:iteration){for(t in 1:(tfinal/7+1)){for(l in 1:length(TT)){for(s in 1:length(samplesize)){
  DF=data.frame(OP[seq(0,tfinal,SI)+1,2,l,s,i],c(0,OP[seq(0,tfinal,SI),4,l,s,i])); colnames(DF)=c("Monitored","True")
  Error[t,l,s,i]=ifelse(DF[t,1]>=TT[l] & DF[t,2]<TT[l],"E",ifelse(DF[t,1]<TT[l] & DF[t,2]>=TT[l],"M",ifelse(DF[t,1]>=TT[l] & DF[t,2]>=TT[l], "C",ifelse(DF[t,1]>=TT[l] & DF[t,2]>=TT[l] & DF[t-1,1]<TT[l],"L",NA))))
}}};print(i)}

for(t in 2:(tfinal/7+1)){for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){
  if(is.na(Error[t,l,s,i])==FALSE& is.na(Error[t-1,l,s,i])==FALSE & Error[t,l,s,i]=="C" & Error[t-1,l,s,i]=="M"){Error[t,l,s,i]="L"}}}}}

ErrorCount=matrix(NA,4*length(TT)*length(samplesize)); dim(ErrorCount)=c(length(TT),length(samplesize),4)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){
  ErrorCount[l,s,1]=ifelse(names(table(Error[,l,s,]))[1]=="C",table(Error[,l,s,])[1],0)
  ErrorCount[l,s,2]=ifelse(names(table(Error[,l,s,]))[1]=="E",table(Error[,l,s,])[1],ifelse(names(table(Error[,l,s,]))[2]=="E",table(Error[,l,s,])[2],0))
  ErrorCount[l,s,3]=ifelse(names(table(Error[,l,s,]))[1]=="L",table(Error[,l,s,])[1],ifelse(names(table(Error[,l,s,]))[2]=="L",table(Error[,l,s,])[2],ifelse(names(table(Error[,l,s,]))[3]=="L",table(Error[,l,s,])[3],0)))
  ErrorCount[l,s,4]=ifelse(names(table(Error[,l,s,]))[1]=="M",table(Error[,l,s,])[1],ifelse(names(table(Error[,l,s,]))[2]=="M",table(Error[,l,s,])[2],ifelse(names(table(Error[,l,s,]))[3]=="M",table(Error[,l,s,])[3],ifelse(names(table(Error[,l,s,]))[4]=="M",table(Error[,l,s,])[4],0))))
};print(l)}

# Difference at TX
DIFF=matrix(NA,80*length(TT)*length(samplesize)*iteration); dim(DIFF)=c(80,length(TT),length(samplesize),iteration)
for(i in 1:iteration){for(t in 1:80){for(l in 1:length(TT)){for(s in 1:length(samplesize)){
  DIFF[t,l,s,i]=ifelse(is.na(Error[t,l,s,i]),NA,c(OP[c(seq(0,tfinal,SI)+1)[t],2,l,s,i]-OP[c(seq(0,tfinal,SI)+1)[t],4,l,s,i]))
}}};print(i)}

Tr=matrix(NA,length(1:tfinal)*length(TT)*length(samplesize)*iteration); dim(Tr)=c(length(1:tfinal),length(TT),length(samplesize),iteration)
for(i in 1:iteration){for(s in 1:length(samplesize)){for(l in 1:length(TT)){for(t in 2:tfinal){Tr[t,l,s,i]=ifelse(OP[t-1,4,l,s,i]>OP[t,4,l,s,i],"T",NA)}}}}

########################## Figure 2
ErrorProp=matrix(NA,3*length(samplesize)*length(TT)); dim(ErrorProp)=c(length(TT),length(samplesize),3)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){ErrorProp[l,s,]=ErrorCount[l,s,]/sum(ErrorCount[l,s,])}}
EP=melt(ErrorProp); colnames(EP)=c("Threshold","SampleSize","Treatment","Proportion")
EP[,1]=as.factor(paste("Threshold =",TT[EP[,1]]))
EP[,2]=as.factor(rep(samplesize,each=length(TT)))
EP[,3]=as.factor(ifelse(EP[,3]==1,"B",ifelse(EP[,3]==2,"A","C")))
EP=subset(EP,Threshold!="Threshold = 2"&Threshold!="Threshold = 2.5"&Threshold!="Threshold = 1.5")

ggplot(EP,aes(fill=Treatment,y=Proportion,x=SampleSize))+geom_bar(position="fill", stat="identity")+
  facet_grid(Threshold~.)+theme(legend.position="bottom")+xlab("Sample Size")+
  scale_fill_manual(values = c("A" = "Blue", "B" = "green3", "C" = "red"),name = "Treatment", labels = c("Early Treatment", "Timely Treatment", "Late Treatment"))+
  guides(fill=guide_legend(title=""))+ theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16))

################################## Figure 4
ee=matrix(NA,length(TT)*length(samplesize)*iteration*4); dim(ee)=c(length(TT),length(samplesize),iteration,4) 
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){ee[l,s,i,]=c(length(which(na.omit(Error[,l,s,i])=="C")),length(which(na.omit(Error[,l,s,i])=="E")),length(which(na.omit(Error[,l,s,i])=="L")),length(which(na.omit(Error[,l,s,i])=="M")))}}}
TxFreq=matrix(NA,length(TT)*length(samplesize)); dim(TxFreq)=c(length(TT),length(samplesize))
TF=matrix(NA,length(TT)*length(samplesize)*iteration); dim(TF)=c(length(TT),length(samplesize),iteration)
for(l in 1:length(TT)){for(s in 1:length(samplesize)){TF[l,s,]=rowSums(ee[l,s,,-4])}}
for(l in 1:length(TT)){for(s in 1:length(samplesize)){TxFreq[l,s]=mean(rowSums(ee[l,s,,-4]))}}
IB=IBOT=matrix(NA,length(TT)*length(samplesize)*iteration); dim(IB)=dim(IBOT)=c(length(TT),length(samplesize),iteration) 
mIB=mIBOT=matrix(NA,length(TT)*length(samplesize)); dim(mIB)=dim(mIBOT)=c(length(TT),length(samplesize)) 
for(l in 1:length(TT)){for(s in 1:length(samplesize)){for(i in 1:iteration){IB[l,s,i]=auc(OP[,1,l,s,i],OP[,4,l,s,i]); IBOT[l,s,i]=auc(OP[,1,l,s,i],ifelse(OP[,4,l,s,i]<TT[l],0,OP[,4,l,s,i]))}; mIB[l,s]=mean(IB[l,s,]); mIBOT[l,s]=mean(IBOT[l,s,])}}

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)); par(mar=c(4,5.5,4,13),xpd=TRUE,cex.lab=1.3); colfunc <- colorRampPalette(c("black", "gray"))
matplot(samplesize,t(TxFreq[c(1,2,3,7),]),cex.main=1.5,cex.lab=1.5,cex.axis=1.5,type="l",pch=1,las=1,main="(A)",xlab="Sample Size",lwd=3,ylab="Frequency of Treatments",lty=1,col=colfunc(4),ylim=c(5.5,7.5))
legend("bottomright",bty='n',inset=c(-0.4,0),legend=TT[c(1,2,3,7)],title="Treatment Threshold",col=colfunc(4),lty=1,lwd=3,cex=1.5)
par(mar=c(5,5,4,2))
matplot(samplesize,t(mIB[c(1,2,3,7),]),type="l",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,pch=1,las=1,xlab="Sample Size",main="(B)",ylab="IB",col=colfunc(4),ylim=c(10,600),lty=1,log='y',lwd=3)
matplot(samplesize,t(mIBOT[c(1,2,3,7),]),type="l",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,pch=2,las=1,xlab="Sample Size",main="(C)",ylab="IB Over Threshold",col=colfunc(4),lty=1,ylim=c(1,100),log='y',lwd=3)

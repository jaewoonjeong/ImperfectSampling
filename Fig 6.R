setwd("")
source("functions.R")

iteration=1
samplesize=c(20) # number of fish to be sampled
SI=c(2,6,18) # sampling interval in days
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

l=1; 
par(mfcol=c(2,3),mar=c(0,1,0.5,1),oma=c(5,5,2.5,0),cex.lab=1.3,cex.axis=1.3)
for(i in 1){for(s in c(1)){for(v in 1:3){
  comparison=data.frame(OP[seq(SI[v],tfinal,SI[v])+1,c(2),1,s,v,i],OP[seq(SI[v],tfinal,SI[v]),c(4),1,s,v,i])
  ErrorColor=ifelse(comparison[,1]>=TT[l] & comparison[,2]>=TT[l],6,ifelse(comparison[,1]>=TT[l]  & comparison[,2]<TT[l],5,ifelse(comparison[,1]<TT[l] & comparison[,2]>=TT[l],4,ifelse(comparison[,1]>=TT[l]  & comparison[,2]>=TT[l] ,7,3))))
  for(t in 2:dim(comparison)[1]){if(ErrorColor[t]==6 & ErrorColor[t-1]==4){ErrorColor[t]=7}}  
  plot(x=times,y=OP[,4,1,s,v,i],type="l",xlab="",ylab=ifelse(s==1,"Abundance",""),xlim=c(0,550),ylim=c(0,2.2),main="",lty=1,lwd=2,col=1,axes=F,las=1,cex.lab=1.5,cex.axis=1.5); abline(h=1,col="yellow")
  mtext(paste("(",LETTERS[v],")"),lwd=2,side=3,line=1,cex=1.5)
  #legend("topleft",cex=1.3,bty='n',legend=c(paste("Sample Size =",samplesize[s]),paste("Sampling interval =",SI[v],"days")))
  ifelse(v==1,axis(2,las=1)&mtext("Abundance",2,line=3),NA); box()
  xx=c(times,rev(times)); yy=c(rep(0,length(times)),rev(OP[,4,1,s,v,i]))
  polygon(xx,yy,col="gray",border=NA,lty=1)
  xx=c(times,rev(times)); yy=c(rep(TT[l],length(times)),rev(ifelse(OP[,4,1,s,v,i]>TT[l],OP[,4,1,s,v,i],TT[l])))
  polygon(xx, yy, col = "black", border = "black", lty=1); abline(h=TT[l],lty=1,col=7)
  plot(main="",x=seq(SI[v],tfinal,SI[v])+1,y=OP[seq(SI[v],tfinal,SI[v])+1,2,1,s,v,i],xlab="",cex=ifelse(ErrorColor==3,1,ifelse(ErrorColor==4,1,1.5)),pch=19,ylab=ifelse(s==1,"Abundance",""),type="o",xlim=c(0,550),ylim=c(0,2.6),las=1,axes=F,
       col=ifelse(ErrorColor==3,"gray",ifelse(ErrorColor==4,"gray",ifelse(ErrorColor==5,4,ifelse(ErrorColor==7,2,3)))))
  #if(sum(which(ErrorColor==4))!=0){points(x=which(ErrorColor==4)*SI[v]+1,y=OP[which(ErrorColor==4)*SI[v]+1,2,l,s,t,i],pch=4,col="black")}
  abline(h=1,col="yellow"); ifelse(v==1,axis(2,las=1)&mtext("Abundance",2,line=3),NA); box()
  mtext(side=1,"Time in days",line=3) ; axis(1)
  legend("topleft",title="Treatment",legend=c("Early","Timely","Late"),col=c(4,3,2),pch=19,cex=1.5)
}}}

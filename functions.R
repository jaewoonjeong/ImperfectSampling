V<-function(T,S){V=NULL
Tvec = T
Svec=S
for(i in seq(1,length(Tvec))){T = Tvec[i]
S = Svec[i]
V = rbind(V,-2.2579 + 0.4557*log(T)+0.6373*log(S))}
V[V<0]<-0
V[V>1]<-1
return(V)}

mup <- function(S) {
  Svec = S
  mup = rep(0,length(Svec))
  for (i in seq(1,length(Svec))){S = Svec[i]
  if(S<=25){mup[i] = 2.16} else if(S>25){mup[i] = exp(13.64-0.515*S)}} 
  return(mup)}
mui = function(S){exp(-.0232-.059*S)}
muc = function(Sal){muc = mua(Sal)
return(muc)}
mua = function(Sal){sigma0a = 4.12
sigma1a = 0.124
mua = exp(((sigma0a + sigma1a * Sal)))/24
mua[mua>1963/24]=1963/24
fmua = 1/mua
return(fmua)}

gammaP = function(T){b0 = 1.216
b1 = -1.38
A = exp(b0+b1*log(T/10))
A[A>23.5]=23.5
gammaP = 1/A
return(gammaP)}

gammaC = function(Temp){Tref = 10
alphac = 11.94
Betac = 0.177
gC = ((Betac * (Temp - Tref + alphac)) / (alphac)) ^ 2
gC[gC<((Betac * (8 - Tref + alphac)) / (alphac)) ^ 2]=((Betac * (8 - Tref + alphac)) / (alphac)) ^ 2
return(gC)}


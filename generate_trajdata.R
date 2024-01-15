generate_group_data <- function(n) {
  set.seed(956)
  J=3 #number of groups
  time=1:6 #Folow-up period
  t=6 #Number of periods
  props_c=c(0.3,0.25,0.45)  #proportion in each group
  nper_group= props_c*n #True sample size ine each group

  #Fonction expit
  expit <- plogis;
  G1 <- rep(1.42,5)
  G2 <- rep(-0.67,5)
  G3 <-rep(0.83,5)

  beta <- rbind(G1,G2,G3)
  B <- 1000
  #Initialization
  sd0<-c(1.75,0.75,2.5)
  V <- rnorm(n)
  probsL0<-c(0.75,0.1,0.001)
  L0<-unlist(sapply(1:J,function(x)rbinom(nper_group[x],1,probsL0[x])))
  probsA0<-c(1,0.5,0.001)
  A0<-unlist(sapply(1:J,function(x)rbinom(nper_group[x],1,probsA0[x])))

  A<-matrix(0,n,t)
  L<-matrix(0,n,t)

  A[,1]<-A0
  L[,1]<-L0

  #Matrix of parameters
  True_param_indiv_temp<-do.call(rbind,lapply(1:J,function(x)
    matrix(rep(beta[x,],nper_group[x]),nrow=nper_group[x],ncol= 5,byrow = TRUE)))

  TPI<-data.frame(True_param_indiv_temp)
  TPI$Group<-unlist(sapply(1:J,function(x)rep(x,nper_group[x])))

  #Generate the exposure and covariates
  for(i in 2:t){
    L[,i]<-rbinom(n,1,expit(0.5*L[,i-1]+0.25*A[,i-1] + V))
    A[,i]<-sapply(1:n,function(x)rbinom(1,1,expit(TPI[x,1]+TPI[x,2]*L[x,i]+TPI[x,3]*L[x,i]^2+TPI[x,4]*A[x,i-1]*L[x,i]+TPI[x,5]*A[x,i-1] + V[x])))
  }

  #True proportion per group
  dat.rshA<-reshape(data.frame(A),varying=time, direction="long",sep="")
  colnames(dat.rshA)<-c("time","A","id")
  datA<-data.frame(A)
  datA$Group<-TPI$Group

  dat.rshL<-reshape(data.frame(L),varying=time, direction="long",sep="")
  colnames(dat.rshL)<-c("time","L","id")

  dat.rshA <-dat.rshA[order(dat.rshA$time),]
  dat.rshL <-dat.rshL[order(dat.rshL$time),]
  dat.rshAL <- dat.rshA
  dat.rshAL$L <- dat.rshL$L
  dat.rshAL$V <- V

  return(dat.rshAL)
}

# Example usage
# dat.rshAL <- generate_group_data()

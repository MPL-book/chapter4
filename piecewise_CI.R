summary(m1.tstage.int,se="M2QM2")

#log of baseline hazard ------

#get times
time=predict(m1.tstage,type = "hazard",cov=c(0,0,0,0,0,0,0))$time
theta=m1.tstage$theta
covar_theta=m1.tstage$covar$M2QM2[c(16:21),c(16:21)]

psi=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=1)
h_0=psi%*%theta
logh0=log(h_0)
LL=rep(NA,1000)
UL=rep(NA,1000)

for(t in 1:1000){
  d_logh0=as.matrix(psi[t,])%*%as.matrix(1/h_0[t])
  var_logh0=t(d_logh0)%*%covar_theta%*%d_logh0
  se_logh0=sqrt(var_logh0)
  
  LL[t]=exp(logh0[t]-1.96*se_logh0)
  UL[t]=exp(logh0[t]+1.96*se_logh0)
  
}

plot(time,UL,type="l",col="grey",lty=2,ylim=c(0,0.5),xlim=c(0,25),ylab="Baseline hazard function",xlab="Time (years)",
     main="Estimate of the baseline hazard function")
points(time,h_0,type="l")

points(time,LL,type="l",col="grey",lty=2)





pop.survival.t1a[415]

Z.t1a=matrix(c(1,0,apply(m1.tstage$data$Z,2,mean)[-c(1:2)]),nrow=1)
ZB.t1a=Z.t1a%*%m1.tstage$beta
pi.t1a=exp(ZB.t1a)/(1+exp(ZB.t1a))
survival.t1a=predict(m1.tstage,type = "survival",cov=c(0,apply(m1.tstage$data$X[,-1],2,mean)))$survival
survival.t1a_t=survival.t1a[415]
X.t1a=matrix(c(0,apply(m1.tstage$data$X[,-1],2,mean)),nrow=1)
time=predict(m1.tstage,type = "survival",cov=c(0,0,0,0,0,0,0))$time[415]
Psi_t=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
mu_t=exp(X.t1a%*%m1.tstage$gamma)

dSpop_dBeta=-pi.t1a%*%(1-pi.t1a)%*%(1-survival.t1a_t)%*%Z.t1a
dSpop_dGamma=-pi.t1a%*%survival.t1a_t%*%(-log(survival.t1a_t))%*%X.t1a
dSpop_dTheta=-pi.t1a%*%survival.t1a_t%*%mu_t%*%Psi_t

dSpop_dEta=matrix(c(dSpop_dBeta,dSpop_dGamma,dSpop_dTheta),nrow=1)

covar_eta=m1.tstage$covar$M2QM2

var_Spop.t1a=dSpop_dEta%*%covar_eta%*%t(dSpop_dEta)
se_Spop.t1a=sqrt(var_Spop.t1a)

pop.survival.t1a[415]
pop.survival.t1a[415]-1.96*se_Spop.t1a
pop.survival.t1a[415]+1.96*se_Spop.t1a




Z.t1b=matrix(c(1,1,apply(m1.tstage$data$Z,2,mean)[-c(1:2)]),nrow=1)
ZB.t1b=Z.t1a%*%m1.tstage$beta
pi.t1b=exp(ZB.t1b)/(1+exp(ZB.t1b))
survival.t1b=predict(m1.tstage,type = "survival",cov=c(1,apply(m1.tstage$data$X[,-1],2,mean)))$survival
survival.t1b_t=survival.t1b[415]
X.t1b=matrix(c(1,apply(m1.tstage$data$X[,-1],2,mean)),nrow=1)
time=predict(m1.tstage,type = "survival",cov=c(0,0,0,0,0,0,0))$time[415]
Psi_t=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
mu_t=exp(X.t1b%*%m1.tstage$gamma)

dSpop_dBeta=-pi.t1b%*%(1-pi.t1b)%*%(1-survival.t1b_t)%*%Z.t1b
dSpop_dGamma=-pi.t1b%*%survival.t1b_t%*%(-log(survival.t1b_t))%*%X.t1b
dSpop_dTheta=-pi.t1b%*%survival.t1b_t%*%mu_t%*%Psi_t

dSpop_dEta=matrix(c(dSpop_dBeta,dSpop_dGamma,dSpop_dTheta),nrow=1)

covar_eta=m1.tstage$covar$M2QM2

var_Spop.t1b=dSpop_dEta%*%covar_eta%*%t(dSpop_dEta)
se_Spop.t1b=sqrt(var_Spop.t1b)


pop.survival.t1b[415]
pop.survival.t1b[415]-1.96*se_Spop.t1b
pop.survival.t1b[415]+1.96*se_Spop.t1b



#MALE VS FEMALE --------


Z.male=matrix(c(apply(m1.tstage$data$Z[,c(1:2)],2,mean),1,apply(m1.tstage$data$Z[,-c(1:3)],2,mean)),nrow=1)
ZB.male=Z.male%*%m1.tstage$beta
pi.male=exp(ZB.male)/(1+exp(ZB.male))
survival.male=predict(m1.tstage,type = "survival",cov=c(mean(m1.tstage$data$X[,1]),1,apply(m1.tstage$data$X[,-c(1:2)],2,mean)))$survival
pop.survival.male=as.numeric(pi.male)*as.numeric(survival.male) + 1-as.numeric(pi.male)
survival.male_t=survival.male[415]
X.male=matrix(c(mean(m1.tstage$data$X[,1]),1,apply(m1.tstage$data$X[,-c(1:2)],2,mean)),nrow=1)
time=predict(m1.tstage,type = "survival",cov=c(0,0,0,0,0,0,0))$time[415]
Psi_t=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
mu_t=exp(X.male%*%m1.tstage$gamma)

dSpop_dBeta=-pi.male%*%(1-pi.male)%*%(1-survival.male_t)%*%Z.male
dSpop_dGamma=-pi.male%*%survival.male_t%*%(-log(survival.male_t))%*%X.male
dSpop_dTheta=-pi.male%*%survival.male_t%*%mu_t%*%Psi_t

dSpop_dEta=matrix(c(dSpop_dBeta,dSpop_dGamma,dSpop_dTheta),nrow=1)

covar_eta=m1.tstage$covar$M2QM2

var_Spop.male=dSpop_dEta%*%covar_eta%*%t(dSpop_dEta)
se_Spop.male=sqrt(var_Spop.male)

pop.survival.male[415]
pop.survival.male[415]-1.96*se_Spop.male
pop.survival.male[415]+1.96*se_Spop.male


Z.female=matrix(c(apply(m1.tstage$data$Z[,c(1:2)],2,mean),0,apply(m1.tstage$data$Z[,-c(1:3)],2,mean)),nrow=1)
ZB.female=Z.female%*%m1.tstage$beta
pi.female=exp(ZB.female)/(1+exp(ZB.female))
survival.female=predict(m1.tstage,type = "survival",cov=c(mean(m1.tstage$data$X[,1]),0,apply(m1.tstage$data$X[,-c(1:2)],2,mean)))$survival
pop.survival.female=as.numeric(pi.female)*as.numeric(survival.female) + 1-as.numeric(pi.female)
survival.female_t=survival.female[415]
X.female=matrix(c(mean(m1.tstage$data$X[,1]),0,apply(m1.tstage$data$X[,-c(1:2)],2,mean)),nrow=1)
time=predict(m1.tstage,type = "survival",cov=c(0,0,0,0,0,0,0))$time[415]
Psi_t=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
mu_t=exp(X.female%*%m1.tstage$gamma)

dSpop_dBeta=-pi.female%*%(1-pi.female)%*%(1-survival.female_t)%*%Z.female
dSpop_dGamma=-pi.female%*%survival.female_t%*%(-log(survival.female_t))%*%X.female
dSpop_dTheta=-pi.female%*%survival.female_t%*%mu_t%*%Psi_t

dSpop_dEta=matrix(c(dSpop_dBeta,dSpop_dGamma,dSpop_dTheta),nrow=1)

covar_eta=m1.tstage$covar$M2QM2

var_Spop.female=dSpop_dEta%*%covar_eta%*%t(dSpop_dEta)
se_Spop.female=sqrt(var_Spop.female)

pop.survival.female[415]
pop.survival.female[415]-1.96*se_Spop.female
pop.survival.female[415]+1.96*se_Spop.female






#logit...... -------

time=predict(m1.tstage,type = "survival",cov=c(0,apply(m1.tstage$data$X[,-1],2,mean)))$time
survival=predict(m1.tstage,type = "survival",cov=c(0,apply(m1.tstage$data$X[,-1],2,mean)))$survival


Z.t1a=matrix(c(apply(m1.tstage$data$Z,2,mean)[1],0,apply(m1.tstage$data$Z,2,mean)[-c(1:2)]),nrow=1)
ZB.t1a=Z.t1a%*%m1.tstage$beta
pi.t1a=exp(ZB.t1a)/(1+exp(ZB.t1a))

pop.survival.t1a=as.numeric(pi.t1a)*as.numeric(survival) + 1-as.numeric(pi.t1a)

plot(time,pop.survival.t1a,type="l")

X.t1a=matrix(c(0,apply(m1.tstage$data$X[,-1],2,mean)),nrow=1)
mu_X.t1a=exp(X.t1a%*%m1.tstage$gamma)

cov_eta=m1.tstage$covar$M2QM2

UL_population.survival.t1a.logit=rep(0,length(time))

LL_population.survival.t1a.logit=rep(0,length(time))

for(p in 2:length(time)){
  t=time[p]
  St=survival[p]
  logitSpop=log(pop.survival.t1a[p]/(1-pop.survival.t1a[p]))
  Psi_t=basis_phmc(t,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
  
  dlogitSt_dbeta=solve(pop.survival.t1a[p]*(1-pop.survival.t1a[p]))%*%(St-1)%*%pi.t1a%*%(1-pi.t1a)%*%Z.t1a
  
  dlogitSt_dgamma=-solve(pop.survival.t1a[p]*(1-pop.survival.t1a[p]))%*%(pi.t1a%*%St%*%(-log(St)))%*%X.t1a
  
  dlogitSt_dtheta=-mu_X.t1a%*%pi.t1a%*%St%*%solve(pop.survival.t1a[p]*(1-pop.survival.t1a[p]))%*%Psi_t
  
  dlogitSt_deta=matrix(c(dlogitSt_dbeta,dlogitSt_dgamma,dlogitSt_dtheta),nrow=1)
  
  var_logitSt=dlogitSt_deta%*%cov_eta%*%t(dlogitSt_deta)
  se_logitSt=sqrt(var_logitSt)
  
  
  
  LL_population.survival.t1a.logit[p]=logitSpop-1.96*se_logitSt
  UL_population.survival.t1a.logit[p]=logitSpop+1.96*se_logitSt
  
}

LL=exp(LL_population.survival.t1a.logit)/(1+exp(LL_population.survival.t1a.logit))
UL=exp(UL_population.survival.t1a.logit)/(1+exp(UL_population.survival.t1a.logit))



plot(time,pop.survival.t1a,ylim=c(0.6,1),type="l",main="Estimated mixture survival (T1a vs T1b)",
     ylab="Estimated survival function",xlab="Time (years)")
points(c(time[1],time[-c(1:16)]),c(pop.survival.t1a[1],LL[-c(1:16)]),type="l",lty=2,col="grey")
points(c(time[1],time[-c(1:16)]),c(pop.survival.t1a[1],UL[-c(1:16)]),type="l",lty=2,col="grey")



survival=predict(m1.tstage,type = "survival",cov=c(1,apply(m1.tstage$data$X[,-1],2,mean)))$survival


Z.t1b=matrix(c(1,1,apply(m1.tstage$data$Z,2,mean)[-c(1:2)]),nrow=1)
ZB.t1b=Z.t1b%*%m1.tstage$beta
pi.t1b=exp(ZB.t1b)/(1+exp(ZB.t1b))

pop.survival.t1b=as.numeric(pi.t1b)*as.numeric(survival) + 1-as.numeric(pi.t1b)


X.t1b=matrix(c(1,apply(m1.tstage$data$X[,-1],2,mean)),nrow=1)
mu_X.t1b=exp(X.t1b%*%m1.tstage$gamma)

cov_eta=m1.tstage$covar$M2QM2

UL_population.survival.t1b.logit=rep(0,length(time))

LL_population.survival.t1b.logit=rep(0,length(time))



for(p in 2:length(time)){
  t=time[p]
  St=survival[p]
  logitSpop=log(pop.survival.t1b[p]/(1-pop.survival.t1b[p]))
  Psi_t=basis_phmc(t,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
  
  dlogitSt_dbeta=solve(pop.survival.t1b[p]*(1-pop.survival.t1b[p]))%*%(St-1)%*%pi.t1b%*%(1-pi.t1b)%*%Z.t1b
  
  dlogitSt_dgamma=-solve(pop.survival.t1b[p]*(1-pop.survival.t1b[p]))%*%(pi.t1b%*%St%*%(-log(St)))%*%X.t1b
  
  dlogitSt_dtheta=-mu_X.t1b%*%pi.t1b%*%St%*%solve(pop.survival.t1b[p]*(1-pop.survival.t1b[p]))%*%Psi_t
  
  dlogitSt_deta=matrix(c(dlogitSt_dbeta,dlogitSt_dgamma,dlogitSt_dtheta),nrow=1)
  
  var_logitSt=dlogitSt_deta%*%cov_eta%*%t(dlogitSt_deta)
  se_logitSt=sqrt(var_logitSt)
  
  
  
  LL_population.survival.t1b.logit[p]=logitSpop-1.96*se_logitSt
  UL_population.survival.t1b.logit[p]=logitSpop+1.96*se_logitSt
  
}

LL=exp(LL_population.survival.t1b.logit)/(1+exp(LL_population.survival.t1b.logit))
UL=exp(UL_population.survival.t1b.logit)/(1+exp(UL_population.survival.t1b.logit))

points(time,pop.survival.t1b,ylim=c(0.6,1),type="l",col="red")
points(c(time[1],time[-c(1:16)]),c(pop.survival.t1b[1],LL[-c(1:16)]),type="l",lty=2,col="pink")
points(c(time[1],time[-c(1:16)]),c(pop.survival.t1b[1],UL[-c(1:16)]),type="l",lty=2,col="pink")

legend(x="bottomright",legend=c("Stage 1a","Stage 1b"),col = c("black","red"),lty=1,bty="n",cex=0.75)



##logit for male vs female



time=predict(m1.tstage,type = "survival",cov=c(0,apply(m1.tstage$data$X[,-1],2,mean)))$time
survival=predict(m1.tstage,type = "survival",cov=c(mean(m1.tstage$data$X[,1]),1,apply(m1.tstage$data$X[,-c(1:2)],2,mean)))$survival

Z.male=matrix(c(apply(m1.tstage$data$Z[,c(1:2)],2,mean),1,apply(m1.tstage$data$Z[,-c(1:3)],2,mean)),nrow=1)
ZB.male=Z.male%*%m1.tstage$beta
pi.male=exp(ZB.male)/(1+exp(ZB.male))

pop.survival.male=as.numeric(pi.male)*as.numeric(survival) + 1-as.numeric(pi.male)


X.male=matrix(c(mean(m1.tstage$data$X[,1]),1,apply(m1.tstage$data$X[,-c(1:2)],2,mean)),nrow=1)
mu_X.male=exp(X.male%*%m1.tstage$gamma)

cov_eta=m1.tstage$covar$M2QM2

UL_population.survival.male.logit=rep(0,length(time))

LL_population.survival.male.logit=rep(0,length(time))


for(p in 2:length(time)){
  t=time[p]
  St=survival[p]
  logitSpop=log(pop.survival.male[p]/(1-pop.survival.male[p]))
  Psi_t=basis_phmc(t,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
  
  dlogitSt_dbeta=solve(pop.survival.male[p]*(1-pop.survival.male[p]))%*%(St-1)%*%pi.male%*%(1-pi.male)%*%Z.male
  
  dlogitSt_dgamma=-solve(pop.survival.male[p]*(1-pop.survival.male[p]))%*%(pi.male%*%St%*%(-log(St)))%*%X.male
  
  dlogitSt_dtheta=-mu_X.male%*%pi.male%*%St%*%solve(pop.survival.male[p]*(1-pop.survival.male[p]))%*%Psi_t
  
  dlogitSt_deta=matrix(c(dlogitSt_dbeta,dlogitSt_dgamma,dlogitSt_dtheta),nrow=1)
  
  var_logitSt=dlogitSt_deta%*%cov_eta%*%t(dlogitSt_deta)
  se_logitSt=sqrt(var_logitSt)
  
  
  
  LL_population.survival.male.logit[p]=logitSpop-1.96*se_logitSt
  UL_population.survival.male.logit[p]=logitSpop+1.96*se_logitSt
  
}

LL=exp(LL_population.survival.male.logit)/(1+exp(LL_population.survival.male.logit))
UL=exp(UL_population.survival.male.logit)/(1+exp(UL_population.survival.male.logit))
pop.survival.male[415]
LL[415]
UL[415]

plot(time,pop.survival.male,ylim=c(0.6,1),type="l",main="Estimated mixture survival (male vs female)",
     ylab="Estimated survival function",xlab="Time (years)")
points(c(time[1],time[-c(1:16)]),c(pop.survival.male[1],LL[-c(1:16)]),type="l",lty=2,col="grey")
points(c(time[1],time[-c(1:16)]),c(pop.survival.male[1],UL[-c(1:16)]),type="l",lty=2,col="grey")


survival=predict(m1.tstage,type = "survival",cov=c(mean(m1.tstage$data$X[,1]),0,apply(m1.tstage$data$X[,-c(1:2)],2,mean)))$survival

Z.female=matrix(c(apply(m1.tstage$data$Z[,c(1:2)],2,mean),0,apply(m1.tstage$data$Z[,-c(1:3)],2,mean)),nrow=1)
ZB.female=Z.female%*%m1.tstage$beta
pi.female=exp(ZB.female)/(1+exp(ZB.female))

pop.survival.female=as.numeric(pi.female)*as.numeric(survival) + 1-as.numeric(pi.female)


X.female=matrix(c(mean(m1.tstage$data$X[,1]),0,apply(m1.tstage$data$X[,-c(1:2)],2,mean)),nrow=1)
mu_X.female=exp(X.female%*%m1.tstage$gamma)

cov_eta=m1.tstage$covar$M2QM2

UL_population.survival.female.logit=rep(0,length(time))

LL_population.survival.female.logit=rep(0,length(time))


for(p in 2:length(time)){
  t=time[p]
  St=survival[p]
  logitSpop=log(pop.survival.female[p]/(1-pop.survival.female[p]))
  Psi_t=basis_phmc(t,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
  
  dlogitSt_dbeta=solve(pop.survival.female[p]*(1-pop.survival.female[p]))%*%(St-1)%*%pi.female%*%(1-pi.female)%*%Z.female
  
  dlogitSt_dgamma=-solve(pop.survival.female[p]*(1-pop.survival.female[p]))%*%(pi.female%*%St%*%(-log(St)))%*%X.female
  
  dlogitSt_dtheta=-mu_X.female%*%pi.female%*%St%*%solve(pop.survival.female[p]*(1-pop.survival.female[p]))%*%Psi_t
  
  dlogitSt_deta=matrix(c(dlogitSt_dbeta,dlogitSt_dgamma,dlogitSt_dtheta),nrow=1)
  
  var_logitSt=dlogitSt_deta%*%cov_eta%*%t(dlogitSt_deta)
  se_logitSt=sqrt(var_logitSt)
  
  
  
  LL_population.survival.female.logit[p]=logitSpop-1.96*se_logitSt
  UL_population.survival.female.logit[p]=logitSpop+1.96*se_logitSt
  
}

LL=exp(LL_population.survival.female.logit)/(1+exp(LL_population.survival.female.logit))
UL=exp(UL_population.survival.female.logit)/(1+exp(UL_population.survival.female.logit))
pop.survival.female[415]
LL[415]
UL[415]

points(time,pop.survival.female,type="l",col="red")
points(c(time[1],time[-c(1:16)]),c(pop.survival.female[1],LL[-c(1:16)]),type="l",lty=2,col="pink")
points(c(time[1],time[-c(1:16)]),c(pop.survival.female[1],UL[-c(1:16)]),type="l",lty=2,col="pink")
legend(x="bottomright",legend=c("Male","Female"),col = c("black","red"),lty=1,bty="n",cex=0.75)

#logit piecewise CI for baseline hazard (and baseline survival?)

theta=m1.tstage$theta
psi=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=1)
cov_theta=m1.tstage$covar$M2QM2[-c(1:15),-c(1:15)]

h0=t(theta)%*%t(psi)

LL=UL=rep(0,length(time))

for(p in 1:length(time)){
  
  t=time[p]
  h0t=h0[,p]
  logith0t=log(h0t/(1-h0t))
  psi_t=psi[p,]
  
  dh0_dtheta=solve(h0t*(1-h0t))%*%t(psi_t)
  
  var_logit=dh0_dtheta%*%cov_theta%*%t(dh0_dtheta)
  se_logit=sqrt(var_logit)
  
  LL.logit=logith0t-1.96*se_logit
  UL.logit=logith0t+1.96*se_logit
  
  LL[p]=exp(LL.logit)/(1+exp(LL.logit))
  UL[p]=exp(UL.logit)/(1+exp(UL.logit))
  
}

plot(time,h0,type="l",ylim=c(0,1))
points(time,UL,type="l",lty=2)
points(time,LL,type="l",lty=2)


Psi=basis_phmc(time,knots=m1.tstage$knots,order=m1.tstage$control$order,which=2)
cov_theta=m1.tstage$covar$M2QM2[-c(1:15),-c(1:15)]

S0=exp(-t(theta)%*%t(Psi))

LL.s=UL.s=rep(0,length(time))

for(p in 1:length(time)){
  
  t=time[p]
  S0t=S0[,p]
  logitS0t=log(S0t/(1-S0t))
  Psi_t=Psi[p,]
  
  dS0_dtheta=solve((1.0001-S0t))%*%(-Psi_t)
  
  var_logit=dS0_dtheta%*%cov_theta%*%t(dS0_dtheta)
  se_logit=sqrt(var_logit)
  
  LL.logit=logitS0t-1.96*se_logit
  UL.logit=logitS0t+1.96*se_logit
  
  LL.s[p]=exp(LL.logit)/(1+exp(LL.logit))
  UL.s[p]=exp(UL.logit)/(1+exp(UL.logit))
  
}


plot(time,S0,type="l",ylim=c(0,1),xlab="Time (years)",ylab="Baseline survival function",
     main="Estimate of the baseline survival function")
points(time,UL.s,type="l",col="grey",lty=2)
points(time,LL.s,type="l",col="grey",lty=2)





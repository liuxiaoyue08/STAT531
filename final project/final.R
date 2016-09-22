 library(pomp)
# library(ggplot2)
blowfly <- read.csv("blowfly4.csv",skip=3, header=T)
# #Analyse t<=400
blowfly <- blowfly[blowfly$day<=400,]
colnames(blowfly)=c("time","pop")
# head(blowfly)
# #plot data
# plot(pop~time,data=blowfly,type='o')
# 
#building pomp model
require(pomp)
blowflies <- pomp(blowfly,times="time",t0=0)
plot(blowflies)

#Adding Deterministic Ricker model
blowfly_skeleton <- Csnippet("DN = r*N*exp(-N);")

blowfly_statenames <- c("N")
blowfly_paramnames <- c("r")

blowflies <- pomp(blowflies,skeleton=blowfly_skeleton,skeleton.type='map',
                 paramnames=blowfly_paramnames,statenames=blowfly_statenames)

traj <- trajectory(blowfies,params=c(N.0=80,r=12),as.data.frame=T)
ggplot(data=traj,aes(x=time,y=N))+geom_line()

#Adding in the process model simulator
stochStep <- Csnippet("
                      e = rnorm(0,sigma);
                      N = r*N*exp(-N+e);
                      ")
pomp(blowflies,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=14),
     paramnames=c("r","sigma"),statenames=c("N","e")) -> blowflies


#simulate the stochastic Ricker model.
sim <- simulate(blowflies,params=c(N.0=1,e.0=0,r=18,sigma=0.5),
                as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

#Adding in the measurement model and parameters
rmeas <- Csnippet("pop = rpois(phi*N);")
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")

pomp(blowflies,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> blowflies

coef(blowflies) <- c(N.0=1,e.0=0,r=18,sigma=0.5,phi=200)

sims <- simulate(blowflies,nsim=3,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+
  facet_wrap(~sim)


#Reformulating the Ricker model
skel <- Csnippet("DN = r*N*exp(-N/k);")
pomp(blowflies,skeleton=skel,skeleton.type='map',statenames=c("N"),paramnames=c("r","k"))->blowflies
traj <- trajectory(blowflies,params=c(N.0=1,r=12, k=5),as.data.frame=TRUE)
ggplot(data=traj,aes(x=time,y=N))+geom_line()

stochStep <- Csnippet("N=r*N*exp(-N/k);")
pomp(blowflies,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),paramnames=c("r","k"),statenames=c("N")) -> blowflies

sim <- simulate(blowflies,params=c(N.0=1,r=12,k=5),as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

rmeas <- Csnippet("pop = rnbinom(phi*N,si);")
dmeas <- Csnippet("lik = dnbinom(pop,m,mu=phi*N,si,give_log);")
pomp(blowflies,rmeasure=rmeas,statenames=c("N"),paramnames=c("phi","si"))->blowflies
coef(blowflies) <- c(N.0=1,r=18,phi=100,si=0.5,k=10)
sims<-simulate(blowflies,nsim=5,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+
  facet_wrap(~sim)






#stochastic Beverton-Holt model
skel <- Csnippet("DN = ((a*N)/(1+(b*N)));")
blowflies <- pomp(blowflies,skeleton=skel,skeleton.type='map',paramnames=c("a", "b"),statenames=c("N"))
traj <- trajectory(blowflies,params=c(N.0=1,a=12,b=5),as.data.frame=TRUE)
ggplot(data=traj,aes(x=time,y=N))+geom_line()

stochStep <- Csnippet("
  e = rlnorm(mu,sigma);
                      N = ((a*N)/(1+(b*N)))*e;
                      ")
pomp(blowflies,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),paramnames=c("a","b", "mu","sigma"),statenames=c("N","e")) -> blowflies
sim <- simulate(blowflies,params=c(N.0=1,e.0=0,a=12,b=5,sigma=0.5,mu=-0.5*(0.5)^2),as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

rmeas <- Csnippet("pop = rpois(phi*N);")
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")
pomp(blowflies,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> blowflies
coef(blowflies) <- c(N.0=1,e.0=0,a=12,b=5,sigma=0.5,mu=-0.5*(0.5)^2,phi=50)
sims <- simulate(blowflies,nsim=5,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+facet_wrap(~sim)


#Proposed model  stochastic delay-difference models


rmeas <- Csnippet("pop = rnbinom(N,1/sigma.y^2);")
dmeas <- Csnippet("lik = dnbinom(pop,m,mu=N,1/sigma.y^2,give_log);")
params=c("P","N0","delta","sigma.P","sigma.d","sigma.y")
states=c("N1","R","S","e","eps")





#stochastic Beverton-Holt model
library(pomp)
library(ggplot2)
blowfly <- read.csv("blowfly4.csv",skip=3, header=T)
#Analyse t<=400
blowfly <- blowfly[blowfly$day<=400,]
colnames(blowfly)=c("time","pop")
head(blowfly)
#plot data
plot(pop~time,data=blowfly,type='o')

#building pomp model
require(pomp)
blowflies <- pomp(blowfly,times="time",t0=0)
plot(blowflies)
skel <- Csnippet("DN = ((a*N)/(1+(b*N)));")
blowflies <- pomp(blowflies,skeleton=skel,skeleton.type='map',paramnames=c("a", "b"),statenames=c("N"))
traj <- trajectory(blowflies,params=c(N.0=1,a=12,b=5),as.data.frame=TRUE)
ggplot(data=traj,aes(x=time,y=N))+geom_line()

stochStep <- Csnippet("
  e = rlnorm(mu,sigma);
  N = ((a*N)/(1+(b*N)))*e;
                      ")
pomp(blowflies,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),paramnames=c("a","b", "mu","sigma"),statenames=c("N","e")) -> blowflies
sim <- simulate(blowflies,params=c(N.0=1,e.0=0,a=12,b=5,sigma=0.5,mu=-0.5*(0.5)^2),as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

rmeas <- Csnippet("pop = rpois(phi*N);")
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")
pomp(blowflies,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> blowflies
coef(blowflies) <- c(N.0=1,e.0=0,a=12,b=5,sigma=0.5,mu=-0.5*(0.5)^2,phi=50)
sims <- simulate(blowflies,nsim=5,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+facet_wrap(~sim)


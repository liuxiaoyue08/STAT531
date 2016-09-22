library(pomp)
library(ggplot2)
set.seed(123456789)
sheeps <- read.csv("annual-sheep-population-1000s-in.csv",skip=3, header=T)
colnames(sheeps)<-c("year","pop")
plot(pop~year,data=sheeps,type='o')

sheep <- pomp(sheeps,times="year",t0=1870)

skel <- Csnippet("DN = ((a*N)/(1+(b*N)));")
sheep <- pomp(sheep,skeleton=skel,skeleton.type='map',paramnames=c("a", "b"),statenames=c("N"))
traj <- trajectory(sheep,params=c(N.0=1,a=12,b=5),as.data.frame=TRUE)
ggplot(data=traj,aes(x=time,y=N))+geom_line()

stochStep <- Csnippet("
                      e = rlnorm(-0.5*sigma*sigma,sigma);
                      N = ((a*N)/(1+(b*N)))*e;
                      ")
pomp(sheep,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),paramnames=c("a","b","sigma"),statenames=c("N","e")) -> sheep
sim <- simulate(sheep,params=c(N.0=1,e.0=0,a=12,b=5,sigma=0.5),as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

rmeas <- Csnippet("pop = rpois(phi*N);")
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")
pomp(sheep,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> sheep
coef(sheep) <- c(N.0=2000,e.0=0,a=2.7,b=0.001,sigma=0.1,phi=1)
sims <- simulate(sheep,nsim=5,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+facet_wrap(~sim)

sims <- simulate(sheep,params=c(N.0=2150,e.0=0,a=2.8,b=0.001,sigma=0.05,phi=1),nsim=1,
                 as.data.frame=TRUE,include.data=TRUE)

ggplot(sims,mapping=aes(x=time,y=pop,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)


simulate(sheep,params=c(N.0=2000,e.0=0,a=2.7,b=0.001,sigma=0.1,phi=1),
         nsim=10000,states=TRUE) -> x

ell <- dmeasure(sheep,y=obs(sheep),x=x,times=time(sheep),log=TRUE,
                params=c(N.0=2000,e.0=0,a=2.7,b=0.001,sigma=0.1,phi=1))
dim(ell)
ell <- apply(ell,1,sum); summary(exp(ell)); logmeanexp(ell,se=TRUE)
pf <- pfilter(sheep,Np=5000,params=c(N.0=2000,e.0=0,a=2.7,b=0.001,sigma=0.1,phi=1))
logLik(pf)
pf <- replicate(10,pfilter(sheep,Np=5000,params=c(N.0=2000,e.0=0,a=2.7,b=0.001,sigma=0.1,phi=1)))
ll <- sapply(pf,logLik); ll
logmeanexp(ll,se=TRUE)


sliceDesign(
  c(N.0=2000,e.0=0,a=2.7,b=0.001,sigma=0.1,phi=1),
  a=rep(seq(from=0.5,to=4,length=40),each=3),
  b=rep(seq(from=0.001,to=1,length=40),each=3),
  sigma=rep(seq(from=0.05,to=0.5,length=40),each=3)) -> p

require(foreach)
require(doMC)

registerDoMC(cores=5)   
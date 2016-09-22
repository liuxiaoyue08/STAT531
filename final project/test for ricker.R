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

ell <- dmeasure(blowflies,y=obs(blowflies),x=sims,times=time(blowflies),log=TRUE,
                params=c(N.0=1,r=18,phi=100,si=0.5,k=10))
dim(ell)
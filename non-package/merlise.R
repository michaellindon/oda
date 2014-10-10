rm(list=ls())
library(oda)

#Generate Design Matrix
set.seed(1)
no=200
p=1000

xo=matrix( rnorm(no*p,mean=0,sd=1), no, p) 

#Scale and Center Design Matrix
xo=scale(xo,center=T,scale=F)
var=apply(xo^2,2,sum)
xo=scale(xo,center=F,scale=sqrt(var/no))

#Generate Data under True Model
p=dim(xo)[2]
b=rep(0,p)
b[1]=1
b[2]=2
b[3]=3
b[4]=4
b[5]=5
yo=xo%*%b+rnorm(no,0,1)
yo=yo-mean(yo)

#Scale Data and Produce xa
xo=scale(xo,center=T,scale=F)
var=apply(xo^2,2,sum)
xo=scale(xo,center=F,scale=sqrt(var/no))
xoxo=t(xo)%*%xo
A=-xoxo
foo=eigen(A)
L=foo$values
d=rep(-min(L),p)
L=L-min(L)
L=sqrt(L)
V=foo$vectors
xa=diag(L)%*%t(V)

lam=rep(1,p) #Ridge Parameter
priorprob=rep(0.01,p) #Bernoulli prior probs
burnin=1000
iterations=4000
dof=1

collapsed_normal=col_normal_gibbs(yo,xo,xa,lam,priorprob,burnin,iterations)
normal=normal_gibbs(yo,xo,xa,d,lam,priorprob,burnin,iterations)
plot(collapsed_normal$prob)
points(normal$prob,col="green")


collapsed_cauchy=col_t_gibbs(yo,xo,xa,priorprob,burnin,iterations,dof)
cauchy=t_gibbs(yo,xo,xa,d,priorprob,burnin,iterations,dof)
plot(collapsed_cauchy$prob)
points(cauchy$prob,col="green")



#EM fun
col_normal_em=col_normal_em(yo,xo,xa,d,lam,priorprob,3) ###
normal_em=normal_em(yo,xo,xa,d,lam,priorprob) ###no greedy yet
cauchy_em=t_em(yo,xo,xa,d,priorprob,dof,0) ###no greedy yet



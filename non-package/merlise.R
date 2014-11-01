library(oda)

#Generate Design Matrix
set.seed(1)
no=10000
p=15

xo=matrix(runif(no*p), no, p) 

#Scale and Center Design Matrix
xo=scale(xo,center=T,scale=F)
var=apply(xo^2,2,sum)
xo=scale(xo,center=F,scale=sqrt(var/no))

#Generate Data under True Model
p=dim(xo)[2]
b=rep(0,p)
b[1:4] = 4:1
yo=xo%*%b+rnorm(no,0,8.3) # + 5  for Steve but mean substracted
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
priorprob=rep(1/p,p) #Bernoulli prior probs
burnin=1000
iterations=10000
dof=1

#collapsed_normal=col_normal_gibbs(yo,xo,xa,lam,priorprob,burnin,iterations)
normal=normal_gibbs(yo,xo,xa,d,lam,priorprob,burnin,iterations)
#plot(collapsed_normal$prob)
                                        #points(normal$prob,col="green")
plot(normal$prob,col="green")


#collapsed_cauchy=col_t_gibbs(yo,xo,xa,priorprob,burnin,iterations,dof)
#cauchy=t_gibbs(yo,xo,xa,d,priorprob,burnin,iterations,dof)
#plot(collapsed_cauchy$prob)
#points(cauchy$prob,col="green")

library(BAS)
temp = data.frame(Y=yo, X=xo)

bas.lm(Y ~ ., prior="ZS-null", data=temp, modelprior=beta.binomial(1/p))

#EM fun
col_normal_em=col_normal_em(yo,xo,xa,d,lam,priorprob,3) ###
normal_em=normal_em(yo,xo,xa,d,lam,priorprob) ###no greedy yet
cauchy_em=t_em(yo,xo,xa,d,priorprob,dof,0) ###no greedy yet

summary(lm(Y ~ ., data=temp))

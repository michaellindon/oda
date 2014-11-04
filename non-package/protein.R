library(oda)
protein <- read.table("protein.data", header=T)

protein.lm <- lm(prot.act4 ~ 
	buf + buf*Ph + buf*NaCl + buf*con + buf*ra +
	 buf*det + buf*MgCl2 + buf*temp + 
	Ph + Ph*NaCl + Ph*con + Ph*ra + Ph*det +
	 Ph*MgCl2 +Ph*temp +
   	NaCl + NaCl*con + NaCl*ra + NaCl*det +
	 NaCl*MgCl2 + NaCl*temp +
	con + con*ra + con*det + 
	 con*MgCl2 +con*temp +
	ra + ra*det + ra*MgCl2 + ra*temp +
	det + det*MgCl2 + det*temp +
	MgCl2 + MgCl2*temp + I(NaCl^2) + I(Ph^2) + I(con^2) + I(temp^2), 
	data=protein,x=T)

summary(protein.lm)
xo=protein.lm$x
xo=xo[,-1]
yo=protein$prot.act4

#Scale Data and Produce xa
no=length(yo)
p=dim(xo)[2]
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
xaxa=t(xa)%*%xa

lam=rep(1,p)
#Fixed prior inclusion probabilities
priorprob=rep(0.1,p)
col_normal=col_normal_em(yo,xo,xaxa,d,lam,priorprob,3)
col_n_gibbs=col_normal_gibbs(yo,xo,xa,lam,priorprob,500,20000)


#Beta-Binomial Priors Beta(eta,zeta)
eta=1
zeta=1
col_normal_bb=col_normal_em_bb(yo,xo,xaxa,d,lam,eta,zeta,3)
col_n_gibbs_bb=col_normal_gibbs_bb(yo,xo,xa,lam,eta,zeta,500,20000)







col_normal=col_normal_em(yo,xo,xaxa,do,lam,priorprob,3)
col_n_gibbs=col_normal_gibbs(yo,xo,xa,lam,priorprob,500,20000)
plot(col_n_gibbs$prob)
abline(v=col_normal$top_model)
col_normalbb=col_normal_em_bb(yo,xo,xaxa,do,rep(1,p),1,p,3)
col_n_gibbs=col_normal_gibbs(yo,xo,xa,lam,priorprob,500,20000)

#Scale Data and Produce xa
no=length(yo)
p=dim(xo)[2]
xo=scale(xo,center=T,scale=F)
var=apply(xo^2,2,sum)
xo=scale(xo,center=F,scale=sqrt(var/no))
xoxo=t(xo)%*%xo
A=-xoxo
diag(A)=0
diag(A)=abs(min(eigen(A)$values))+0.001
d=diag(A)+diag(xoxo)
xa=chol(A)
V=eigen(A)$vectors
L=eigen(A)$values
L=sqrt(L)
xa=diag(L)%*%t(V)


p=dim(xo)[2]
lam=rep(1,p)
priorprob=rep(0.5,p)
m_gibbs=mixture_gibbs(yo,xo,xa,lam,priorprob,500,200000,1)
n_gibbs=normal_gibbs(yo,xo,xa,d,lam,priorprob,500,20000)
col_n_gibbs=col_normal_gibbs(yo,xo,xa,lam,priorprob,500,20000)
plot(n_gibbs$prob)

n_var=normal_var(yo,xo,xa,d,lam,priorprob,0.001,100)
# n_var$lb_trace[2:dim(n_var$lb_trace)[1],1]-n_var$lb_trace[1:(dim(n_var$lb_trace)[1]-1),1]>=0
mean(n_var$lb_trace[2:dim(n_var$lb_trace)[1],1]-n_var$lb_trace[1:(dim(n_var$lb_trace)[1]-1),1]>=0)
plot(n_var$prob_trace[1,],ylim=c(0,1),type="l")
for(i in 2:p) lines(n_var$prob_trace[i,],ylim=c(0,1),col=i)
plot(n_var$lb_trace)

m_var=mixture_var(yo,xo,xa,d,lam,priorprob,1,0.001,100)
m_var$lb_trace[2:dim(m_var$lb_trace)[1],1]-m_var$lb_trace[1:(dim(m_var$lb_trace)[1]-1),1]>=0
mean(m_var$lb_trace[2:dim(m_var$lb_trace)[1],1]-m_var$lb_trace[1:(dim(m_var$lb_trace)[1]-1),1]>=0)
plot(m_var$prob_trace[1,],ylim=c(0,1),type="l")
for(i in 2:p) lines(m_var$prob_trace[i,],ylim=c(0,1),col=i)
plot(m_var$lb_trace)

em=normal_em(yo,xo,xa,d,rep(1,p),rep(0.3,p))
em$lpd_trace[2:dim(em$lpd_trace)[1],1]-em$lpd_trace[1:(dim(em$lpd_trace)[1]-1),1]>=0
mean(em$lpd_trace[2:dim(em$lpd_trace)[1],1]-em$lpd_trace[1:(dim(em$lpd_trace)[1]-1),1]>=0)
plot(em$B_trace[1,],ylim=c(0,1),type="l")
for(i in 2:p) lines(em$B_trace[i,],ylim=c(0,1),col=i)
plot(em$lpd_trace)

lasso=lasso_em(yo,xo,xa,d,1000)
lasso$lpd_trace[2:dim(lasso$lpd_trace)[1],1]-lasso$lpd_trace[1:(dim(lasso$lpd_trace)[1]-1),1]>=0
mean(lasso$lpd_trace[2:dim(lasso$lpd_trace)[1],1]-lasso$lpd_trace[1:(dim(lasso$lpd_trace)[1]-1),1]>=0)
plot(lasso$B)
which(round(lasso$B,2)!=0)

stephens <- varbvsoptimize(xo,yo,4,1/4,log(4.5),alpha0=NULL,mu0=NULL,verbose=FALSE)
which(stephens$alpha>0.5)
bols=solve(t(xo)%*%xo)%*%t(xo)%*%yo
plot(n_gibbs$prob)
abline(v=sort(abs(bols),index.return=TRUE,decreasing=TRUE)$ix[1:10],col="orange")
bridge=solve(t(xo)%*%xo+2*no*diag(p))%*%t(xo)%*%yo
plot(n_gibbs$prob)
abline(v=sort(abs(bridge),index.return=TRUE,decreasing=TRUE)$ix[1:10],col="red")
em=normal_em(yo,xo,xa,d,rep(1,p),rep(0.5,p))
plot(n_gibbs$prob)
abline(v=which(em$prob>0.5),col="blue")
abline(v=em$top_model,col="red")

library(BAS, lib.loc="~/Rlib64")

protein.bma = bas.lm(prot.act4 ~ 
	buf + buf*Ph + buf*NaCl + buf*con + buf*ra +
	 buf*det + buf*MgCl2 + buf*temp + 
	Ph + Ph*NaCl + Ph*con + Ph*ra + Ph*det +
	 Ph*MgCl2 +Ph*temp +
   	NaCl + NaCl*con + NaCl*ra + NaCl*det +
	 NaCl*MgCl2 + NaCl*temp +
	con + con*ra + con*det + 
	 con*MgCl2 +con*temp +
	ra + ra*det + ra*MgCl2 + ra*temp +
	det + det*MgCl2 + det*temp +
	MgCl2 + MgCl2*temp + I(NaCl^2) + I(Ph^2) + I(con^2) + I(temp^2),
 data=protein, n.models=2^20, prior="ZS-null",initprobs="Uniform", update=5000)

plot(protein.bma$logmarg)
plot(cumsum(exp(protein.bma$logmarg - max(protein.bma$logmarg))))
plot(protein.bma)
image(protein.bma,top=50)
summary(protein.bma)

#> max(protein.bma$logmarg)
#[1] 44.75594

#> max(protein.bma$logmarg)
#[1] 41.75033

Xnames = protein.bma$namesx

protein.step = lm(formula = prot.act4 ~ buf + Ph + NaCl + con + ra + det + MgCl2 +     temp + I(NaCl^2) + I(con^2) + buf:Ph + buf:NaCl + buf:con +     buf:ra + buf:det + buf:MgCl2 + buf:temp + Ph:con + Ph:ra +     Ph:det + Ph:MgCl2 + NaCl:con + NaCl:ra + NaCl:det + NaCl:MgCl2 +     con:ra + con:det + ra:det + ra:temp + det:MgCl2 + det:temp +     MgCl2:temp, data = protein, x = T)


protein.bma2 = bas.lm(prot.act4 ~ buf + Ph + NaCl + con + ra + det + MgCl2 +     temp + I(NaCl^2) + I(con^2) + buf:Ph + buf:NaCl + buf:con +     buf:ra + buf:det + buf:MgCl2 + buf:temp + Ph:con + Ph:ra +     Ph:det + Ph:MgCl2 + NaCl:con + NaCl:ra + NaCl:det + NaCl:MgCl2 +     con:ra + con:det + ra:det + ra:temp + det:MgCl2 + det:temp +     MgCl2:temp, data = protein, initprobs="eplogp", prior="ZS-null", n.models=1^15, update=5000)


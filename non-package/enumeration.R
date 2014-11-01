lpostprob <- function(gammanew,xo,yo,lam,no,p,a,yoyo,priorodds)
{
  Lamg=diag(lam[gammanew==1]);
  ncolumn <- sum(gammanew == 1);
  if (ncolumn == 0)
  {
    b=0.5*yoyo;
    lpprob = -(0.5*(no-1))*log(2*pi)+(lgamma(a)-a*log(b)) ;
  } else
  {
    xog <- xo[,gammanew==1];
    b=0.5*(yoyo-t(yo)%*%xog%*%solve(t(xog)%*%xog+Lamg)%*%t(xog)%*%yo);
    lpprob =sum(gammanew*log(priorodds)) -(0.5*(no-1))*log(2*pi)+lgamma(a) -a*log(b)+0.5*log(det(as.matrix(Lamg))) -0.5*log(det(as.matrix(Lamg+t(xog)%*%xog)))   ;
  }
  return((lpprob));
}

xo <- as.matrix(xo); yo <- as.vector(yo);
yo=yo-mean(yo);
p=dim(xo)[2]
no=dim(xo)[1]
lam=rep(1,p)
a=(no-1)/2
yoyo=t(yo)%*%yo
priorodds=priorprob/(1-priorprob)

tf <- c(TRUE, FALSE)
models <- expand.grid(replicate(p,tf,simplify=FALSE))
names(models) <- NULL
models=as.matrix(models)
lpp.all=apply(models,1,lpostprob,xo,yo,lam,no,p,a,yoyo,priorodds)

results=cbind(lpp.all, models)
order=sort(results[,1],index=TRUE,decreasing=TRUE)
results[order$ix,][1:10,]
results[order$ix,1]=results[order$ix,1]-results[order$ix[1],1]
results[order$ix,1]=exp(results[order$ix,1])
nconstant=sum(results[,1])
results[,1]=results[,1]/nconstant
postprob=results[order$ix,]  
round(postprob[1:10,],3)
inclusionprob=rep(0,dim(postprob)[2]-1)
for(i in 1:dim(postprob)[2]-1){
  inclusionprob[i]=sum(postprob[,i+1]*postprob[,1]);
}
round(inclusionprob,4)

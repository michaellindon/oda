blm<-function( yo,xo, betaprior, lam, tdof, modelprior,  priorprob, beta1, beta2, burnin, niter,  collapsed){
  no=length(yo)
  p=dim(as.matrix(xo))[[2]]
  na=p
  
  if(modelprior=="betabinomial" || modelprior=="BetaBinomial" || modelprior=="BETABINOMIAL")
    {
    modelprior=2
    print("Beta Binomial")
  }else if(modelprior=="bernoulli" || modelprior=="Bernoulli" || modelprior=="BERNOULLI")
    {
    modelprior=1
    print("Bernoulli")
  }else{
    #uniform default
  modelprior=0
  print("Uniform")
  }
  
  if(betaprior=="normal" || betaprior=="Normal" || betaprior=="NORMAL" || betaprior=="gaussian" || betaprior=="Gaussian" || betaprior=="GAUSSIAN")
    {
    scalemixture=0
    print("Normal")
  }else if(betaprior=="t"){
    scalemixture=1
    print("t")
  }else{
    #cauchy default
    scalemixture =1 
    tdof=1
    print("Cauchy")
  }
  if(collapsed==0) print("Uncollapsed")
  if(collapsed==1) print("Collapsed")
  
  coefficients=rep(0,p)
  pip=rep(0,p)
  B_mcmc=rep(0,p*niter)
  lam_mcmc=rep(0,p*niter)
  prob_mcmc=rep(0,p*niter)
  gamma_mcmc=rep(0,p*niter)
  phi_mcmc=rep(0,niter)
  intercept_mcmc=rep(0,niter)
  xo_scale=rep(0,p)
  gibbs=.C("lm_gibbs", as.double(yo), as.double(xo), as.double(lam), as.integer(modelprior),as.double(priorprob), as.double(beta1), as.double(beta2), as.integer(burnin), as.integer(niter), as.integer(scalemixture), as.double(tdof), as.integer(collapsed), as.integer(no), as.integer( na), as.integer(p),as.double(B_mcmc),as.double(prob_mcmc),as.integer(gamma_mcmc),as.double(phi_mcmc),as.double(lam_mcmc),as.double(coefficients),as.double(pip),as.double(intercept_mcmc),as.double(xo_scale),package="oda")
return(list(xo_scale=gibbs[[24]],model_mcmc=matrix(gibbs[[18]],p,niter),model=cbind(1,matrix(gibbs[[2]],no,p)),intercept_mcmc=gibbs[[23]],coefficients=c(mean(gibbs[[23]]),gibbs[[21]]),pip=gibbs[[22]],lam_mcmc=matrix(gibbs[[20]],p,niter),phi_mcmc=gibbs[[19]],B_mcmc=matrix(gibbs[[16]],p,niter),prob_mcmc=matrix(gibbs[[17]],p,niter),niter=niter,no=no,p=p))
}




bglm<-function( yo,xo, betaprior, lam, tdof, modelprior,  priorprob, beta1, beta2, burnin, niter){
  no=length(yo)
  p=dim(as.matrix(xo))[[2]]
  na=p
  
  if(modelprior=="betabinomial" || modelprior=="BetaBinomial" || modelprior=="BETABINOMIAL")
  {
    modelprior=2
    print("Beta Binomial")
  }else if(modelprior=="bernoulli" || modelprior=="Bernoulli" || modelprior=="BERNOULLI")
  {
    modelprior=1
    print("Bernoulli")
  }else{
    #uniform default
    modelprior=0
    print("Uniform")
  }
  
  if(betaprior=="normal" || betaprior=="Normal" || betaprior=="NORMAL" || betaprior=="gaussian" || betaprior=="Gaussian" || betaprior=="GAUSSIAN")
  {
    scalemixture=0
    print("Normal")
  }else if(betaprior=="t"){
    scalemixture=1
    print("t")
  }else{
    #cauchy default
    scalemixture =1 
    tdof=1
    print("Cauchy")
  }

  
  coefficients=rep(0,p)
  pip=rep(0,p)
  B_mcmc=rep(0,p*niter)
  lam_mcmc=rep(0,p*niter)
  prob_mcmc=rep(0,p*niter)
  gamma_mcmc=rep(0,p*niter)
  phi_mcmc=rep(0,niter)
  intercept_mcmc=rep(0,niter)
  xo_scale=rep(0,p)
  gibbs=.C("glm_gibbs", as.double(yo), as.double(xo), as.double(lam), as.integer(modelprior),as.double(priorprob), as.double(beta1), as.double(beta2), as.integer(burnin), as.integer(niter), as.integer(scalemixture), as.double(tdof), as.integer(no), as.integer( na), as.integer(p),as.double(B_mcmc),as.double(prob_mcmc),as.integer(gamma_mcmc),as.double(phi_mcmc),as.double(lam_mcmc),as.double(coefficients),as.double(pip),as.double(intercept_mcmc),as.double(xo_scale),package="oda")
  return(list(xo_scale=gibbs[[23]],coefficients=c(mean(gibbs[[22]]),gibbs[[20]]),pip=gibbs[[21]],lam_mcmc=matrix(gibbs[[19]],p,niter),phi_mcmc=matrix(gibbs[[18]],p,niter),B_mcmc=matrix(gibbs[[15]],p,niter),gamma_mcmc=matrix(gibbs[[17]],p,niter),prob_mcmc=matrix(gibbs[[16]],p,niter)))
}

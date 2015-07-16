#include "oda.h"


double log_posterior_density(std::vector<int>& gamma, std::vector<double>& priorprob, double phi, double a, double b);

extern "C" void em(double * ryo, double * rxo,  double * rlam, int * rmodelprior, double * rpriorprob, double * rbeta1, double * rbeta2, int * rmaxniter, int * rscalemixture, double * ralpha,  int * rno, int * rna, int * rp, double * B_trace, double * prob_trace, int * gamma_trace, double * phi_trace, double * lam_trace, double * lpd_trace, double * xo_scale, double * rtol)
{
	GetRNGstate();
	//MCMC Variables//
	double tol=*rtol;
	int niter=*rmaxniter;

	//Dimensions//
	int p=*rp;
	int no=*rno;
	int na=*rna;

	//Phi Variables//
	double a,b;
	double phi=1.0;

	//Yo Variables//
	std::vector<double> yo(ryo, ryo+no); 
	double yobar=center_yo(yo);
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);
	std::vector<double> xo(rxo, rxo+no*p);
	standardize_xo(xo,xo_scale,no,p);
	std::copy(xo.begin(),xo.end(),rxo);

	std::vector<double> xoyo(p);
	dgemv_( &transT, &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);

	std::vector<double> xoxo(p*p);
	dgemm_( &transT, &transN, &p, &p, &no, &unity, &*xo.begin(), &no, &*xo.begin(), &no, &inputscale0, &*xoxo.begin(), &p );

	//Construct Xa//
	std::vector<double> xa(p*(p+1)/2); //Triangular Packed Storage
	std::vector<double> d(p);
	chol_xa(xa,xoxo,d,p);


	//Reserve Memory for Submatrices//
	std::vector<double> xog; xog.reserve(no*p);
	std::vector<double> xogyo; xogyo.reserve(p);
	std::vector<double> xogxog_Lamg; xogxog_Lamg.reserve(p*p);
	std::vector<double> xag; xag.reserve(na*p);

	//Ya Variables//
	std::vector<double> xaya(p);

	//Beta Variables//
	std::vector<double> Bols(p);
	std::vector<double> B(p,0.0);
	std::vector<double> Bg; Bg.reserve(p);
	for(int i=0; i<p; ++i){
		B[i]=Rf_rnorm(0,1);
	}


	//Lambda Variables//
	int scalemixture=*rscalemixture;
	double alpha=*ralpha;
	std::vector<double> lam(rlam,rlam+p);
	std::vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Gamma Variables//
	std::vector<int> gamma(p,0);
	for(int i=0; i<p; ++i){
		if(Rf_runif(0,1)<0.5){
			gamma[i]=1; //Note gamma is initialized at zero
		}
	}
	int p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);
	bool gamma_diff=true;
	int modelprior=*rmodelprior;

	//Probability Variables//
	std::vector<double> prob(p);
	std::vector<double> odds(p);
	std::vector<double> priorprob(rpriorprob,rpriorprob+p);

	//Theta Variables//
	double theta=0.5;
	double beta1=*rbeta1;
	double beta2=*rbeta2;

	//Store Initial Values//
	std::copy(B.begin(),B.end(),B_trace);
	std::copy(prob.begin(),prob.end(),prob_trace);
	std::copy(gamma.begin(),gamma.end(),gamma_trace);
	std::copy(lam.begin(),lam.end(),lam_trace);
	phi_trace[0]=phi;

	//Run EM Algorithm//
	int t=1;
	do
	{
		if(p_gamma) submatrices_uncollapsed(gamma_diff,B,xog,xag,lamg,Bg,gamma,lam,xo,xa,p_gamma,p,no,na);
		if(p_gamma){
			a=0.5*(no+p_gamma-3);
			std::vector<double> residual=yo;
			dgemv_(&transN , &no, &p_gamma, &nunity, &*xog.begin(), &no, &*Bg.begin(), &inc, &inputscale1, &*residual.begin(), &inc);
			b=0.5*ddot_(&no, &*residual.begin(), &inc, &*residual.begin(), &inc);
			for(size_t i=0; i!=Bg.size(); ++i) b+=0.5*Bg[i]*Bg[i]*lamg[i]; 
		}else{
			a=0.5*(no-3);
			b=0.5*yoyo;
		}
		lpd_trace[t-1]=log_posterior_density(gamma,priorprob,phi,a,b);
		phi=a/b;

		if(p_gamma!=0){
			dgemv_(&transN , &na, &p_gamma, &unity, &*xag.begin(), &na, &*Bg.begin(), &inc, &inputscale0, &*xaya.begin(), &inc);
			//dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
			dtpmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &*xaya.begin(), &inc);
		}else{
			for(size_t i=0; i!=xaya.size(); ++i) xaya[i]=0;
			//dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
		}

		//Compute Probabilities//
		if(modelprior==1)
		{
			bernoulli_probabilities(prob,odds,Bols,d,xoyo,xaya,priorprob,lam,phi);
		}else if(modelprior==2)
		{
			betabinomial_probabilities(prob,odds,Bols,d,xoyo,xaya,theta,lam,phi);
		}else
		{
			uniform_probabilities(prob,odds,Bols,d,xoyo,xaya,lam,phi);
		}

		for (int i = 0; i < p; ++i)
		{
			if(prob[i]*sqrt(phi*(lam[i]+d[i])/(2*M_PI))>(1-prob[i])){
				gamma[i]=1;
			}else{
				gamma[i]=0;
			}
		}
		p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);


		for (int i = 0; i < p; ++i)
		{
			if(gamma[i]==1){
				B[i]=(d[i]/(d[i]+lam[i]))*Bols[i];
			}else{
				B[i]=0;
			}
		}
		//
		//Draw Theta//
//		if(modelprior==2) theta=Rf_rbeta(beta1+p_gamma,p-p_gamma+beta2);

		//Draw Lambda//
//		if(scalemixture) draw_lambda_t(lam,gamma,alpha,B,phi);


		//Store Draws//
		std::copy(gamma.begin(),gamma.end(),(gamma_trace+p*t));
		std::copy(prob.begin(),prob.end(),(prob_trace+p*t));
		std::copy(B.begin(),B.end(),(B_trace+p*t));
		std::copy(lam.begin(),lam.end(),(lam_trace+p*t));
		phi_trace[t]=phi;

		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_trace,t,p);
		t++;
	}while (t<niter && ((lpd_trace[t-1]-lpd_trace[t-2])*(lpd_trace[t-1]-lpd_trace[t-2]))>tol);
	PutRNGstate();
}


double log_posterior_density(std::vector<int>& gamma, std::vector<double>& priorprob, double phi, double a, double b){
	int sum=0;
	for (size_t i = 0; i < gamma.size(); ++i){
		sum+=gamma[i]*(log(priorprob[i])-log(1-priorprob[i]));
	}

	return(a*log(phi/(2*M_PI))-phi*b+sum);
}

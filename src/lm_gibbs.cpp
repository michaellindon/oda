#include "oda.h"
#include <iostream>

extern "C" void lm_gibbs(double * ryo, double * rxo,  double * rlam, int * rmodelprior, double * rpriorprob, double * rbeta1, double * rbeta2, int * rburnin, int * rniter, int * rscalemixture, double * ralpha, int * rcollapsed, int * rno, int * rna, int * rp, double * B_mcmc, double * prob_mcmc, int * gamma_mcmc, double * phi_mcmc, double * B_rb, double * prob_rb)
{
	//MCMC Variables//
	int burnin=*rburnin;
	int niter=*rniter;
	int collapsed=*rcollapsed;

	//Dimensions//
	int p=*rp;
	int no=*rno;
	int na=*rna;

	//Phi Variables//
	double b=1.0;
	double phi=1.0;

	//Yo Variables//
	std::vector<double> yo(ryo, ryo+no); 
	std::vector<double> xo(rxo, rxo+no*p);
	double yobar=scale(yo,xo,no,p);
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);

	std::vector<double> xoyo(p);
	dgemv_( &transT, &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);

	std::vector<double> xoxo(p*p);
	dgemm_( &transT, &transN, &p, &p, &no, &unity, &*xo.begin(), &no, &*xo.begin(), &no, &inputscale0, &*xoxo.begin(), &p );

	//Construct Xa//
	std::vector<double> xa(xoxo);
	std::vector<double> d(p);
	chol_xa(xa,xoxo,d,p);
	

	//Reserve Memory for Submatrices//
	std::vector<double> xog; xog.reserve(no*p);
	std::vector<double> xogyo; xogyo.reserve(p);
	std::vector<double> xogxog_Lamg; xogxog_Lamg.reserve(p*p);
	std::vector<double> xag; xag.reserve(na*p);

	//Ya Variables//
	std::vector<double> mu(na);
	std::vector<double> Z; Z.reserve(na);
	std::vector<double> xaya(p);

	//Beta Variables//
	std::vector<double> Bols(p);
	std::vector<double> B(p,0.0);
	std::vector<double> Bg; Bg.reserve(p);

	//Lambda Variables//
	int scalemixture=*rscalemixture;
	double alpha=*ralpha;
	std::vector<double> lam(rlam,rlam+p);
	std::vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Gamma Variables//
	std::vector<int> gamma(p,0);
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
	std::copy(B.begin(),B.end(),B_mcmc);
	std::copy(prob.begin(),prob.end(),prob_mcmc);
	std::copy(gamma.begin(),gamma.end(),gamma_mcmc);
	phi_mcmc[0]=phi;

	//Run Gibbs Sampler//
	for (int t = 1; t < niter; ++t)
	{
		switch(collapsed)
		{
			case false:
				if(p_gamma) submatrices_uncollapsed(gamma_diff,B,xog,xag,lamg,Bg,gamma,lam,xo,xa,p_gamma,b,p,no,na);
				phi=draw_uncollapsed_phi(p_gamma,no,yo,xog,Bg,lamg,yoyo);
				draw_uncollapsed_xaya(xaya,xa,xag,Bg,phi,Z,na,p,p_gamma);
				break;

			case true:
				switch(scalemixture)
				{
					case false:
						if(gamma_diff && p_gamma) submatrices_collapsed(gamma_diff,mu,xag,xogxog_Lamg,xogyo,lamg,Bg,gamma,xoyo,lam,xa,xoxo,p_gamma,b,yoyo,p,na);
						break;
					case true: 
						if(p_gamma) submatrices_collapsed(gamma_diff,mu,xag,xogxog_Lamg,xogyo,lamg,Bg,gamma,xoyo,lam,xa,xoxo,p_gamma,b,yoyo,p,na);
						break;
				}
				phi=draw_collapsed_phi(b,p_gamma,no, yoyo,xogyo,Bg);
				draw_collapsed_xaya(xaya,xa,xag,mu,phi,Z,xogxog_Lamg,na,p,p_gamma);

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

		//Draw Gamma//
		draw_gamma(gamma,p_gamma,prob);


		//Draw Theta//
		if(modelprior==2) theta=Rf_rbeta(beta1+p_gamma,p-p_gamma+beta2);


		//Draw Beta//
		draw_beta(gamma,B,Bols,d,lam,phi);


		//Draw Lambda//
		if(scalemixture) draw_lambda_t(lam,gamma,alpha,B,phi);


		//Store Draws//
		std::copy(gamma.begin(),gamma.end(),(gamma_mcmc+p*t));
		std::copy(prob.begin(),prob.end(),(prob_mcmc+p*t));
		std::copy(B.begin(),B.end(),(B_mcmc+p*t));
		phi_mcmc[t]=phi;

		//Rao Blackwell//
		if(t>=burnin) rao_blackwell(B_rb,prob_rb,B,prob,burnin,niter);

		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_mcmc,t,p);

	}
}

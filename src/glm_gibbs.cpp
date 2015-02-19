#include "oda.h"

extern "C" void glm_gibbs(double * rZ, double * rxo,  double * rlam, int * rmodelprior, double * rpriorprob, double * rbeta1, double * rbeta2, int * rburnin, int * rniter, int * rscalemixture, double * ralpha,  int * rno, int * rna, int * rp, double * B_mcmc, double * prob_mcmc, int * gamma_mcmc, double * phi_mcmc, double * lam_mcmc, double * B_rb, double * prob_rb)
{
	//MCMC Variables//
	int burnin=*rburnin;
	int niter=*rniter;

	//Dimensions//
	int p=*rp;
	int no=*rno;
	int na=*rna;

	//Phi Variables//
	double phi=1.0;

	//Yo Variables//
	std::vector<double> Z(rZ, rZ+no); 
	std::vector<double> xo(rxo, rxo+no*p);
	scale_xo(xo,no,p);
	std::vector<double> xoyo(p);
	double yobar=0;

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

	//Lambda Variables//
	int scalemixture=*rscalemixture;
	double alpha=*ralpha;
	std::vector<double> lam(rlam,rlam+p);
	std::vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Gamma Variables//
	std::vector<int> gamma(p,1);
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
	std::copy(lam.begin(),lam.end(),lam_mcmc);

	//Run Gibbs Sampler//
	for (int t = 1; t < niter; ++t)
	{

		//Form Submatrices//
		if(p_gamma) submatrices_uncollapsed(gamma_diff,B,xog,xag,lamg,Bg,gamma,lam,xo,xa,p_gamma,p,no,na);

		//Draw xoyo//
		draw_xoyo(Z,xoyo,yobar,xo,xog,Bg,phi,no,p,p_gamma);

		//Draw xaya//
		draw_uncollapsed_xaya(xaya,xa,xag,Bg,phi,na,p,p_gamma);

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
		std::copy(lam.begin(),lam.end(),(lam_mcmc+p*t));

		//Rao Blackwell//
		if(t>=burnin) rao_blackwell(B_rb,prob_rb,B,prob,burnin,niter);

		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_mcmc,t,p);

	}
}

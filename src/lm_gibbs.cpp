#include <Rcpp.h>
#include "oda.h"

// [[Rcpp::export]]
Rcpp::List lm_gibbs(Rcpp::NumericVector ryo, Rcpp::NumericMatrix rxo, Rcpp::NumericMatrix rxa, Rcpp::NumericMatrix rxoxo, Rcpp::NumericVector rd, Rcpp::NumericVector rlam, SEXP rmodelprior, Rcpp::NumericVector rpriorprob, SEXP rbeta1, SEXP rbeta2, SEXP rburnin, SEXP rniter, SEXP rscalemixture, SEXP ralpha, SEXP rcollapsed)
{

	//MCMC Variables//
	int burnin=Rcpp::as<int >(rburnin);
	int niter=Rcpp::as<int >(rniter);
	bool collapsed=Rcpp::as<bool >(rcollapsed);

	//Dimensions//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();

	//Create Model Matrices//
	std::vector<double> d(rd.begin(),rd.end());
	std::vector<double> xo(rxo.begin(), rxo.begin()+no*p);
	std::vector<double> xa(rxa.begin(), rxa.begin()+na*p);
	std::vector<double> xog; xog.reserve(no*p);
	std::vector<double> xag; xag.reserve(na*p);
	std::vector<double> xoxo(rxoxo.begin(),rxoxo.end());
	std::vector<double> xogyo; xogyo.reserve(p);
	std::vector<double> xogxog_Lamg; xogxog_Lamg.reserve(p*p);

	//Phi Variables//
	double b=1.0;
	double phi=1.0;

	//Yo Variables//
	std::vector<double> yo(ryo.begin(), ryo.end()); 
	double yobar=std::accumulate(yo.begin(),yo.end(),0.0)/yo.size();
	for(std::vector<double>::iterator it=yo.begin(); it!=yo.end(); ++it) *it-=yobar;
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);
	std::vector<double> xoyo(p);
	dgemv_( &transT, &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);

	//Ya Variables//
	std::vector<double> mu(na);
	std::vector<double> Z; Z.reserve(na);
	std::vector<double> xaya(p);

	//Beta Variables//
	std::vector<double> Bols(p);
	std::vector<double> B(p,0.0);
	std::vector<double> Bg; Bg.reserve(p);

	//Lambda Variables//
	int scalemixture=Rcpp::as<bool >(rscalemixture);
	double alpha=Rcpp::as<double >(ralpha);
	std::vector<double> lam(rlam.begin(),rlam.end());
	std::vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Gamma Variables//
	std::vector<int> gamma(p,0);
	int p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);
	bool gamma_diff=true;
	std::string modelprior=Rcpp::as<std::string >(rmodelprior);

	//Probability Variables//
	std::vector<double> prob(p);
	std::vector<double> odds(p);
	std::vector<double> priorprob(rpriorprob.begin(),rpriorprob.end());

	//Theta Variables//
	double theta=0.5;
	double beta1=Rcpp::as<double >(rbeta1);
	double beta2=Rcpp::as<double >(rbeta2);

	//Allocate Space for MCMC Draws//
	std::vector<double> B_mcmc(p*niter);
	std::vector<double> prob_mcmc(p*niter);
	std::vector<int>  gamma_mcmc(p*niter);
	std::vector<double> phi_mcmc(niter);

	//Store Initial Values//
	std::copy(B.begin(),B.end(),B_mcmc.begin());
	std::copy(prob.begin(),prob.end(),prob_mcmc.begin());
	std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin());
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
		if(modelprior=="bernoulli")
		{
			bernoulli_probabilities(prob,odds,Bols,d,xoyo,xaya,priorprob,lam,phi);
		}else if(modelprior=="betabinomial")
		{
			betabinomial_probabilities(prob,odds,Bols,d,xoyo,xaya,theta,lam,phi);
		}else
		{
			uniform_probabilities(prob,odds,Bols,d,xoyo,xaya,lam,phi);
		}

		//Draw Gamma//
		draw_gamma(gamma,p_gamma,prob);


		//Draw Theta//
		if(modelprior=="betabinomial") theta=Rf_rbeta(beta1+p_gamma,p-p_gamma+beta2);


		//Draw Beta//
		draw_beta(gamma,B,Bols,d,lam,phi);


		//Draw Lambda//
		if(scalemixture) draw_lambda_t(lam,gamma,alpha,B,phi);


		//Store Draws//
		std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin()+p*t);
		std::copy(prob.begin(),prob.end(),prob_mcmc.begin()+p*t);
		std::copy(B.begin(),B.end(),B_mcmc.begin()+p*t);
		phi_mcmc[t]=phi;


		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_mcmc,t,p);

	}


	//Return Data to R//
	return Rcpp::List::create
		(
		 Rcpp::Named("B_mcmc") = B_mcmc,
		 Rcpp::Named("phi_mcmc") = phi_mcmc,
		 Rcpp::Named("prob_mcmc") = prob_mcmc,
		 Rcpp::Named("gamma_mcmc") = gamma_mcmc
		) ;
}

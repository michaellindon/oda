#include <Rcpp.h>
#include "oda.h"

// [[Rcpp::export]]
Rcpp::List oda_gibbs(Rcpp::NumericVector ryo, Rcpp::NumericMatrix rxo, Rcpp::NumericMatrix rxa, Rcpp::NumericMatrix rxoxo, Rcpp::NumericVector rd, Rcpp::NumericVector rlam, SEXP rmodelprior, Rcpp::NumericVector rpriorprob, SEXP rbeta1, SEXP rbeta2, SEXP rburnin, SEXP rniter, SEXP rscalemixture, SEXP ralpha)
{

	//Dimensions//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();

	//Construct Vectors//
	std::vector<double> xo(rxo.begin(), rxo.begin()+no*p);
	std::vector<double> xa(rxa.begin(), rxa.begin()+na*p);
	std::vector<double> yo(ryo.begin(), ryo.end()); 
	double yobar=std::accumulate(yo.begin(),yo.end(),0.0)/yo.size();
	for(std::vector<double>::iterator it=yo.begin(); it!=yo.end(); ++it) *it-=yobar;
	std::vector<double> lam(rlam.begin(),rlam.end());
	std::vector<double> d(rd.begin(),rd.end());
	std::vector<double> xoxo(rxoxo.begin(),rxoxo.end());
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);
	int scalemixture=Rcpp::as<bool >(rscalemixture);
	double alpha=Rcpp::as<double >(ralpha);
	std::string modelprior=Rcpp::as<std::string >(rmodelprior);
	double beta1=Rcpp::as<double >(rbeta1);
	double beta2=Rcpp::as<double >(rbeta2);
	std::vector<double> priorprob(rpriorprob.begin(),rpriorprob.end());

	//Create Matrices//
	std::vector<double> xoyo(p);
	dgemv_( &transT, &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);
	std::vector<double> xogyo; xogyo.reserve(p);
	std::vector<double> xogxog_Lamg; xogxog_Lamg.reserve(p*p);
	std::vector<double> xag; xag.reserve(na*p);
	std::vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Phi Variables//
	double a=0.5*((double)no-1.0);
	double b=1.0;
	double phi=1.0;
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);

	//Ya Variables//
	std::vector<double> mu(na);
	std::vector<double> Z; Z.reserve(na);
	std::vector<double> xaya(p);

	//Beta Variables//
	std::vector<double> B(p,0.0);
	std::vector<double> Bg; Bg.reserve(p);

	//Gamma Variables//
	std::vector<int> gamma(p,0);
	std::vector<double> Bols(p);
	std::vector<double> prob(p);
	std::vector<double> odds(p);
	bool gamma_diff=true;
	int p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);

	//Theta Variables//
	double theta=0.5;

	//Allocate Space for MCMC Draws//
	std::vector<double> B_mcmc(p*niter);
	std::vector<double> prob_mcmc(p*niter);
	std::vector<int>  gamma_mcmc(p*niter);
	std::vector<double> phi_mcmc(niter);
	std::copy(B.begin(),B.end(),B_mcmc.begin());
	std::copy(prob.begin(),prob.end(),prob_mcmc.begin());
	std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin());
	phi_mcmc[0]=phi;



	//Run Gibbs Sampler//
	for (int t = 1; t < niter; ++t)
	{
		//Form Submatrices If gamma[t]!=Gamma[t-1]//
		if(gamma_diff)
		{
			if(p_gamma!=0)
			{
				submatrices_collapsed(mu,xag,xogxog_Lamg,xogyo,lamg,Bg,gamma,xoyo,lam,xa,xoxo,p_gamma,b,yoyo,p,na);
			}else
			{
				a=0.5*((double)no-1.0);
				b=0.5*yoyo;
			}
		}


		//Draw Phi//
		phi=Rf_rgamma(a,(1.0/b)); //rgamma uses scale


		//Draw Ya//
		draw_collapsed_xaya(xaya,xa,xag,mu,phi,Z,xogxog_Lamg,na,p,p_gamma);


		//Draw Gamma//
		if(modelprior=="bernoulli")
		{
			bernoulli_probabilities(prob,odds,Bols,d,xoyo,xaya,priorprob,lam,phi);
		}else
		{
			betabinomial_probabilities(prob,odds,Bols,d,xoyo,xaya,theta,lam,phi);
		}
		draw_gamma(gamma,prob);


		//Draw Theta//
		p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);
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

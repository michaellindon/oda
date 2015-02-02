#include <Rcpp.h>
#include "oda.h"

// [[Rcpp::export]]
Rcpp::List col_normal_gibbs(Rcpp::NumericVector ryo, Rcpp::NumericMatrix rxo, Rcpp::NumericMatrix rxa, Rcpp::NumericMatrix rxoxo, Rcpp::NumericVector rd, Rcpp::NumericVector rlam, Rcpp::NumericVector rpriorprob, SEXP rburnin, SEXP rniter){

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
	std::vector<double> priorprob(rpriorprob.begin(),rpriorprob.end());
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);

	//Create Matrices//
	std::vector<double> xoyo(p);
	dgemv_( &transT, &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);
	std::vector<double> xogyo; xogyo.reserve(p);
	std::vector<double> xogxog_Lamg; xogxog_Lamg.reserve(p*p);
	std::vector<double> xag; xag.reserve(na*p);
	std::vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Phi Variables//
	double a=(double)0.5*(no-1);
	double b=1.0;
	double phi=1.0;
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);

	//Ya Variables//
	std::vector<double> mu(na);
	std::vector<double> Z; Z.reserve(na);
	std::vector<double> xaya(p);

	//Beta Variables//
	std::vector<double> Bg; Bg.reserve(p);

	//Gamma Variables//
	std::vector<int> gamma(p);
	std::vector<double> U(p);
	std::vector<double> Bols(p);
	std::vector<double> prob(p);
	std::vector<double> odds(p);
	std::vector<double> priorodds(p);
	std::vector<double> ldl(p);
	std::vector<double> dli(p);
	for(size_t i=0; i!=priorodds.size(); ++i){
		priorodds[i]=priorprob[i]/(1-priorprob[i]);
		ldl[i]=sqrt(lam[i]/(d[i]+lam[i]));
		dli[i]=1/(d[i]+lam[i]);
	}
	bool gamma_diff=true;
	int p_gamma;

	//Allocate Space for MCMC Draws//
	std::vector<double> prob_mcmc(p*niter);
	std::vector<int>  gamma_mcmc(p*niter);
	std::vector<double> phi_mcmc(niter);
	std::copy(prob.begin(),prob.end(),prob_mcmc.begin());
	std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin());
	phi_mcmc[0]=phi;



	//Begin Gibbs Sampling Algorithm//
	for (int t = 1; t < niter; ++t)
	{
		//If gamma[t]!=Gamma[t-1] Form Submatrices//
		if(gamma_diff)
		{
			p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);
			if(p_gamma!=0){
				submatrices_collapsed(mu,xag,xogxog_Lamg,xogyo,lamg,Bg,gamma,xoyo,lam,xa,xoxo,p_gamma,b,yoyo,p,na);
			}else{
				b=0.5*yoyo;
			}
		}

		//Draw Phi//
		phi=Rf_rgamma(a,(1.0/b)); //rgamma uses scale

		//Draw Ya//
		draw_collapsed_xaya(xaya,xa,xag,mu,phi,Z,xogxog_Lamg,na,p,p_gamma);

		//Draw Gamma//
		fixed_probabilities(prob, odds, Bols, d, xoyo, xaya, priorodds, ldl, dli, phi);
		draw_gamma(gamma,prob);

		//Store Draws//
		std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin()+p*t);
		std::copy(prob.begin(),prob.end(),prob_mcmc.begin()+p*t);
		phi_mcmc[t]=phi;

		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_mcmc,t,p);

	}


	//Return Data to R//
	return Rcpp::List::create(
			Rcpp::Named("phi_mcmc") = phi_mcmc,
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("gamma_mcmc") = gamma_mcmc
			) ;
}


//Notes: Lam is O(pxp) space complexity, terrible, just use lam O(p) instead.
//Better to pass xoxo and d instead of xoxo & xaxa and computing d inside the code, as this requires O(pxp) instead of O(p) [xa is needed anyway] and xa.t()*xa takes forever
//dtritri otherwise need to do dpotrs and then dpotrf again

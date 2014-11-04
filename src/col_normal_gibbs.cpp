#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List col_normal_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter){

	//Define Variables//
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	double a=(double)0.5*(no-1);
	double b;
	double phi=1;
	Mat<double> xag;
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> xoxog(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> Lamg(p,p);
	Mat<double> E(na,na);
	Mat<double> L(na,na);
	Mat<double> xog;
	Mat<double> ya_mcmc(na,niter,fill::zeros);
	Mat<double> prob_mcmc(p,niter,fill::zeros);
	Mat<uword>  gamma_mcmc(p,niter,fill::ones);
	Col<double> phi_mcmc(niter,fill::ones);
	Col<double> mu(na);
	Col<double> ya(na,fill::zeros);
	Col<double> Z(na);
	Col<double> U(p);
	Col<double> d(p);
	Col<double> Bols(p,fill::zeros);
	Col<double> xoyo(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p);
	Col<double> ldl(p);
	Col<double> dli(p);
	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p,fill::ones);


	//Copy RData Into Matrix Classes//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	yo-=mean(yo);


	//Create Matrices//
	xoyo=xo.t()*yo;
	xoxo=xo.t()*xo;
	xaxa=xa.t()*xa;
	D=xaxa+xoxo;
	d=D.diag();
	Lam=diagmat(lam);
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);
	double yoyo=dot(yo,yo);

	ya_mcmc.col(0)=ya;
	phi_mcmc(0)=phi;
	gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;


	//Run Gibbs Sampler//
	for (int t = 1; t < niter; ++t)
	{
		//Form Submatrices
		inc_indices=find(gamma);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);
		xoxog=xoxo.submat(inc_indices,inc_indices);

		//Draw Phi//
		b=0.5*(yoyo-dot(xog.t()*yo, solve(xoxog+Lamg,xog.t()*yo)));
		phi=R::rgamma(a,(1/b)); //rgamma uses scale

		//Draw Ya//
		mu=xag*solve(xoxog+Lamg,xog.t()*yo);
		E=xag*solve(xoxog+Lamg,xag.t());
		for(int i=0; i < na; ++i) E(i,i)+=1;
		E=E/phi;
		L=chol(E);
		for (int i = 0; i < na; ++i) Z(i)=R::rnorm(0,1);
		ya=mu+L.t()*Z;

		Bols=(1/d)%(xoyo+xa.t()*ya);
		odds=priorodds%ldl%trunc_exp(0.5*phi*dli%d%d%Bols%Bols);
		prob=odds/(1+odds);

		//Draw Gamma//
		for (int i = 0; i < p; ++i) U(i)=R::runif(0,1);
		for (int i = 0; i < p; ++i)
		{
			if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN
			if(U(i)<prob(i)){
				gamma(i)=1;
			}else{
				gamma(i)=0;
			}
		}

		//Store Draws//
		gamma_mcmc.col(t)=gamma;
		prob_mcmc.col(t)=prob;
		ya_mcmc.col(t)=ya;
		phi_mcmc(t)=phi;
	}


	return Rcpp::List::create(
			Rcpp::Named("phi_mcmc") = phi_mcmc,
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("prob") = mean(prob_mcmc.cols(burnin-1,niter-1),1),
			Rcpp::Named("gamma_mcmc") = gamma_mcmc,
			Rcpp::Named("ya_mcmc") = ya_mcmc
			) ;
}

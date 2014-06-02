#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List col_mixture_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter, SEXP ralpha){

	//Define Variables//
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	double a=(double)0.5*(no-1);
	double b;
	double phi;
	double alpha=Rcpp::as<double >(ralpha);
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
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> P1(no,no);
	Mat<double> Px(no,no);
	Mat<double> ya_mcmc(na,niter,fill::zeros);
	Mat<double> prob_mcmc(p,niter,fill::zeros);
	Mat<double> B_mcmc(p,niter,fill::zeros);
	Mat<double> lam_mcmc(p,niter,fill::zeros);
	Mat<uword>  gamma_mcmc(p,niter,fill::ones);
	Col<double> phi_mcmc(niter,fill::ones);
	Col<double> one(no,fill::ones);
	Col<double> mu(na);
	Col<double> ya(na);
	Col<double> Z(na);
	Col<double> d(p);
	Col<double> Bols(p,fill::zeros);
	Col<double> B(p);
	Col<double> xoyo(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p);
	Col<double> ldl(p);
	Col<double> dli(p);
	Col<uword> gamma(p,fill::ones);
	Col<uword> inc_indices(p,fill::ones);


	//Copy RData Into Matrix Classes//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);


	//Create Matrices//
	xoyo=xo.t()*yo;
	xoxo=xo.t()*xo;
	xaxa=xa.t()*xa;
	D=xaxa+xoxo;
	d=D.diag();
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);


	//Initialize Parameters at MLE//
	P1=one*(one.t()*one).i()*one.t();
//	Px=xo*(xoxo).i()*xo.t();
	//phi=(no-1)/dot(yo,((Ino-P1-Px)*yo));
  phi=1;
//	Bols=(xoxo).i()*xo.t()*yo;
	ya=xa*Bols;

	//Run Gibbs Sampler//
	ya_mcmc.col(0)=ya;
	phi_mcmc(0)=phi;
	gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;
	for (int t = 1; t < niter; ++t)
	{
		//Form Submatrices
		inc_indices=find(gamma);
		Lam=diagmat(lam);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);
		xoxog=xoxo.submat(inc_indices,inc_indices);

		//Draw Phi//
		b=0.5*dot(yo,(Ino-P1-xog*(xoxog+Lamg).i()*xog.t())*yo);
		phi=R::rgamma(a,(1/b)); //rgamma uses scale

		//Draw Ya//
		mu=xag*(xoxog+Lamg).i()*xog.t()*yo;
		E=Ina+xag*(xoxog+Lamg).i()*xag.t();
		E=E/phi;
		L=chol(E);
		for (int i = 0; i < na; ++i) Z(i)=R::rnorm(0,1);
		ya=mu+L.t()*Z;

		Bols=(1/d)%(xoyo+xa.t()*ya);
		odds=priorodds%ldl%trunc_exp(0.5*phi*dli%d%d%Bols%Bols);
		prob=odds/(1+odds);
		for (int i = 0; i < p; ++i)
		{
			//Draw Gamma, Beta, lambda//
			if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN
			if(R::runif(0,1)<prob(i)){
				gamma(i)=1;
				Z(i)=R::rnorm(0,1);
				B(i)=dli(i)*d(i)*Bols(i)+sqrt(dli(i)/phi)*Z(i);
				lam(i)=R::rgamma(0.5*(alpha+1),2/(alpha+phi*B(i)*B(i))); //rgamma uses scale
			}else{
				gamma(i)=0;
				B(i)=0;
				lam(i)=R::rgamma(0.5*alpha,2/alpha); //rgamma uses scale
			}
		}
		ldl=sqrt(lam/(d+lam));
		dli=1/(d+lam);

		//Store Draws//
		lam_mcmc.col(t)=lam;
		B_mcmc.col(t)=B;
		gamma_mcmc.col(t)=gamma;
		prob_mcmc.col(t)=prob;
		ya_mcmc.col(t)=ya;
		phi_mcmc(t)=phi;
	}


	return Rcpp::List::create(
			Rcpp::Named("B_mcmc") = B_mcmc,
			Rcpp::Named("B") = B,
			Rcpp::Named("lam_mcmc") = lam_mcmc,
			Rcpp::Named("lam") = lam,
			Rcpp::Named("phi_mcmc") = phi_mcmc,
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("prob") = mean(prob_mcmc.cols(burnin-1,niter-1),1),
			Rcpp::Named("gamma_mcmc") = gamma_mcmc,
			Rcpp::Named("ya_mcmc") = ya_mcmc
			) ;
}

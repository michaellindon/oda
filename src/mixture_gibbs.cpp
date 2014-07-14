#include <chrono>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List mixture_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter, SEXP ralpha){

	//Define Variables//
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
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
	Mat<double> L(na,na);
	Mat<double> xog;
	Mat<double> P1(no,no);
	Mat<double> Px(no,no);
	Mat<double> ya_mcmc(na,niter,fill::zeros);
	Mat<double> prob_mcmc(p,niter,fill::zeros);
	Mat<double> B_mcmc(p,niter,fill::zeros);
	Mat<double> lam_mcmc(p,niter,fill::zeros);
	Mat<uword>  gamma_mcmc(p,niter,fill::ones);
	Col<double> phi_mcmc(niter,fill::ones);
	Col<double> one(no,fill::ones);
	Col<double> ya(na);
	Col<double> Z(na);
	Col<double> d(p);
	Col<double> yoc(no);
	Col<double> Bols(p);
	Col<double> B(p);
	Col<double> Bg;
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
	P1=one*(one.t()*one).i()*one.t();
	Px=xo*(xoxo).i()*xo.t();
	yoc=yo-mean(yo);


	//Initialize Parameters at MLE//
	phi=1;
	Bols=(xoxo).i()*xo.t()*yo;
	ya=xa*Bols;
	B=Bols;

	//Run Gibbs Sampler//
	ya_mcmc.col(0)=ya;
	phi_mcmc(0)=phi;
	gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;
	lam_mcmc.col(0)=lam;
	B_mcmc.col(0)=B;
	auto start = std::chrono::steady_clock::now();
	for (int t = 1; t < niter; ++t)
	{

		//Form Submatrices
		Lam=diagmat(lam);
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);
		xoxog=xoxo.submat(inc_indices,inc_indices);

		//Draw Phi//
		b=0.5*dot(yoc-xog*Bg,yoc-xog*Bg)+0.5*dot(Bg,Lamg*Bg);
		phi=R::rgamma((double)0.5*(no+sum(gamma)-1),(1/b)); //rgamma uses scale

		//Draw Ya//
		for (int i = 0; i < na; ++i) Z(i)=R::rnorm(0,1);
		ya=xag*Bg+sqrt(1/phi)*Z;


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
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed=end-start;

	std::cout <<  elapsed.count() << " sec - Total Runtime" << std::endl;
	std::cout <<  elapsed.count()/niter << " sec - Per Iteration (avg)" << std::endl;


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

#include <chrono>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List normal_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd,NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter){

	//Define Variables//
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();

	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> xoyo=xo.t()*yo;
	Col<double> xaya(p);
	Mat<double> xat=xa.t();
	Mat<double> D=diagmat(d);
	Mat<double> Lam=diagmat(lam);
	Col<double> priorodds=priorprob/(1-priorprob);
	Col<double> odds=priorodds;
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);
	Col<double> prob(p,fill::ones);	prob*=0.5;
	Col<double> ya(na);
	Col<double> Z(na);
	Col<double> Bols(p);
	Col<double> B(p,fill::zeros);
	double b;
	double phi=1;

	//Create Sub Matrices//
	Col<uword> gamma(p,fill::ones);
	Col<uword> inc_indices(p,fill::ones);
	Col<double> Bg;
	Mat<double> xog(no,p);
	Mat<double> xag(na,p);
	Mat<double> Lamg(p,p);

	//Create Trace Matrices//
	Mat<double> ya_mcmc(na,niter,fill::zeros);
	Mat<double> prob_mcmc(p,niter,fill::zeros);
	Mat<double> B_mcmc(p,niter,fill::zeros);
	Mat<uword>  gamma_mcmc(p,niter,fill::ones);
	Col<double> phi_mcmc(niter,fill::ones);

	//Run Gibbs Sampler//
	ya_mcmc.col(0)=ya;
	phi_mcmc(0)=phi;
	gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;
	B_mcmc.col(0)=B;
	auto start = std::chrono::steady_clock::now();
	for (int t = 1; t < niter; ++t)
	{

		//Form Submatrices//
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);

		//Draw Phi//
		b=0.5*dot(yo-xog*Bg,yo-xog*Bg)+0.5*dot(Bg,Lamg*Bg);
		phi=R::rgamma((double)0.5*(no+sum(gamma)-1),(1/b)); //rgamma uses scale

		//Draw Ya//
		for (int i = 0; i < na; ++i) Z(i)=R::rnorm(0,1);
		ya=xag*Bg+sqrt(1/phi)*Z;
		xaya=xat*ya;

		for(int i=0; i<p; ++i){
			Bols(i)=(1/d(i))*(xoyo(i)+xaya(i));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
			prob(i)=odds(i)/(1+odds(i));
		}

		for(int i=0; i<p; ++i)  if(prob(i)!=prob(i)) prob(i)=1;   //Catch NaN


		//Draw Gamma, Beta//
		for (int i = 0; i < p; ++i)
		{
			if(R::runif(0,1)<prob(i)){
				gamma(i)=1;
				Z(i)=R::rnorm(0,1);
				B(i)=dli(i)*d(i)*Bols(i)+sqrt(dli(i)/phi)*Z(i);
			}else{
				gamma(i)=0;
				B(i)=0;
			}
		}

		//Store Draws//
		B_mcmc.col(t)=B;
		gamma_mcmc.col(t)=gamma;
		prob_mcmc.col(t)=prob;
		ya_mcmc.col(t)=ya;
		phi_mcmc(t)=phi;
	}

	//Report Runtime//
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed=end-start;
	std::cout <<  elapsed.count() << " sec - Total Runtime" << std::endl;
	std::cout <<  elapsed.count()/niter << " sec - Per Iteration (avg)" << std::endl;


	return Rcpp::List::create(
			Rcpp::Named("B_mcmc") = B_mcmc,
			Rcpp::Named("B") = mean(B_mcmc.cols(burnin-1,niter-1),1),
			Rcpp::Named("phi_mcmc") = phi_mcmc,
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("prob") = mean(prob_mcmc.cols(burnin-1,niter-1),1),
			Rcpp::Named("gamma_mcmc") = gamma_mcmc,
			Rcpp::Named("ya_mcmc") = ya_mcmc
			) ;
}

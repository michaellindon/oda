#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List probit_normal_gibbs(NumericVector rbit, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd,NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter){

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
	arma::colvec bit(rbit.begin(), rbit.size(), false);

	//Pre-Processing//
	Mat<double> xot=xo.t();
	Mat<double> xat=xa.t();
	Col<double> xoyo(p);
	Col<double> xaya(p);
	Mat<double> D=diagmat(d);
	Mat<double> Lam=diagmat(lam);
	Col<double> logpriorodds=log(priorprob/(1-priorprob));
	Col<double> odds=logpriorodds;
	Col<double> logldl=0.5*log(lam/(d+lam));
	Col<double> dli=1/(d+lam);
	Col<double> prob(p,fill::ones);	prob*=0.5;
	Col<double> yo(no,fill::zeros);
	Col<double> ya(na);
	Col<double> U(no);
	Col<double> mu_yo(no,fill::zeros);
	Col<double> Z(na);
	Col<double> Bols(p);
	Col<double> B(p,fill::zeros);
	double phi=1;
	double Fa;

	//Create Sub Matrices//
	Col<uword> gamma(p,fill::ones);
	Col<uword> inc_indices(p,fill::ones);
	Col<double> Bg;
	Mat<double> xog(no,p);
	Mat<double> xag(na,p);
	Mat<double> Lamg(p,p);

	//Create Trace Matrices//
	//Mat<double> yo_mcmc(no,niter,fill::zeros);
	//Mat<double> ya_mcmc(na,niter,fill::zeros);
	Mat<double> prob_mcmc(p,niter,fill::zeros);
	Mat<double> B_mcmc(p,niter,fill::zeros);
	//Mat<uword>  gamma_mcmc(p,niter,fill::ones);

	//Run Gibbs Sampler//
	//ya_mcmc.col(0)=ya;
	//gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;
	B_mcmc.col(0)=B;
	for (int t = 1; t < niter; ++t)
	{

		//Form Submatrices//
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);

		//Draw Yo//
		for (int i = 0; i < no; ++i) U(i)=R::runif(0,1);
		mu_yo=xog*Bg;
		for (int i = 0; i < no; ++i)
		{
			if(bit(i)==1){
				Fa=R::pnorm(-mu_yo(i),0,1,1,0);
				yo(i)=mu_yo(i)+R::qnorm(Fa+U(i)*(1-Fa),0,1,1,0);
			}else{
				yo(i)=mu_yo(i)+R::qnorm(U(i)*R::pnorm(-mu_yo(i),0,1,1,0),0,1,1,0);
			}
		}	
		xoyo=xot*yo;

		//Phi=1//

		//Draw Ya//
		for (int i = 0; i < na; ++i) Z(i)=R::rnorm(0,1);
		ya=xag*Bg+sqrt(1/phi)*Z;
		xaya=xat*ya;

		for(int i=0; i<p; ++i){
			Bols(i)=(1/d(i))*(xoyo(i)+xaya(i));
			odds(i)=logpriorodds(i)+logldl(i)+(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
      odds(i)=exp(odds(i));
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
	//	gamma_mcmc.col(t)=gamma;
		prob_mcmc.col(t)=prob;
	//	yo_mcmc.col(t)=yo;
	//	ya_mcmc.col(t)=ya;
	}



	return Rcpp::List::create(
			Rcpp::Named("B_mcmc") = B_mcmc,
			Rcpp::Named("B") = mean(B_mcmc.cols(burnin-1,niter-1),1),
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("prob") = mean(prob_mcmc.cols(burnin-1,niter-1),1)
	//		Rcpp::Named("gamma_mcmc") = gamma_mcmc,
	//		Rcpp::Named("yo_mcmc") = yo_mcmc,
	//		Rcpp::Named("ya_mcmc") = ya_mcmc
			) ;
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int t=1;
	double a=(no-1)/2;
	double b;
	double delta;
	Mat<double> xag;
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> xoxog(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> Lamg(p,p);
	Mat<double> E(na,na);
	Mat<double> xog;
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> P1(no,no);
	Mat<double> prob_trace(p,1000,fill::zeros);
	Mat<uword>  gamma_trace(p,1000,fill::ones);
	Col<double> one(no,fill::ones);
	Col<double> mu(na);
	Col<double> ya(na);
	Col<double> d(p);
	Col<double> Bols(p);
	Col<double> xoyo(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p,fill::ones);
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

	//Create Matrices//
	xoyo=xo.t()*yo;
	xoxo=xo.t()*xo;
	xaxa=xa.t()*xa;
	D=xaxa+xoxo;
	d=D.diag();
	P1=one*(one.t()*one).i()*one.t();
	prob=priorprob;
	Lam=diagmat(lam);
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);

	//Randomize Initial Gammas//
	for (int i = 0; i < p; ++i) if(R::runif(0,1)>0.5) gamma(i)=1;

	//Run Gibbs Sampler//
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	do{
		//Form Submatrices
		inc_indices=find(gamma);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);
		xoxog=xoxo.submat(inc_indices,inc_indices);

		//E-Step//
		//Phi//
		b=0.5*dot(yo,(Ino-P1-xog*(xoxog+Lamg).i()*xog.t())*yo);

		//Ya//
		mu=xag*(xoxog+Lamg).i()*xog.t()*yo;
		E=Ina+xag*(xoxog+Lamg).i()*xag.t();

		//Gamma//
		for (int i = 0; i < p; ++i)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+dot(xa.col(i),mu));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*(a/b)*dli(i)*d(i)*d(i)*Bols(i)*Bols(i)+0.5*dli(i)*dot(xa.col(i),E*xa.col(i)));
			prob(i)=odds(i)/(1+odds(i));
			//if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN

			//M-Step//
			//Choose Median Probability Model//
			if(0.5<prob(i)){
				gamma(i)=1;
			}else{
				gamma(i)=0;
			}
		}

		//Store Values//
		gamma_trace.col(t)=gamma;
		prob_trace.col(t)=prob;

		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		t=t+1;
	} while(delta>0.00001);

	gamma_trace.resize(p,t);
	prob_trace.resize(p,t);

	return Rcpp::List::create(
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("gamma") = gamma,
			Rcpp::Named("gamma_trace") = gamma_trace
			) ;

}

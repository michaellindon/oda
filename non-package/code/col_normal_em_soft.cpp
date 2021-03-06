#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List col_normal_em_soft(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int t=1;
	double a=(double)0.5*(no-3);
	double b;
	double delta;
	double phi=1;
	Mat<double> Aooa(no,no);
	Mat<double> A(no+na,no+na);
	Mat<double> Aoo(no,no);
	Mat<double> Aaa(na,na);
	Mat<double> Aoa(no,na);
	Mat<double> Aao(na,no);
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> xcxcLami(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> E(na,na);
	Mat<double> xc(no+na,p);
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> Inc=eye(no+na,no+na);
	Mat<double> P=eye(p,p);
	Mat<double> P1(no,no);
	Mat<double> ya_trace(p,1000,fill::ones);
	Mat<double> prob_trace(p,1000,fill::zeros);
	Col<double> phi_trace(1000,fill::ones);
	Col<double> one(no,fill::ones);
	Col<double> mu(na);
	Col<double> ya(na,fill::zeros);
	Col<double> d(p);
	Col<double> Bols(p);
	Col<double> xoyo(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p,fill::ones);
	Col<double> ldl(p);
	Col<double> dli(p);

	//Copy RData Into Matrix Classes//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);

	//Create Matrices//
	xc=join_cols(xo,xa);
	xoyo=xo.t()*yo;
	xoxo=xo.t()*xo;
	xaxa=xa.t()*xa;
	D=xaxa+xoxo;
	d=D.diag();
	xcxcLami=diagmat(1/(d+lam));
	Lam=diagmat(lam);
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);
	P1=one*(one.t()*one).i()*one.t();

	//Randomize Initial Probabilities//
	for (int i = 0; i < p; ++i) prob(i)=R::runif(0,1);
	P.diag()=prob;

	//Run EM//
	prob_trace.col(0)=prob;
	do{
		//Phi Maximization Step//
		A=Inc-xc*P*xcxcLami*xc.t();
		Aoo=A.submat(0,0,no-1,no-1)-P1;
		Aaa=A.submat(no,no,no+na-1,no+na-1);
		Aoa=A.submat(0,no,no-1,no+na-1);
		Aao=A.submat(no,0,no+na-1,no-1);
		Aooa=Aoo-Aoa*Aaa.i()*Aao;
		b=0.5*(dot(yo,Aooa*yo));
		phi=((double)a)/b;

		//Ya Maximization Step//
		ya=-solve(Aaa,Aao*yo);

		//Probability Expectation Step//
		for (int i = 0; i < p; i++)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+dot(xa.col(i),ya));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*(d(i)*d(i)*Bols(i)*Bols(i)));
			prob(i)=odds(i)/(1+odds(i));
		}
		P.diag()=prob;

		//Store Values//
		prob_trace.col(t)=prob;
		ya_trace.col(t-1)=ya;
		phi_trace(t-1)=phi;

		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		t=t+1;
	} while(delta>0.001);

	phi_trace.resize(t-1);
	ya_trace.resize(p,t-1);
	prob_trace.resize(p,t);

	return Rcpp::List::create(
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("ya") = ya,
			Rcpp::Named("ya_trace") = ya_trace,
			Rcpp::Named("phi") = phi,
			Rcpp::Named("phi_trace") = phi_trace
			) ;

}
